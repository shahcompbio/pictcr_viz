import React, { useState, useRef, useEffect } from "react";

import * as d3 from "d3";
import { useD3, DownloadIcon } from "@shahlab/planetarium";
import _ from "lodash";
import { drag as d3Drag } from "d3-drag";
import { ease as d3Ease } from "d3-ease";

import makeStyles from "@mui/styles/makeStyles";

import { Layout } from "@shahlab/planetarium";
import Tooltip from "@mui/material/Tooltip";
import Typography from "@mui/material/Typography";
import Table from "@mui/material/Table";
import TableBody from "@mui/material/TableBody";
import TableCell from "@mui/material/TableCell";
import TableRow from "@mui/material/TableRow";

import { INFO } from "../config";

import { jsPDF } from "jspdf";
import canvg from "canvg";

const TOP_NUM = 3;

const moveScrollerToSpot = (xPos) =>
  d3
    .select("#scroller")
    .transition()
    .ease(d3.easeQuad)
    .duration(100)
    .attr("x", xPos);

const ScrollBar = ({ inputRef, inputRefMapping, scrollBarWidth, height }) => {
  const [xScale, setXScale] = useState(null);

  const getBandWidth = (d, scrollBarWidthMapping, fullBodyWidth) =>
    (scrollBarWidthMapping[d["refName"]] * scrollBarWidth) /
    (fullBodyWidth - 342);

  const getPlotWidths = (inputRefMapping, inputRef) => {
    const fullBodyWidth = d3.select("body").node().scrollWidth;
    var scrollBarWidthMapping = {};
    const plotWidths = inputRefMapping
      .map((d) => d["refName"])
      .map((d, i) => {
        const plotWidth = d3
          .select(inputRef.current[d])
          .node()
          .getBoundingClientRect().width;

        scrollBarWidthMapping[d] = plotWidth;
        return plotWidth;
      });
    var scrollBarWidthXMapping = {};
    inputRefMapping
      .map((d) => d["refName"])
      .map((d, i) => {
        const boundingBox = d3
          .select(inputRef.current[d])
          .node()
          .getBoundingClientRect();

        scrollBarWidthXMapping[d] = boundingBox;
        return boundingBox;
      });

    const barWidths = inputRefMapping.reduce(
      (final, d) => {
        const bandWidth = getBandWidth(d, scrollBarWidthMapping, fullBodyWidth);
        const xStart = final["total"];
        final["total"] = xStart + bandWidth;
        final[d["refName"]] = { xStart: xStart, bandWidth: bandWidth };
        return final;
      },
      { total: 0 }
    );
    return {
      barWidths: barWidths,
      plotWidths: plotWidths,
      plotWidthXMapping: scrollBarWidthXMapping,
      scrollBarWidthMapping: scrollBarWidthMapping,
    };
  };
  useEffect(() => {
    var ticking = false;
    var lastScrollLeft = 0;

    var calcLeft = function (rect, win) {
      return Math.max(rect.width + rect.x, 0) / rect.width;
    };
    var calcRight = function (rect, win) {
      return Math.max(win.innerWidth - rect.x, 0) / rect.width;
    };

    const onScroll = () => {
      if (!ticking) {
        window.requestAnimationFrame(function () {
          var documentScrollLeft = window.scrollX;
          if (lastScrollLeft != documentScrollLeft) {
            var plotIndex = 0;
            const {
              plotWidths,
              barWidths,
              plotWidthXMapping,
              scrollBarWidthMapping,
            } = getPlotWidths(inputRefMapping, inputRef);
            //in between lastScrollLeft and window.innerWidth
            const windowWidth = window.innerWidth - lastScrollLeft;

            const objectsInView = Object.keys(plotWidthXMapping)
              .filter((d) => d !== "total")
              .sort(function (a, b) {
                return inputRefMapping.indexOf(a) - inputRefMapping.indexOf(b);
              })
              .reduce((final, d) => {
                const rect = plotWidthXMapping[d];
                const percentage =
                  rect.x < 0
                    ? calcLeft(rect, window)
                    : rect.x + rect.width > window.innerWidth
                    ? calcRight(rect, window)
                    : 1;
                const location =
                  rect.x < 0
                    ? "left"
                    : rect.x + rect.height > window.innerWidth
                    ? "right"
                    : "center";
                final[d] = { percentage: percentage, location: location };
                return final;
              }, {});
            const p = plotWidths.reduce((final, curr) => {
              final = [final, { curr: plotWidths[curr["startX"]] }];
              return final;
            }, []);

            const sortedxStarts = Object.keys(barWidths)
              .filter((d) => d !== "total")
              .map((d) => ({
                xStart: barWidths[d]["xStart"],
                bandWidthStart: barWidths[d]["bandWidth"],
              }))
              .sort((a, b) => a["xStart"] - b["xStart"]);

            const sortedBarWidths = Object.keys(barWidths)
              .filter((d) => d !== "total")
              .sort(function (a, b) {
                return inputRefMapping.indexOf(a) - inputRefMapping.indexOf(b);
              });
            inputRefMapping.map((d) => {
              if (objectsInView[d["refName"]]["percentage"] >= 0.75) {
                d3.select("#text-" + d["refName"])
                  .attr("font-weight", "bold")
                  .attr("font-size", 20);
              } else {
                d3.select("#text-" + d["refName"])
                  .attr("font-weight", "normal")
                  .attr("font-size", 16);
              }
            });
            const finalWidth = inputRefMapping
              .filter((d) => d !== "total")
              .reduce(
                (final, d) => {
                  // /  console.log(d);
                  const name = d["refName"];
                  if (objectsInView[d["refName"]]["percentage"] !== 0) {
                    //if this is the first visible plot
                    if (!final.hasOwnProperty("start")) {
                      final["start"] = d;
                      const element = barWidths;

                      moveScrollerToSpot(
                        objectsInView[name]["location"] === "center"
                          ? barWidths[name]["xStart"]
                          : barWidths[name]["xStart"] +
                              (1 - objectsInView[name]["percentage"]) *
                                barWidths[name]["bandWidth"]
                      );
                      final["width"] =
                        objectsInView[name]["location"] === "center"
                          ? barWidths[name]["bandWidth"]
                          : objectsInView[name]["percentage"] *
                            barWidths[name]["bandWidth"];
                    } else {
                      //change the scroller width
                      final["width"] =
                        final["width"] +
                        (objectsInView[name]["location"] === "center"
                          ? barWidths[name]["bandWidth"]
                          : objectsInView[name]["percentage"] *
                            barWidths[name]["bandWidth"]);
                    }
                  }
                  return final;
                },
                { width: 0 }
              )["width"];
            d3.select("#scroller").attr("width", finalWidth);

            lastScrollLeft = documentScrollLeft;
          }

          ticking = false;
        });
        ticking = true;
      }
    };
    window.addEventListener("scroll", onScroll);

    return () => {
      window.removeEventListener("scroll", onScroll);
    };
  }, []);

  const drawScroll = (svg, inputRefMapping) => {
    const { barWidths, plotWidths, scrollBarWidthMapping, getBandWidth } =
      getPlotWidths(inputRefMapping, inputRef);

    const xScale = d3
      .scaleOrdinal()
      .domain([...plotWidths])
      .range([0, scrollBarWidth]);

    const maxXPos = plotWidths.reduce((final, d) => final + d, 0);

    const enter = svg.selectAll("g").data(inputRefMapping).enter();
    const gElement = enter.append("g").attr("class", "scrollWrapper");
    const getScrollerWidth = () => {
      var plotIndex = 0;
      const lastScrollLeft = window.scrollX + window.screen.availWidth;

      const filters = d3.select("#filters-grid").node().getBoundingClientRect();

      [filters.width, ...plotWidths].reduce((total, currWidth, i) => {
        if (total !== null) {
          if (total + currWidth > lastScrollLeft) {
            total = null;
          } else {
            plotIndex = i - 1;
            total = total + currWidth;
          }
        }
        return total;
      }, 0);
      const lastFullyVisibleItem = inputRefMapping[plotIndex];

      const sortedBarWidths = Object.keys(barWidths)
        .filter((d) => d !== "total")
        .map((d) => ({
          xStart: barWidths[d]["xStart"],
          bandWidth: barWidths[d]["bandWidth"],
        }))
        .sort((a, b) => a["xStart"] - b["xStart"]);

      return sortedBarWidths
        .filter((d, i) => i <= plotIndex)
        .reduce((f, d) => f + d["bandWidth"], 0);
    };
    gElement
      .append("rect")
      .attr("x", (d, i) => barWidths[d["refName"]]["xStart"])
      .attr("y", 11)
      .attr("width", (d) => barWidths[d["refName"]]["bandWidth"])
      .attr("height", 22)
      .attr("fill", "#f0eded");

    gElement
      .append("text")
      .attr(
        "x",
        (d, i) =>
          barWidths[d["refName"]]["xStart"] +
          barWidths[d["refName"]]["bandWidth"] / 2
      )
      .attr("y", 15)
      .attr("dy", ".35em")
      .attr("text-anchor", "middle")
      .attr("dominant-baseline", "central")
      .attr("fill", "black")
      //.attr("font-weight", "light")
      .attr("id", (d) => "text-" + d["refName"])
      .text((d) => d["plotName"]);

    const scrollToByRef = (refName) =>
      inputRef.current[refName].scrollIntoView({
        behavior: "smooth",
        block: "center",
        inline: "center",
      });
    const snapToGrid = (event) => {
      const currXPos = d3.select("#scroller").node().getBoundingClientRect().x;
      const sortedBarWidths = Object.keys(barWidths)
        .filter((d) => d !== "total")
        .map((d) => ({
          xStart: barWidths[d]["xStart"],
          bandWidth: barWidths[d]["bandWidth"],
        }))
        .sort((a, b) => a["xStart"] - b["xStart"]);
      //sortedBarWidths.map;
    };
    const setScrollerPosition = (event, d) => {
      moveScrollerToSpot(event.x);
    };
    const setScrollWidth = (event) => {
      snapToGrid(event);
      getScrollerWidth();
    };
    const scroller = svg.append("rect").attr("id", "scroller");

    const dragger = d3Drag()
      .on("start", setScrollerPosition)
      .on("drag", setScrollerPosition)
      .on("end", setScrollWidth);

    scroller
      .attr("x", 0)
      .attr("y", 2.5)
      //.attr("fill", "#8686db")
      .attr("height", 40)
      .attr("width", getScrollerWidth())
      //  .attr("width", (d) => getBandWidth(d))
      .attr("rx", 5)
      .attr("opacity", 0.5);

    d3.select("#scroller").call(dragger);
    d3.selectAll(".scrollWrapper")
      .selectAll("rect")
      .on("click", setScrollerPosition);
  };

  const ref = useD3(
    (svg) => {
      drawScroll(svg, inputRefMapping);
    },
    scrollBarWidth,
    height,
    [inputRefMapping]
  );
  return (
    <div
      id="scrollbar-wrapper"
      style={{
        position: "fixed",
        width: "100vw",
        height: "10px",
        marginTop: 50,
        background: "white",
      }}
    >
      Scroll To:
      <svg ref={ref} style={{ display: "block", margin: "auto" }} />
    </div>
  );
};

export default ScrollBar;
