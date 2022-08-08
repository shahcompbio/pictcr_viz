import React, { useState, useRef, Fragment } from "react";

import * as d3 from "d3";
import { useD3, DownloadIcon } from "@shahlab/planetarium";
import _ from "lodash";
import { rgb } from "d3-color";

import makeStyles from "@mui/styles/makeStyles";
import InfoIcon from "@mui/icons-material/Info";
import { Layout } from "@shahlab/planetarium";
import Tooltip from "@mui/material/Tooltip";
import Grid from "@mui/material/Grid";
import Typography from "@mui/material/Typography";
import { StyledTitle } from "./util/util";
import Title from "./Title.js";

import { INFO } from "../config";

import { jsPDF } from "jspdf";
import canvg from "canvg";

const b = {
  w: 100,
  h: 30,
  s: 5,
  t: 10,
};

const ANGLE = 0.01;
const useStyles = makeStyles((theme) => ({
  tooltip: {
    minWidth: 220,
  },
  tooltipWrapper: {
    position: "relative",
    width: "10px",
    height: "10px",
  },
}));
const parentIsTopTen = (d) => d.parent.data.name === "Top Ten";

const Sunburst = ({
  data,
  colors,
  width,
  height,
  subsetParam,
  type,
  otherSubsetParam,
  cloneColorScale,
  cloneHEXColorScale,
  selectedCloneColor = null,
}) => {
  console.log(data);
  const [isTooltipOpen, setIsTooltipOpen] = useState(false);
  const tooltipRef = useRef(null);
  const [hoverItem, setHoverItem] = useState(null);
  const [hoverParent, setRootParent] = useState(null);

  const classes = useStyles();

  const totalCount = data.length;

  const allSubsets = _.groupBy(data, (datum) => datum[subsetParam]);

  const subsetCounts = Object.keys(allSubsets).sort(
    (a, b) => allSubsets[b].length - allSubsets[a].length
  );
  const topTenSubsetCounts = subsetCounts.slice(0, 10);
  const topTenCount = topTenSubsetCounts.reduce(
    (final, d) => (final = final + allSubsets[d].length),
    0
  );

  const resetSubsetcounts = subsetCounts.slice(10, subsetCounts.length);

  const restCount = resetSubsetcounts
    .filter((d) => allSubsets[d].length !== 1)
    .reduce((final, d) => (final = final + allSubsets[d].length), 0);

  const singletonCounts = resetSubsetcounts
    .map((d) => ({
      name: d,
      children: allSubsets[d],
    }))
    .filter((d) => d.children.length === 1);

  const restCounts = resetSubsetcounts
    .map((d) => ({
      name: d,
      children: allSubsets[d],
    }))
    .filter((d) => d.children.length !== 1);

  const colorScale = d3
    .scaleOrdinal()
    .domain([topTenSubsetCounts, "Other", "Top Ten"])
    .range(
      colors.slice(0, Math.min(topTenSubsetCounts.length, colors.length + 2))
    );

  const topTen = topTenSubsetCounts.map((d) => ({
    name: d,
    children: allSubsets[d],
  }));

  const topTenObj = {
    name: "Top Ten",
    children: [...topTen],
    cellCount: topTenCount,
  };
  const singleObj = {
    name: "Singleton",
    children: [...singletonCounts],
    cellCount: singletonCounts.length,
  };
  const otherObj = {
    name: "Other",
    children: [...restCounts],
    cellCount: restCounts.reduce(
      (final, curr) => final + curr.children.length,
      0
    ),
  };
  const hierarchy = {
    name: "flair",
    children: [topTenObj, otherObj, singleObj],
  };

  const format = d3.format(",d");
  const radius = width / 9;
  const arc = d3
    .arc()
    .startAngle((d) => d.x0)
    .endAngle((d) => d.x1)
    .padAngle((d) => Math.min((d.x1 - d.x0) / 2, 0.005))
    .padRadius(radius * 1.5)
    .innerRadius((d) => d.y0 * radius)
    .outerRadius((d) => Math.max(d.y0 * radius, d.y1 * radius - 1));

  const partition = (data) => {
    const root = d3
      .hierarchy(data)
      .sum((d) => 1)
      .sort((a, b) => b.value - a.value);
    return d3.partition().size([2 * Math.PI, root.height + 1])(root);
  };

  const appendTooSmall = (small, d) => {
    const parent = d.parent.data.name;
    if (small.hasOwnProperty(parent)) {
      if (small[parent]["start"] < arc.startAngle()(d)) {
        small[parent]["start"] = arc.startAngle()(d);
        small[parent]["startNode"] = d;
        small[parent]["cellCount"] = small[parent]["cellCount"] + 1;
      }
      if (small[parent]["end"] > arc.startAngle()(d)) {
        small[parent]["end"] = arc.endAngle()(d);
        small[parent]["endNode"] = d;
        small[parent]["cellCount"] = small[parent]["cellCount"] + 1;
      }
    } else {
      small[parent] = {
        cellCount: 1,
        start: arc.startAngle()(d),
        end: arc.endAngle()(d),
        startNode: d,
        endNode: d,
      };
    }
    return small;
  };

  const drawBreadCrumbs = (svg, height) => {
    var trail = svg
      .append("svg")
      .attr("width", 200)
      .attr("height", 300)
      .attr("y", 10)
      .attr("id", "trail");
    trail.append("text").attr("id", "endlabel").style("fill", "#000");
  };

  const drawArea = (svg, setIsTooltipOpen, setHoveredItem, colors) => {
    const root = partition(hierarchy);
    root.each((d) => (d.current = d));

    const g = svg
      .append("g")
      .attr("transform", `translate(${width / 2},${height / 2})`);

    var small = {};

    d3.select(tooltipRef.current).style("pointer-events", "none");

    const path = g
      .append("g")
      .selectAll("path")
      .data(
        root.descendants().filter(function (d) {
          const isTooSmall =
            Math.abs(arc.startAngle()(d) - arc.endAngle()(d)) > ANGLE;
          if (!isTooSmall && d.parent && d.height === 1) {
            small = appendTooSmall(small, d);
          }
          return isTooSmall;
        })
      )
      .join("path")
      .attr("fill", (d) => {
        while (d.depth > 1) d = d.parent;
        return selectedCloneColor
          ? selectedCloneColor
          : colorScale(d.data.name);
      })
      .attr("fill-opacity", (d) => {
        return arcVisible(d.current) ? (d.children ? 0.6 : 0.4) : 0;
      })
      .attr("d", (d) => arc(d.current))
      .on("mouseover", mouseover)
      .on("mouseout", mouseout);

    const smallList = Object.keys(small)
      .filter((d) => d !== "flair")
      .map((d) => {
        var curr = small[d]["startNode"].current;
        curr["cellCount"] = small[d]["cellCount"];
        curr["x0"] = small[d]["endNode"].current["x0"];
        curr["y1"] = small[d]["endNode"].current["y1"];
        return curr;
      })
      .map((d) => (d.current = d));

    const smallListTitleOmit = smallList.map((d) => d.data.name);

    const joinedPath = g
      .append("g")
      .selectAll("path")
      .data(smallList)
      .join("path")
      .attr("fill", (d) => {
        while (d.depth > 1) d = d.parent;
        return colorScale(d.data.name);
      })
      .attr("fill-opacity", (d) => {
        return arcVisible(d.current) ? (d.children ? 0.6 : 0.4) : 0;
      })
      .attr("d", (d) => arc(d))
      .on("mouseover", mouseover)
      .on("mouseout", mouseout);

    path
      .filter((d) => d.children)
      .style("cursor", "pointer")
      .on("click", clicked);

    path.append("title").text(
      (d) =>
        `${d
          .ancestors()
          .map((d) => d.data.name)
          .reverse()
          .join("/")}\n${format(d.value)}`
    );

    const label = g
      .append("g")
      .attr("pointer-events", "none")
      .attr("text-anchor", "middle")
      .style("user-select", "none")
      .selectAll("text")
      .data(
        root.descendants().filter(function (d) {
          return Math.abs(arc.startAngle()(d) - arc.endAngle()(d)) > ANGLE;
        })
      )
      .join("text")
      .attr("font-size", 10)
      .attr("dy", "0.35em")
      .attr("fill-opacity", (d) => +labelVisible(d.current))
      .attr("transform", (d) => labelTransform(d.current))
      .text((d) => {
        if (smallListTitleOmit.indexOf(d.data.name) !== -1) {
          return "";
        } else {
          return d.data.name;
        }
      });

    const parent = g
      .append("circle")
      .datum(root)
      .attr("r", radius)
      .attr("fill", "none")
      .attr("pointer-events", "all")
      .on("click", clicked);

    g.append("text")
      .attr("dx", -35)
      .attr("dy", 5)
      .attr("pointer-events", "none")
      .text("zoom out")
      .attr("display", "none")
      .attr("id", "zoomText");

    function breadcrumbPoints(d, i) {
      var points = [];
      points.push("0,0");
      points.push(b.w + ",0");
      points.push(b.w + b.t + "," + b.h / 2);
      points.push(b.w + "," + b.h);
      points.push("0," + b.h);
      if (i > 0) {
        points.push(b.t + "," + b.h / 2);
      }
      return points.join(" ");
    }
    function updateBreadcrumbs(nodeArray, percentageString, data, isClicked) {
      var g = d3
        .select("#trail")
        .selectAll("g")
        .data(nodeArray, function (d) {
          return d.data.name + d.depth;
        })
        .attr("id", isClicked ? "clicked" : "");

      // Add breadcrumb and label for entering nodes.
      var entering = g
        .enter()
        .append("g")
        .attr("transform", (d, i) => "translate(10, " + i * (b.h + b.s) + ")");

      entering
        .append("polygon")
        .attr("points", breadcrumbPoints)
        .style("fill", function (d) {
          if (d.depth === 1) {
            return colorScale(d.data.name);
          } else {
            var color = rgb(colorScale(d.parent.data.name));
            color.opacity = 0.6;
            return color;
          }
        });

      entering
        .append("text")
        .attr("x", (b.w + b.t) / 2)
        .attr("y", (d) => {
          if (d.depth === 2 && !parentIsTopTen(d)) {
            return -b.h;
          } else {
            return b.h / 2;
          }
        })
        .attr("dy", "0.35em")
        .attr("text-anchor", "middle")
        .attr("font-size", 12)
        .text(function (d) {
          if (parentIsTopTen(d)) {
            return "Clone " + d.data.name;
          } else if (
            d.parent &&
            (d.parent.data.name === "Other" ||
              d.parent.data.name === "Singleton")
          ) {
            return "";
          } else {
            return d.data.name + " Clones";
          }
        });

      g.exit().remove();
      const yPos =
        nodeArray.length !== 1 && !parentIsTopTen(nodeArray[1])
          ? -(b.h + 10) + (nodeArray.length + 0.5) * (b.h + b.s)
          : (nodeArray.length + 0.5) * (b.h + b.s);

      //  if (typeof percentageString !== "function") {
      d3.select("#trail")
        .select("#endlabel")
        .attr("x", b.w / 2 + 10)
        .attr("y", yPos)
        .attr("dy", "0.35em")
        .attr("text-anchor", "middle")
        .text(percentageString + " cells");
      //  }
      // Make the breadcrumb trail visible, if it's hidden.
      d3.select("#trail").attr("display", "all");
    }

    function mouseout(d, i) {
      if (!d3.select("#trail #clicked")) {
        d3.select("#trail").attr("display", "none");
      }
      setHoverItem(null);
      setIsTooltipOpen(false);
      setRootParent(null);
    }
    function getAncestors(node) {
      var path = [];
      var current = node;
      while (current.parent) {
        path.unshift(current);
        current = current.parent;
      }
      return path;
    }
    function mouseover(d, i) {
      const arcLocation = arc.centroid(d.current, i);
      setIsTooltipOpen(true);
      setHoveredItem(d.data.name);
      setRootParent(d.parent.data.name);

      d3.select(tooltipRef.current)
        .transition()
        .style("left", arcLocation[0] + width / 2 + "px")
        .style("top", arcLocation[1] + height / 2 + "px");

      var sizeCount = d.data.cellCount
        ? d.data.cellCount
        : d.cellCount
        ? d.cellCount
        : d.data.children.length;

      var sequenceArray = getAncestors(d);
      updateBreadcrumbs(sequenceArray, sizeCount, d);
    }

    function clicked(p) {
      parent.datum(p.parent || root);

      root.each(
        (d) =>
          (d.target = {
            x0:
              Math.max(0, Math.min(1, (d.x0 - p.x0) / (p.x1 - p.x0))) *
              2 *
              Math.PI,
            x1:
              Math.max(0, Math.min(1, (d.x1 - p.x0) / (p.x1 - p.x0))) *
              2 *
              Math.PI,
            y0: Math.max(0, d.y0 - p.depth),
            y1: Math.max(0, d.y1 - p.depth),
          })
      );

      const t = g.transition().duration(750);

      joinedPath
        .transition(t)
        .tween("data", (d) => {
          const i = d3.interpolate(d.current, d.target);
          return (t) => (d.current = i(t));
        })
        .filter(function (d) {
          return +this.getAttribute("fill-opacity") || arcVisible(d.target);
        })
        .attr("fill-opacity", (d) =>
          arcVisible(d.target) ? (d.children ? 0.6 : 0.4) : 0
        )
        .attrTween("d", (d) => () => arc(d.current));

      path
        .transition(t)
        .tween("data", (d) => {
          const i = d3.interpolate(d.current, d.target);
          return (t) => (d.current = i(t));
        })
        .filter(function (d) {
          return +this.getAttribute("fill-opacity") || arcVisible(d.target);
        })
        .attr("fill-opacity", (d) =>
          arcVisible(d.target) ? (d.children ? 0.6 : 0.4) : 0
        )
        .attrTween("d", (d) => () => arc(d.current));

      label
        .filter(function (d) {
          return +this.getAttribute("fill-opacity") || labelVisible(d.target);
        })
        .transition(t)
        .attr("fill-opacity", (d) => +labelVisible(d.target))
        .attrTween("transform", (d) => () => labelTransform(d.current));

      //d3.select("#zoomText").
      var sequenceArray = getAncestors(p);
      var sizeCount = p.data.cellCount
        ? p.data.cellCount
        : p.cellCount
        ? p.cellCount
        : p.data.children.length;

      updateBreadcrumbs(sequenceArray, sizeCount, p, true);
    }

    function arcVisible(d) {
      return d.y1 <= 3 && d.y0 >= 1 && d.x1 > d.x0;
    }

    function labelVisible(d) {
      return d.y1 <= 3 && d.y0 >= 1 && (d.y1 - d.y0) * (d.x1 - d.x0) > 0.09;
    }

    function labelTransform(d) {
      const x = (((d.x0 + d.x1) / 2) * 180) / Math.PI;
      const y = ((d.y0 + d.y1) / 2) * radius;
      return `rotate(${x - 90}) translate(${y},0) rotate(${x < 180 ? 0 : 180})`;
    }
  };

  const ref = useD3(
    (svg) => {
      drawArea(svg, setIsTooltipOpen, setHoverItem, colors);

      drawBreadCrumbs(svg, height);
    },
    width,
    height,
    [data]
  );

  return (
    <div style={{ position: "relative", marginTop: -145 }}>
      <Title title="Clone Distribution" />
      <div id="tooltipDiv" ref={tooltipRef} className={classes.tooltipWrapper}>
        <Tooltip
          PopperProps={{
            disablePortal: true,
            style: { pointerEvents: "none" },
          }}
          classes={{ tooltip: classes.tooltip }}
          arrow
          title={
            <TooltipText
              hoverParent={hoverParent}
              allSubsets={allSubsets}
              hoverItem={hoverItem}
              otherSubsetParam={otherSubsetParam}
              hierarchy={hierarchy}
              totalCount={totalCount}
            />
          }
          open={isTooltipOpen}
          disableFocusListener
          disableHoverListener
          disableTouchListener
        >
          <div
            style={{
              position: "absolute",
              x: 100,
              y: 100,
              width: 5,
              height: 5,
              fill: "blue",
            }}
          />
        </Tooltip>
      </div>
      <svg ref={ref} />
    </div>
  );
};
const TooltipText = ({
  allSubsets,
  hoverItem,
  otherSubsetParam,
  hoverParent,
  totalCount,
  hierarchy,
}) => {
  return (
    <Fragment>
      {hoverItem && (
        <span>
          <Typography color="inherit" style={{ fontSize: 15 }}>
            {hoverParent &&
            hoverParent !== "flair" &&
            hoverParent !== "Singleton"
              ? "Clone " + hoverItem
              : hoverItem + " Clones"}
          </Typography>
          {hoverParent && hoverParent !== "flair" && (
            <span style={{ fontSize: 13 }}>
              <span>
                <b>
                  {Object.keys(allSubsets[hoverItem]).reduce(
                    (final, curr) => final + allSubsets[hoverItem][curr].length,
                    0
                  )}
                </b>
                <em>{" cells"}</em>
              </span>
              <div>
                <b>
                  {Math.round(
                    (allSubsets[hoverItem].length / totalCount) * 100000
                  ) / 1000}
                </b>
                <em>{"%"}</em>
              </div>
            </span>
          )}
          {hoverParent === "flair" && (
            <span>
              <b>
                {hierarchy["children"]
                  .filter((d) => d["name"] === hoverItem)[0]
                  ["children"].reduce(
                    (final, d) => final + d["children"].length,
                    0
                  )}
              </b>
              <em>{" / " + totalCount + " cells"}</em>
            </span>
          )}
        </span>
      )}
    </Fragment>
  );
};
export default Sunburst;
