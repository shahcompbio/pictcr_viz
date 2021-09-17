import React, { useState, useRef, Fragment } from "react";

import * as d3 from "d3";
import { useD3, DownloadIcon } from "@shahlab/planetarium";
import _ from "lodash";

import { makeStyles } from "@material-ui/core/styles";

import { Layout } from "@shahlab/planetarium";
import Tooltip from "@material-ui/core/Tooltip";
import Typography from "@material-ui/core/Typography";
import Table from "@material-ui/core/Table";
import TableBody from "@material-ui/core/TableBody";
import TableCell from "@material-ui/core/TableCell";
import TableRow from "@material-ui/core/TableRow";

import { INFO } from "../config";

import { jsPDF } from "jspdf";
import canvg from "canvg";

const TOP_NUM = 3;

const useStyles = makeStyles((theme) => ({
  tooltip: {
    minWidth: 220,
  },
}));
const download = async (ref, width, height) => {
  var doc = new jsPDF("L", "px", [width, height]);

  const newCanvas = document.createElement("canvas");
  const context = newCanvas.getContext("2d");

  //const currCanvas = ref.current;
  //const currContext = canvas.getContext("2d");

  let scale = window.devicePixelRatio;
  newCanvas.style.width = width + "px";
  newCanvas.style.height = height + "px";
  newCanvas.width = width * scale;
  newCanvas.height = height * scale;

  context.scale(scale, scale);

  const plot = d3.select(ref.current).node();
  //const profile = canvg.fromString(context, plot.outerHTML, { useCORS: true });

  // Render only first frame, ignoring animations.
  //await profile.render();

  //const png = newCanvas.toBuffer();

  //  const imgProfileBackground = context.drawImage(img, 0, 0, width, height);
  //  context.stroke(path);
  //const currContext = canvas.getContext("2d");

  const plotOuter = plot.outerHTML;
  const profile = canvg.fromString(context, plotOuter);
  profile.start();

  const imgProfileBackground = newCanvas.toDataURL("image/png", 1.0);
  //const docWidth = Math.round((width * 25.4) / 96);
  //const docHeight = Math.round((height * 25.4) / 96);
  doc.addImage(
    imgProfileBackground,
    "PNG",
    0,
    0,
    width / scale,
    height / scale
  );

  doc.save("test.pdf");
};
const Doughnut = ({
  data,
  colors,
  width,
  height,
  subsetParam,
  type,
  otherSubsetParam,
  Select,
}) => {
  const [isTooltipOpen, setIsTooltipOpen] = useState(false);
  const tooltipRef = useRef(null);
  const [hoverItem, setHoverItem] = useState(null);

  const classes = useStyles();

  const totalCount = data.length;

  const allSubsets = _.groupBy(data, (datum) => datum[subsetParam]);

  const subsetCounts = Object.keys(allSubsets)
    .sort((a, b) => allSubsets[b].length - allSubsets[a].length)
    .slice(0, 10);

  const colorScale = d3
    .scaleOrdinal()
    .domain(subsetCounts)
    .range(colors.slice(0, Math.min(subsetCounts.length, colors.length)));

  const subsets = subsetCounts.map((d) => ({
    key: d,
    value: allSubsets[d],
  }));

  const radius = Math.min(width, height) / 2;
  const arc = d3
    .arc()
    .innerRadius(radius * 0.3)
    .outerRadius(radius * 0.65);

  const arcLabel = d3
    .arc()
    .innerRadius(radius * 0.9)
    .outerRadius(radius * 0.9);

  const drawArea = (svg, setIsTooltipOpen, setHoveredItem) => {
    const pie = d3
      .pie()
      .value((d) => d["value"].length)
      .sort(null);
    const arcs = pie(subsets);

    const path = svg
      .append("g")
      .attr("transform", "translate(" + width / 2 + "," + height / 2 + ")")
      .attr("stroke", "white")
      .selectAll("path")
      .data(arcs, (d) => d.data.key);

    /*  path
      .join("path")
      .attr("fill", (d) => colorScale(d["data"].key))
      .attr("d", arc)
      .on("mouseover", mouseover)
      .on("mouseout", mouseout);*/

    path
      .join("path", (update) =>
        update.on("mouseover", mouseover).on("mouseout", mouseout)
      )
      .attr("fill", (d) => colorScale(d["data"].key))
      .attr("d", arc);

    svg
      .append("g")
      .attr("transform", "translate(" + width / 2 + "," + height / 2 + ")")
      .attr("font-family", "MyFontLight")
      .attr("font-size", 12)
      .attr("text-anchor", "middle")
      .selectAll("text")
      .data(arcs)
      .join("text")
      .attr("id", (d) => "label-" + d.data.key)
      .attr("transform", (d, i) => {
        return `translate(${arcLabel.centroid(d, i)})`;
      })
      .call((text) =>
        text
          .append("tspan")
          .attr("y", "-0.4em")
          .attr("font-weight", "bold")
          .text((d) => d.data.name)
      )
      .call((text) =>
        text
          .filter((d) => d.endAngle - d.startAngle > 0.25)
          .append("tspan")

          .attr("fill-opacity", 0.7)
          .text(
            (d) =>
              d.data.key.toLocaleString() +
              " (" +
              d3.format(".0%")(d.data.value.length / totalCount) +
              ")"
          )
      );
    function mouseover(element, index, elements) {
      d3.select(elements[index])
        .transition()
        .duration(500)
        .attr("transform", function (d, i) {
          var x;
          var y;
          if (element.data._translate) {
            x = element.data._translate.x;
            y = element.data._translate.y;
          } else if (!element.data._expanded) {
            element.data._expanded = true;
            var a =
              element.startAngle +
              (element.endAngle - element.startAngle) / 2 -
              Math.PI / 2;
            x = Math.cos(a) * 10;
            y = Math.sin(a) * 10;
            element.data._translate = { x: x, y: y };
          }
          const arcLocation = arcLabel.centroid(element, i);
          setIsTooltipOpen(true);
          setHoveredItem(element.data.key);

          d3.select(tooltipRef.current)
            .transition()
            .style("left", arcLocation[0] + width / 2 + "px")
            .style("top", arcLocation[1] + height / 2 + "px");

          return "translate(" + x + "," + y + ")";
        });
      d3.select("#label-" + element.data.key).style("font-weight", "bold");
    }
    function mouseout(d, index, elements) {
      d3.select(elements[index])
        .transition()
        .delay(200)
        .duration(500)
        .attr("transform", function (e) {
          setIsTooltipOpen(false);
          d.data._expanded = false;
        });
      d3.select("#label-" + d.data.key).style("font-weight", "normal");
    }
  };

  const ref = useD3(
    (svg) => {
      drawArea(svg, setIsTooltipOpen, setHoverItem);
    },
    width,
    height,
    [data, subsetParam]
  );
  return (
    <Layout
      addIcon={
        <DownloadIcon download={async () => download(ref, width, height)} />
      }
      title={INFO[type]["title"]}
      infoText={INFO[type]["text"]}
      addIcon={Select}
    >
      <div style={{ position: "relative" }}>
        <div
          id="tooltipDiv"
          ref={tooltipRef}
          style={{ position: "relative", width: "10px", height: "10px" }}
        >
          <Tooltip
            PopperProps={{
              disablePortal: true,
            }}
            classes={{ tooltip: classes.tooltip }}
            arrow
            title={
              <TooltipText
                allSubsets={allSubsets}
                hoverItem={hoverItem}
                otherSubsetParam={otherSubsetParam}
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
    </Layout>
  );
};
const TooltipText = ({ allSubsets, hoverItem, otherSubsetParam }) => (
  <Fragment>
    {hoverItem && (
      <span>
        <Typography color="inherit">{hoverItem}</Typography>
        <b>{allSubsets[hoverItem].length}</b>
        <em>{" - data points"}</em>
        <div>
          <Typography style={{ fontSize: 15 }} gutterBottom>
            Top {TOP_NUM}{" "}
            {otherSubsetParam === "subtype" ? "Subtypes" : "Clone IDs"}
          </Typography>
          <Table size="small">
            <TableBody>
              {_.chain(allSubsets[hoverItem])
                .groupBy(otherSubsetParam)
                .orderBy(
                  function (o) {
                    return o.length;
                  },
                  ["desc", "desc"]
                )
                .filter(function (o, i) {
                  return i < 3;
                })
                .value()
                .map(function (d) {
                  return (
                    <TableRow
                      key={d[0][otherSubsetParam]}
                      style={{ paddingBottom: "5px" }}
                    >
                      <TableCell
                        style={{
                          borderBottom: "none",
                          paddingRight: 0,
                          color: "#eeeeeede",
                        }}
                      >
                        <b>Clone {d[0][otherSubsetParam]}</b>
                      </TableCell>
                      <TableCell
                        style={{ borderBottom: "none", color: "white" }}
                      >
                        {d.length} data points -{" "}
                        {d3.format(".001%")(
                          d.length / allSubsets[hoverItem].length
                        )}
                      </TableCell>
                    </TableRow>
                  );
                })}
            </TableBody>
          </Table>
        </div>
      </span>
    )}
  </Fragment>
);

export default Doughnut;
