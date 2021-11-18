import React, { useState, useRef, Fragment } from "react";

import * as d3 from "d3";
import { useD3, DownloadIcon } from "@shahlab/planetarium";
import _ from "lodash";

import makeStyles from "@mui/styles/makeStyles";

import { Layout } from "@shahlab/planetarium";
import Tooltip from "@mui/material/Tooltip";
import Typography from "@mui/material/Typography";

import { INFO } from "../config";

import { jsPDF } from "jspdf";
import canvg from "canvg";

const b = {
  w: 75,
  h: 30,
  s: 3,
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
const Sunburst = ({
  data,
  colors,
  width,
  height,
  subsetParam,
  type,
  otherSubsetParam,
  cloneColorScale,
  selectedCloneColor = null,
}) => {
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
  const resetSubsetcounts = subsetCounts.slice(10, subsetCounts.length);

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
  const topTenObj = { name: "Top Ten", children: [...topTen] };
  const singleObj = { name: "Singleton", children: [...singletonCounts] };
  const otherObj = { name: "Other", children: [...restCounts] };
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
      }
      if (small[parent]["end"] > arc.startAngle()(d)) {
        small[parent]["end"] = arc.endAngle()(d);
        small[parent]["endNode"] = d;
      }
    } else {
      small[parent] = {
        start: arc.startAngle()(d),
        end: arc.endAngle()(d),
        startNode: d,
        endNode: d,
      };
    }
    return small;
  };

  const drawBreadCrumbs = (svg) => {
    var trail = svg
      .append("svg")
      .attr("width", 100)
      .attr("height", 300)
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

    function breadcrumbPoints(d, i) {
      var points = [];
      points.push("0,0");
      points.push(b.w + ",0");
      points.push(b.w + b.t + "," + b.h / 2);
      points.push(b.w + "," + b.h);
      points.push("0," + b.h);
      if (i > 0) {
        // Leftmost breadcrumb; don't include 6th vertex.
        points.push(b.t + "," + b.h / 2);
      }
      // debugger;
      return points.join(" ");
    }
    function updateBreadcrumbs(nodeArray, percentageString) {
      // Data join; key function combines name and depth (= position in sequence).
      var g = d3
        .select("#trail")
        .selectAll("g")
        .data(nodeArray, function (d) {
          return d.data.name + d.depth;
        });

      // Add breadcrumb and label for entering nodes.
      var entering = g
        .enter()
        .append("g")
        .attr("transform", function (d, i) {
          return "translate(0, " + i * (b.h + b.s) + ")";
        });

      entering
        .append("polygon")
        .attr("points", breadcrumbPoints)
        .style("fill", function (d) {
          console.log(d);
          return colorScale(d.data.name);
        });

      entering
        .append("text")
        .attr("x", (b.w + b.t) / 2)
        .attr("y", b.h / 2)
        .attr("dy", "0.35em")
        .attr("text-anchor", "middle")
        .text(function (d) {
          return d.data.name;
        });

      // Set position for entering and updating nodes.

      // Remove exiting nodes.
      g.exit().remove();

      // Now move and update the percentage at the end.
      d3.select("#trail")
        .select("#endlabel")
        .attr("x", b.w / 2)
        .attr("y", (nodeArray.length + 0.5) * (b.h + b.s))
        .attr("dy", "0.35em")
        .attr("text-anchor", "middle")
        .text("Count:" + percentageString);

      // Make the breadcrumb trail visible, if it's hidden.
      d3.select("#trail").style("visibility", "");
    }

    function mouseout(d, i) {
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
      console.log(totalCount);
      var sizeCount = d.data.children.length;

      var sequenceArray = getAncestors(d);

      updateBreadcrumbs(sequenceArray, sizeCount);
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

      drawBreadCrumbs(svg);
    },
    width,
    height,
    [data]
  );
  return (
    <Layout
      addIcon={
        <DownloadIcon download={async () => download(ref, width, height)} />
      }
      title={INFO[type]["title"]}
      infoText={INFO[type]["text"]}
    >
      <div style={{ position: "relative" }}>
        <div
          id="tooltipDiv"
          ref={tooltipRef}
          className={classes.tooltipWrapper}
        >
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
    </Layout>
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
            {hoverItem}
          </Typography>
          {hoverParent && hoverParent !== "flair" && (
            <span style={{ fontSize: 13 }}>
              <span>
                <b>{allSubsets[hoverItem].length}</b>
                <em>{" - data points"}</em>
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
              <em>{" - data points"}</em>
            </span>
          )}
        </span>
      )}
    </Fragment>
  );
};
export default Sunburst;
