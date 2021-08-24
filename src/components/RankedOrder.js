import React from "react";
import * as d3 from "d3";
import _ from "lodash";

import {
  useCanvas,
  VerticalLegend,
  drawCanvasAxis,
} from "@shahlab/planetarium";

import { Grid } from "@material-ui/core";

const COLOR_ARRAY = [
  "#5E4FA2",
  "#3288BD",
  "#66C2A5",
  "#FEE08B",
  "#FDAE61",
  "#F46D43",
  "#D53E4F",
  "#c9cc76",
  "#9E0142",
  "#C6AEFF",
  "#BDD8FF",
  "#BDFFB2",
  "#FFC8AE",
  "#FF9FBB",
  "#b2dbd6",
  "#ffd470",
];
const PADDING = 50;
const format = d3.format(".3f");

const RankedOrder = ({ width, height, data }) => {
  // need to calculate rank
  const subsets = _.groupBy(data, "subtype");
  const subsetValues = Object.keys(subsets).sort();
  // group by clonotype, get counts, create records, then sort
  const records = subsetValues.map((subsetName) => {
    const subsetData = subsets[subsetName];

    const cloneGroups = _.groupBy(subsetData, "clone_id");
    const cloneValues = Object.keys(cloneGroups).sort();

    return cloneValues
      .map((cloneName) => ({
        clone_id: cloneName,
        subtype: subsetName,
        frequency: cloneGroups[cloneName].length / subsetData.length,
      }))
      .sort((a, b) => b.frequency - a.frequency);
  });

  const maxRank = records
    .map((grouping) => grouping.length)
    .reduce((rsf, record) => Math.max(rsf, record), 0);
  const maxFreq = records
    .map((grouping) => grouping[0].frequency)
    .reduce((rsf, record) => Math.max(rsf, record), 0);
  const minFreq = records
    .map((grouping) => grouping[grouping.length - 1].frequency)
    .reduce((rsf, record) => Math.min(rsf, record), 1);
  const rankScale = d3
    .scaleLog()
    .domain([1, maxRank])
    .range([PADDING, width - PADDING]);
  const freqScale = d3
    .scaleLog()
    .domain([minFreq, maxFreq])
    .nice()
    .range([height - PADDING, PADDING]);
  const subsetColors = d3
    .scaleOrdinal()
    .domain(subsetValues)
    .range(
      COLOR_ARRAY.slice(0, Math.min(subsetValues.length, COLOR_ARRAY.length))
    );

  const subsetLabels = subsetValues.map((value) => ({
    value,
    label: value,
    color: subsetColors(value),
  }));

  const canvasRef = useCanvas(
    (canvas) => {
      const context = canvas.getContext("2d");

      // draw axis
      drawCanvasAxis({
        context,
        xScale: rankScale,
        yScale: freqScale,
        ticks: 3,
        label: "Clone Frequency",
      });

      drawCanvasAxis({
        context,
        xScale: rankScale,
        yScale: freqScale,
        ticks: 3,
        label: "Rank",
        orientation: "horizontal",
      });

      // draw plot
      context.beginPath();
      context.lineWidth = 2;
      context.globalAlpha = 1;

      subsetValues.forEach((subsetName, index) => {
        const subsetData = records[index];
        context.beginPath();

        context.fillStyle = subsetColors(subsetName);
        context.strokeStyle = subsetColors(subsetName);
        context.moveTo(rankScale(1), freqScale(subsetData[0]["frequency"]));
        subsetData.forEach((record, index) => {
          // context.arc(
          //   rankScale(index + 1),
          //   freqScale(record["frequency"] / total),
          //   POINT_RADIUS,
          //   0,
          //   Math.PI * 2,
          //   true
          // );
          // context.fill();
          context.lineTo(rankScale(index + 1), freqScale(record["frequency"]));
          context.stroke();
        });
      });
    },
    width,
    height,
    []
  );

  return (
    <Grid container direction="row" style={{ padding: 0 }}>
      <Grid item>
        <canvas ref={canvasRef} />
      </Grid>
      <Grid item>
        <VerticalLegend
          width={200}
          height={height}
          labels={subsetLabels}
          setHighlighted={() => {}}
        />
      </Grid>
    </Grid>
  );
};

export default RankedOrder;
