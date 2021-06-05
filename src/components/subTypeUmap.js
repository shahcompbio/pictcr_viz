import React from "react";
import * as d3 from "d3";
import { quantileSorted } from "d3-array";
import _ from "lodash";
import { INFO } from "../config.js";

import {
  Layout,
  VerticalLegend,
  useCanvas,
  isValueHighlighted as isHighlighted,
} from "@shahlab/planetarium";

import Grid from "@material-ui/core/Grid";

import { CONSTANTS } from "../config";

const PADDING = 10;

const LEGEND_WIDTH = 180;
const AXIS_SPACE = 20;
const AXIS_FONT = "normal 10px Helvetica";
const AXIS_COLOR = "#000000";

const LABEL_FONT = "500 12px Helvetica";

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

const NULL_POINT_COLOR = "#e8e8e8";
const POINT_RADIUS = 2;
const PERCENTILE_RANGE = [0.25, 0.75];

const DataWrapper = ({
  chartName,
  data,
  chartDim,
  selectedSubtype,
  hoveredSubtype,
  setSelectedSubtype,
}) => {
  const { xParam, yParam, subtypeParam } = CONSTANTS;

  const setHighlighted = (event, value) => {
    if (event === "mouseenter") {
      setSelectedSubtype({ hover: value });
    } else if (event === "mousedown") {
      setSelectedSubtype({
        hover: null,
        selected: value,
      });
    } else if (event === "mouseout") {
      setSelectedSubtype({
        hover: null,
      });
    }
  };
  return (
    <Layout
      title={INFO["SUBTYPEUMAP"]["title"]}
      infoText={INFO["SUBTYPEUMAP"]["text"]}
    >
      <UMAP
        width={chartDim["width"]}
        height={chartDim["height"]}
        chartName={chartName}
        data={data}
        highlighted={hoveredSubtype || selectedSubtype}
        xParam={xParam}
        yParam={yParam}
        subsetParam={subtypeParam}
        setHighlighted={setHighlighted}
      />
    </Layout>
  );
};

const UMAP = ({
  width,
  height,
  data,
  highlighted,
  xParam,
  yParam,
  subsetParam,
  setHighlighted,
}) => {
  const canvasWidth = width - LEGEND_WIDTH;
  const canvasHeight = height;

  const chartWidth = canvasWidth - AXIS_SPACE;
  const chartHeight = canvasHeight - AXIS_SPACE - PADDING - PADDING;

  const yData = data.map((d) => parseFloat(d[yParam]));
  const xData = data.map((d) => parseFloat(d[xParam]));

  const yMin = Math.min(...yData);
  const yMax = Math.max(...yData);
  const xMin = Math.min(...xData);
  const xMax = Math.max(...xData);

  const xScale = d3
    .scaleLinear()
    .domain([xMin, xMax])
    .range([PADDING, PADDING + chartWidth]);

  const yScale = d3
    .scaleLinear()
    .domain([yMax, yMin])
    .range([PADDING, PADDING + chartHeight]);

  const subsetGroups = _.groupBy(data, subsetParam);
  const subsetValues = Object.keys(subsetGroups).sort();

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
      drawUMAPAxis(context, canvasHeight - AXIS_SPACE, xParam, yParam);
      drawPoints(
        context,
        data,
        xScale,
        yScale,
        xParam,
        yParam,
        subsetParam,
        highlighted,
        subsetColors
      );

      drawSubsetLabels(
        context,
        subsetGroups,
        xScale,
        yScale,
        xParam,
        yParam,
        highlighted
      );
    },
    canvasWidth,
    canvasHeight,
    [highlighted]
  );

  return (
    <Grid container direction="row" style={{ padding: 0 }}>
      <Grid item>
        <canvas ref={canvasRef} />
      </Grid>
      <Grid item>
        <VerticalLegend
          width={LEGEND_WIDTH}
          height={chartHeight}
          labels={subsetLabels}
          setHighlighted={setHighlighted}
        />
      </Grid>
    </Grid>
  );
};

const drawPoints = (
  context,
  data,
  xScale,
  yScale,
  xParam,
  yParam,
  subsetParam,
  highlighted,
  colorScale
) => {
  context.beginPath();
  context.lineWidth = 1;
  context.globalAlpha = 1;

  data.forEach((point) => {
    context.fillStyle = isHighlighted(point[subsetParam], highlighted)
      ? colorScale(point[subsetParam])
      : NULL_POINT_COLOR;

    context.beginPath();
    context.arc(
      xScale(point[xParam]),
      yScale(point[yParam]),
      POINT_RADIUS,
      0,
      Math.PI * 2,
      true
    );
    context.fill();
  });
  if (highlighted) {
    data
      .filter((row) => row[subsetParam] === highlighted)
      .forEach((point) => {
        context.fillStyle = colorScale(point[subsetParam]);
        context.beginPath();
        context.arc(
          xScale(point[xParam]),
          yScale(point[yParam]),
          POINT_RADIUS,
          0,
          Math.PI * 2,
          true
        );
        context.fill();
      });
  }
};

const drawUMAPAxis = (context, chartHeight, xParam, yParam) => {
  context.beginPath();
  context.font = AXIS_FONT;
  context.globalAlpha = 1;

  const START_X = AXIS_SPACE / 2;
  const START_Y = chartHeight + AXIS_SPACE / 2;

  context.fillStyle = AXIS_COLOR;
  context.strokeStyle = AXIS_COLOR;
  context.moveTo(START_X, START_Y);
  context.lineTo(START_X, START_Y - 50);
  context.stroke();

  context.beginPath();
  context.moveTo(START_X, START_Y);
  context.lineTo(START_X + 50, START_Y);
  context.stroke();

  context.textAlign = "left";
  context.textBaseline = "middle";
  context.fillText(xParam, START_X + 52, START_Y);
  context.save();
  context.rotate((270 * Math.PI) / 180);
  context.fillText(yParam, -(START_Y - 52), START_X);
  context.restore();
};

const drawSubsetLabels = (
  context,
  subsetGroups,
  xScale,
  yScale,
  xParam,
  yParam,
  highlighted
) => {
  const subsetValues = Object.keys(subsetGroups);

  subsetValues.forEach((subset) => {
    const subsetData = subsetGroups[subset];

    const filteredData = filterOutliers(subsetData, xParam, yParam);
    const { xMin, xMax, yMin, yMax } = getBoxBounds(
      filteredData,
      xParam,
      yParam
    );

    const [x1, x2, y1, y2] = [
      xScale(xMin),
      xScale(xMax),
      yScale(yMin),
      yScale(yMax),
    ];
    const width = Math.abs(x2 - x1);
    const height = Math.abs(y2 - y1);

    context.font = LABEL_FONT;
    const textWidth = context.measureText(subset).width;
    context.globalAlpha = isHighlighted(subset, highlighted) ? 0.8 : 0.2;
    context.fillStyle = "white";
    context.fillRect(
      x1 + (width - textWidth) / 2 - 1,
      y1 - height / 2 - 7,
      textWidth + 2,
      14
    );

    context.globalAlpha = isHighlighted(subset, highlighted) ? 1 : 0.2;
    context.textAlign = "center";
    context.textBaseline = "middle";
    context.fillStyle = "black";
    context.fillText(subset, x1 + width / 2, y1 - height / 2);
  });
};

const filterOutliers = (data, xParam, yParam) => {
  const xValues = data
    .map((datum) => parseFloat(datum[xParam]))
    .sort((a, b) => a - b);
  const yValues = data
    .map((datum) => parseFloat(datum[yParam]))
    .sort((a, b) => a - b);

  const [xMin, xMax] = PERCENTILE_RANGE.map((range) =>
    quantileSorted(xValues, range)
  );
  const [yMin, yMax] = PERCENTILE_RANGE.map((range) =>
    quantileSorted(yValues, range)
  );

  const xiqr = Math.abs(xMax - xMin) * 1.5;
  const yiqr = Math.abs(yMax - yMin) * 1.5;

  return data.filter((datum) => {
    const x = datum[xParam];
    const y = datum[yParam];

    return (
      xMin - xiqr <= x &&
      x <= xMax + xiqr &&
      yMin - yiqr <= y &&
      y <= yMax + yiqr
    );
  });
};

const getBoxBounds = (data, xParam, yParam) => {
  const xValues = data.map((datum) => datum[xParam]);
  const yValues = data.map((datum) => datum[yParam]);

  const xMin = Math.min(...xValues);
  const xMax = Math.max(...xValues);
  const yMin = Math.min(...yValues);
  const yMax = Math.max(...yValues);

  return { xMin, xMax, yMin, yMax };
};

export default DataWrapper;
