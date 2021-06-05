import React, { useState } from "react";
import * as d3 from "d3";
import * as d3Hexbin from "d3-hexbin";
import { INFO } from "../config.js";

import {
  Layout,
  VerticalLegend,
  useCanvas,
  isValueHighlighted as isHighlighted,
} from "@shahlab/planetarium";

import Grid from "@material-ui/core/Grid";
import Slider from "@material-ui/core/Slider";

import { CONSTANTS } from "../config";

import _ from "lodash";

const PADDING = 10;
const AXIS_SPACE = 20;
const LEGEND_WIDTH = 220;

const LINE_GRAPH_SPACE = 50;

const NULL_POINT_COLOR = "#d2d7d3";
const UNHIGHLIGHTED_POINT_COLOR = "grey";
const POINT_RADIUS = 2;

const AXIS_FONT = "normal 10px Helvetica";
const AXIS_COLOR = "#000000";

const DataWrapper = ({
  chartName,
  data,
  chartDim,
  clonotypeLabels,
  selectedClonotype,
  hoveredClonotype,
  setSelectedClonotype,
}) => {
  const { xParam, yParam, clonotypeParam } = CONSTANTS;

  const setHighlighted = (event, value) => {
    if (event === "mouseenter") {
      setSelectedClonotype({ hover: value });
    } else if (event === "mousedown") {
      setSelectedClonotype({
        hover: null,
        selected: value,
      });
    } else if (event === "mouseout") {
      setSelectedClonotype({
        hover: null,
      });
    }
  };
  return (
    <Layout title={INFO["UMAP"]["title"]} infoText={INFO["UMAP"]["text"]}>
      <UMAP
        width={chartDim["width"]}
        height={chartDim["height"]}
        chartName={chartName}
        data={data}
        highlighted={hoveredClonotype || selectedClonotype}
        xParam={xParam}
        yParam={yParam}
        subsetParam={clonotypeParam}
        subsetLabels={clonotypeLabels}
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
  subsetLabels,
  setHighlighted,
}) => {
  const [radiusRatio, setRadiusRatio] = useState(1);
  const canvasWidth = width - LEGEND_WIDTH - PADDING;
  const canvasHeight = height;

  const chartWidth =
    canvasWidth - AXIS_SPACE - PADDING - PADDING - LINE_GRAPH_SPACE;
  const chartHeight =
    canvasHeight - AXIS_SPACE - PADDING - PADDING - LINE_GRAPH_SPACE;

  const chartX = PADDING;
  const chartY = PADDING + LINE_GRAPH_SPACE;

  const yData = data.map((d) => parseFloat(d[yParam]));
  const xData = data.map((d) => parseFloat(d[xParam]));

  const yMin = Math.min(...yData);
  const yMax = Math.max(...yData);
  const xMin = Math.min(...xData);
  const xMax = Math.max(...xData);

  const xScale = d3
    .scaleLinear()
    .domain([xMin, xMax])
    .range([chartX, chartX + chartWidth]);
  const yScale = d3
    .scaleLinear()
    .domain([yMax, yMin])
    .range([chartY, chartY + chartHeight]);

  const subsetValues = subsetLabels.map(({ value }) => value);
  const subsetColors = subsetLabels.map(({ color }) => color);

  const colorScale = d3.scaleOrdinal().domain(subsetValues).range(subsetColors);

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
        subsetValues,
        colorScale,
        radiusRatio
      );

      drawLineGraph(
        context,
        data,
        subsetValues,
        xScale,
        yScale,
        xParam,
        yParam,
        subsetParam,
        colorScale,
        highlighted,
        chartX + chartWidth + 3
      );
    },
    canvasWidth,
    canvasHeight,
    [highlighted, radiusRatio]
  );

  return (
    <Grid container direction="row" style={{ padding: 0 }}>
      <Grid item>
        <canvas ref={canvasRef} />
      </Grid>
      <Grid
        container
        direction="column"
        style={{ padding: 0, width: LEGEND_WIDTH }}
      >
        <Grid item>
          <VerticalLegend
            width={LEGEND_WIDTH}
            height={chartHeight / 2}
            labels={subsetLabels}
            setHighlighted={setHighlighted}
          />
        </Grid>
        <Grid item>
          Radius Adjustment
          <Slider
            min={0}
            max={3}
            step={0.05}
            value={radiusRatio}
            disabled={highlighted}
            onChange={(event, newValue) => {
              setRadiusRatio(newValue);
            }}
          />
        </Grid>
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
  subsetLabels,
  colorScale,
  radiusRatio
) => {
  context.lineWidth = 1;
  context.globalAlpha = 1;

  const [subsetData, backgroundData] = _.partition(
    data,
    (datum) => subsetLabels.indexOf(datum[subsetParam]) !== -1
  );

  context.fillStyle = NULL_POINT_COLOR;
  backgroundData.forEach((point) => {
    context.beginPath();
    context.arc(
      xScale(point[xParam]),
      yScale(point[yParam]),
      1.5,
      0,
      Math.PI * 2,
      true
    );
    context.fill();
  });

  const yData = data.map((d) => parseFloat(d[yParam]));
  const xData = data.map((d) => parseFloat(d[xParam]));

  const yMin = Math.min(...yData);
  const yMax = Math.max(...yData);
  const xMin = Math.min(...xData);
  const xMax = Math.max(...xData);

  const groupedData = _.groupBy(subsetData, subsetParam);

  const freqData = subsetLabels.reduce((records, subset) => {
    const bins = d3Hexbin
      .hexbin()
      .x((d) => d[xParam])
      .y((d) => d[yParam])
      .radius((xMax - xMin) / 16)
      .extent([
        [xMin, xMax],
        [yMin, yMax],
      ])(groupedData[subset]);

    const freqData = bins.reduce((records, bin) => {
      const freq = bin.length;

      return [...records, ...bin.map((record) => ({ ...record, freq }))];
    }, []);

    return [...records, ...freqData];
  }, []);

  const radiusScale = d3
    .scaleLinear()
    .domain([1, Math.max(...freqData.map((record) => record["freq"]))])
    .range([POINT_RADIUS, POINT_RADIUS * 3]);

  freqData.forEach((point) => {
    context.globalAlpha = isHighlighted(point[subsetParam], highlighted)
      ? 1
      : 0.5;
    context.fillStyle = isHighlighted(point[subsetParam], highlighted)
      ? colorScale(point[subsetParam])
      : UNHIGHLIGHTED_POINT_COLOR;

    context.beginPath();
    context.arc(
      xScale(point[xParam]),
      yScale(point[yParam]),
      Math.max(POINT_RADIUS, radiusRatio * radiusScale(point["freq"])),
      0,
      Math.PI * 2,
      true
    );
    context.fill();
  });

  if (highlighted) {
    context.globalAlpha = 1;
    freqData
      .filter((datum) => datum[subsetParam] === highlighted)
      .forEach((point) => {
        context.fillStyle = isHighlighted(point[subsetParam], highlighted)
          ? colorScale(point[subsetParam])
          : UNHIGHLIGHTED_POINT_COLOR;

        context.beginPath();
        context.arc(
          xScale(point[xParam]),
          yScale(point[yParam]),
          radiusScale(point["freq"]),
          0,
          Math.PI * 2,
          true
        );
        context.fill();
      });
  }
};

const drawUMAPAxis = (context, startY, xParam, yParam) => {
  context.beginPath();
  context.font = AXIS_FONT;
  context.globalAlpha = 1;
  context.lineWidth = 1;

  const START_X = AXIS_SPACE / 2;
  const START_Y = startY + AXIS_SPACE / 2;

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

const drawLineGraph = (
  context,
  data,
  subsetLabels,
  xScale,
  yScale,
  xParam,
  yParam,
  subsetParam,
  colorScale,
  highlighted,
  startX
) => {
  const kde = (kernel, thresholds, data) =>
    thresholds.map((t) => [t, d3.mean(data, (d) => kernel(t - d))]);

  function epanechnikov(bandwidth) {
    return (x) =>
      Math.abs((x /= bandwidth)) <= 1 ? (0.75 * (1 - x * x)) / bandwidth : 0;
  }

  const lineHeightScale = d3
    .scaleLinear()
    .domain([0, 0.8])
    .range([0, LINE_GRAPH_SPACE]);

  var xLine = d3
    .line()
    .curve(d3.curveBasis)
    .x(function (d) {
      return xScale(d[0]);
    })
    .y(function (d) {
      return LINE_GRAPH_SPACE - lineHeightScale(d[1]);
    })
    .context(context);

  var yLine = d3
    .line()
    .curve(d3.curveBasis)
    .x(function (d) {
      return startX + lineHeightScale(d[1]);
    })
    .y(function (d) {
      return yScale(d[0]);
    })
    .context(context);

  const groupedData = _.groupBy(data, subsetParam);

  subsetLabels.forEach((subset) => {
    context.lineWidth = 2;
    context.globalAlpha = isHighlighted(subset, highlighted) ? 1 : 0.5;
    context.strokeStyle = isHighlighted(subset, highlighted)
      ? colorScale(subset)
      : UNHIGHLIGHTED_POINT_COLOR;

    const xDensity = kde(
      epanechnikov(1),
      xScale.ticks(100),
      groupedData[subset].map((row) => parseFloat(row[xParam]))
    );

    context.beginPath();
    xLine(xDensity);
    context.stroke();

    const yDensity = kde(
      epanechnikov(1),
      yScale.ticks(100),
      groupedData[subset].map((row) => parseFloat(row[yParam]))
    );

    context.beginPath();
    yLine(yDensity);
    context.stroke();
  });

  if (highlighted) {
    context.lineWidth = 2;
    context.globalAlpha = 1;
    context.strokeStyle = colorScale(highlighted);

    const xDensity = kde(
      epanechnikov(1),
      xScale.ticks(100),
      groupedData[highlighted].map((row) => parseFloat(row[xParam]))
    );

    context.beginPath();
    xLine(xDensity);
    context.stroke();

    const yDensity = kde(
      epanechnikov(1),
      yScale.ticks(100),
      groupedData[highlighted].map((row) => parseFloat(row[yParam]))
    );

    context.beginPath();
    yLine(yDensity);
    context.stroke();
  }
};

export default DataWrapper;
