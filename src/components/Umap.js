import React, { useState, useRef, useEffect } from "react";
import * as d3 from "d3";
import * as d3Hexbin from "d3-hexbin";

import {
  VerticalLegend,
  useCanvas,
  useLasso,
  drawUMAPAxis,
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

const UMAP = ({
  width,
  height,
  data,
  xParam,
  yParam,
  subsetParam,
  idParam,
  colorScale,
  highlightIDs,
  onLasso,
  onLegendClick,
  disable,
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

  const [lassoData, drawLasso, addLassoHandler, resetLasso] = useLasso(
    data,
    xScale,
    yScale,
    xParam,
    yParam
  );

  const prevHighlightRef = useRef();

  useEffect(() => {
    // resetLasso if highlightIDs is suddenly null and lassoData still exists
    if (highlightIDs === null && prevHighlightRef.current !== null) {
      if (lassoData !== null) {
        resetLasso();
      }
    }
  }, [highlightIDs, lassoData]);

  useEffect(() => {
    prevHighlightRef.current = highlightIDs;
  }, [highlightIDs]);

  const [hoverClone, setHoverClone] = useState(null);

  const subsettedIDs =
    hoverClone !== null
      ? data
          .filter((datum) => datum[subsetParam] === hoverClone)
          .map((datum) => datum[idParam])
      : highlightIDs !== null
      ? highlightIDs
      : data.map((datum) => datum[idParam]);

  const canvasRef = useCanvas(
    (canvas) => {
      const context = canvas.getContext("2d");
      drawUMAPAxis({
        context,
        xPos: AXIS_SPACE,
        yPos: canvasHeight - AXIS_SPACE,
        xLabel: xParam,
        yLabel: yParam,
      });
      drawPoints(
        context,
        data,
        xScale,
        yScale,
        xParam,
        yParam,
        subsetParam,
        subsettedIDs,
        colorScale.domain(),
        colorScale,
        radiusRatio
      );

      drawLineGraph(
        context,
        data,
        colorScale.domain(),
        xScale,
        yScale,
        xParam,
        yParam,
        subsetParam,
        colorScale,
        hoverClone, // highlight
        chartX + chartWidth + 3
      );

      drawLasso(context);
      addLassoHandler(canvas, disable, onLasso);
    },
    canvasWidth,
    canvasHeight,
    [highlightIDs, radiusRatio, lassoData, hoverClone]
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
            ticks={colorScale
              .domain()
              .map((value) => ({ value, label: `Clone ${value}` }))}
            colorScale={colorScale}
            onHover={setHoverClone}
            title={null}
            onClick={onLegendClick}
            fontFamily={{
              regular: "MyFontLight",
              bold: "MyFontBold",
              labelOffset: 3,
            }}
            disable={disable}
          />
        </Grid>
        <Grid item>
          Radius Adjustment
          <Slider
            min={0}
            max={3}
            step={0.05}
            value={radiusRatio}
            disabled={highlightIDs !== null}
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

  const [subsetData, backgroundData] = _.partition(data, (datum) =>
    subsetLabels.includes(datum[subsetParam])
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
    context.globalAlpha =
      highlighted === null || highlighted.includes(point["cell_id"]) ? 1 : 0.5;
    context.fillStyle =
      highlighted === null || highlighted.includes(point["cell_id"])
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
      .filter((datum) => highlighted.includes(datum["cell_id"]))
      .forEach((point) => {
        context.fillStyle = colorScale(point[subsetParam]);

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

export default UMAP;
