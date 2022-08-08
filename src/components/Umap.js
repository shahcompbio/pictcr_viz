import React, { useState, useRef, useEffect } from "react";
import * as d3 from "d3";
import * as d3Hexbin from "d3-hexbin";

import {
  VerticalLegend,
  useCanvas,
  useLasso,
  useGL,
  drawUMAPAxis,
  isHighlighted,
} from "@shahlab/planetarium";
import ReglUmap from "./util/ReglUmap";
import Grid from "@mui/material/Grid";
import Slider from "@mui/material/Slider";
import Title from "./Title.js";

import { drawAxis, getLassoObj } from "./util/util";

import { useData } from "../provider/dataContext";

import { CONSTANTS } from "../config";

import _ from "lodash";

const PADDING = 10;
const AXIS_SPACE = 20;
const LEGEND_WIDTH = 320;

const LINE_GRAPH_SPACE = 50;

const NULL_POINT_COLOR = "#d2d7d3";
const UNHIGHLIGHTED_POINT_COLOR = "grey";
const POINT_RADIUS = 2;

const GREY_VEC4 = "[0.810, 0.786, 0.786, 1.0]";

const COLOR_ARRAY = [
  "[0.369,0.31,0.635,1.0]",
  "[0.196,0.533,0.741,1.0]",
  "[0.4,0.761,0.647,1.0]",
  "[0.996,0.878,0.545,1.0]",
  "[0.957,0.427,0.263,1.0]",
  "[0.835,0.243,0.31,1.0]",
  "[0.788,0.8,0.463,1.0]",
  "[0.62,0.004,0.259,1.0]",
  "[0.776,0.682,1.0,1.0]",
  "[0.741,0.847,1.0,1.0]",
  "[0.741,1.0,0.698,1.0]",
  "[1.0,0.784,0.682,1.0]",
  "[1.0,0.624,0.733,1.0]",
  "[0.698,0.859,0.839,1.0]",
  "[1.0,0.831,0.439,1.0]",
];

const layerNames = ["umapPhenotypeCanvas"];
const getColorScale = ({ data, subsetParam, isCategorical }) => {
  if (isCategorical) {
    const subsetGroups = _.groupBy(data, subsetParam);
    const subsetValues = Object.keys(subsetGroups).sort();

    return d3
      .scaleOrdinal()
      .domain(subsetValues)
      .range(
        COLOR_ARRAY.slice(0, Math.min(subsetValues.length, COLOR_ARRAY.length))
      );
  } else {
    const subsetData = data
      .filter((d) => d.hasOwnProperty(subsetParam))
      .map((d) => parseFloat(d[subsetParam]));

    const subsetMax = Math.max(...subsetData);
    return d3
      .scaleSequential(d3.interpolateViridis)
      .domain([0, subsetMax])
      .nice();
  }
};
const UMAP = ({
  width,
  height,
  data,
  xParam,
  yParam,
  subsetParam,
  idParam,
  colorScale,
  cloneHEXColorScale,
  highlightIDs,
  onLasso,
  onLegendClick,
  disable,
}) => {
  const [{ metadata }, dispatch] = useData();
  const [canvas, setCanvas] = useState(null);
  const [radiusRatio, setRadiusRatio] = useState(2);
  const canvasWidth = width - LEGEND_WIDTH - PADDING;
  const canvasHeight = height;

  const chartWidth =
    canvasWidth - AXIS_SPACE - PADDING - PADDING - LINE_GRAPH_SPACE;
  const chartHeight =
    canvasHeight - AXIS_SPACE - PADDING - PADDING - LINE_GRAPH_SPACE;

  const [wrapperRef] = useGL(chartWidth, chartHeight, layerNames);

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

  useEffect(() => {
    const layers = layerNames.map((layer) => {
      const canvasLayer = d3.select("#" + layer);
      canvasLayer.attr("position", "absolute");
      return canvasLayer.node();
    });
    setCanvas(layers[0]);
  }, [wrapperRef]);

  const canvasRef = useCanvas(
    (canvas) => {
      const context = canvas.getContext("2d");
      drawUMAPAxis({
        context,
        xPos: AXIS_SPACE / 2,
        yPos: chartHeight + AXIS_SPACE,
        xLabel: xParam,
        yLabel: yParam,
      });

      drawLineGraph(
        context,
        data,
        colorScale.domain(),
        xScale,
        yScale,
        xParam,
        yParam,
        subsetParam,
        cloneHEXColorScale,
        hoverClone, // highlight
        chartX + chartWidth + 3
      );

      drawLasso(context);
      addLassoHandler(canvas, disable, onLasso);
    },
    canvasWidth,
    canvasHeight,
    [highlightIDs, radiusRatio, lassoData, hoverClone, data]
  );
  const getDataWithAttributes = (
    metadata,
    data,
    colorScale,
    subsetParam,
    xParam,
    yParam,
    highlighted,
    radiusRatio,
    hoverClone
  ) => {
    console.log(hoverClone);
    console.log(allLabels);
    console.log(metadata);
    const allLabels = colorScale.domain();
    const [subsetData, backgroundData] = _.partition(metadata, (datum) =>
      hoverClone
        ? allLabels.includes(datum[subsetParam]) &&
          hoverClone.includes(datum[subsetParam])
        : allLabels.includes(datum[subsetParam])
    );
    const subsetLabels = _.uniq(subsetData.map((d) => d[subsetParam]));

    const backgroundDataWithAttr = backgroundData.map((d) => ({
      ...d,
      color: GREY_VEC4,
      pointSize: radiusRatio,
    }));

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

    const freqDataWithAttr = freqData.map((point) => {
      const opacity =
        highlighted === null || highlighted.includes(point["cell_id"]);

      const colour =
        highlighted === null || highlighted.includes(point["cell_id"])
          ? colorScale(point[subsetParam])
          : UNHIGHLIGHTED_POINT_COLOR;

      return { ...point, color: colour, pointSize: radiusScale(point["freq"]) };
    });

    return [...backgroundDataWithAttr, ...freqDataWithAttr];
  };

  const modifiedData = getDataWithAttributes(
    metadata,
    data,
    colorScale,
    subsetParam,
    xParam,
    yParam,
    null,
    radiusRatio,
    hoverClone
  );
  return (
    <Grid
      container
      direction="row"
      style={{
        padding: 0,
        position: "relative",
        height: height,
        width: width,
      }}
    >
      <Grid
        container
        direction="column"
        style={{ padding: 0, width: LEGEND_WIDTH, paddingTop: "30px" }}
      >
        <Title title="Clonotype UMAP" />
        <Grid item>
          <VerticalLegend
            width={LEGEND_WIDTH}
            height={chartHeight / 2}
            ticks={colorScale
              .domain()
              .map((value) => ({ value, label: `Clone ${value}` }))}
            colorScale={cloneHEXColorScale}
            onHover={setHoverClone}
            title={null}
            onClick={onLegendClick}
            fontFamily={{
              regular: "Noto Sans",
              bold: "Noto Sans",
              fontSize: "20px",
              labelOffset: 1,
            }}
            squareSize={15}
            disable={disable}
          />
        </Grid>
        <Grid item sx={{ width: "80%" }}>
          Radius Adjustment
          <Slider
            min={0}
            max={10}
            step={0.05}
            sx={{
              color: "black",
            }}
            value={radiusRatio}
            disabled={highlightIDs !== null}
            onChange={(event, newValue) => {
              setRadiusRatio(newValue);
            }}
          />
        </Grid>
      </Grid>
      <Grid
        item
        style={{ position: "", padding: 0, width: chartWidth + 50 }}
        ref={wrapperRef}
      >
        {canvas && modifiedData && modifiedData.length > 0 && (
          <ReglUmap
            canvasRef={canvas}
            pointSize={radiusRatio}
            data={modifiedData}
            lassoDataObj={
              lassoData && lassoData.length > 0
                ? getLassoObj(lassoData, idParam)
                : {}
            }
            isCategorical={false}
            width={chartWidth}
            height={chartHeight}
            xParam={xParam}
            yParam={yParam}
            xScale={xScale}
            yScale={yScale}
            idParam={idParam}
            subsetParam={subsetParam}
          />
        )}
        <canvas ref={canvasRef} style={{ zIndex: 100, position: "absolute" }} />
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
  colorHexScale,
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

  subsetLabels.forEach((subset, i) => {
    if (i !== subsetLabels.length - 1) {
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
    } else {
      context.lineWidth = 2;
      context.globalAlpha = isHighlighted(subset, highlighted) ? 1 : 0.5;

      context.strokeStyle = isHighlighted(subset, highlighted)
        ? colorScale(subset)
        : UNHIGHLIGHTED_POINT_COLOR;
      context.filleStyle = "rgba(255, 255, 255, 0.01)";

      const xDensity = kde(
        epanechnikov(1),
        xScale.ticks(100),
        groupedData[subset].map((row) => parseFloat(row[xParam]))
      );
      context.closePath();
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
    }
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
    context.closePath();

    const yDensity = kde(
      epanechnikov(1),
      yScale.ticks(100),
      groupedData[highlighted].map((row) => parseFloat(row[yParam]))
    );

    context.beginPath();
    yLine(yDensity);
    context.stroke();
    context.closePath();
  }
};

export default UMAP;
