import React, { useEffect, useState } from "react";
import _ from "lodash";
import { INFO } from "../config.js";

import {
  useCanvas,
  useLasso,
  VerticalLegend,
  Layout,
  drawUMAPAxis,
  Select,
  //  WebglUMAP,
} from "@shahlab/planetarium";
import ColorSelect from "./VertSelect";
import ReglUmap from "./util/ReglUmap";
import { useGL } from "./util/useGL";

import { CircularProgress, Grid } from "@mui/material";
import { drawAxis, getLassoObj } from "./util/util";
import TtestResults from "./TtestResults";
import * as d3 from "d3";
import { useData } from "../provider/dataContext";

const PADDING = 10;
const AXIS_SPACE = 20;
const NUM_LEGEND_WIDTH = 250;
const AXIS_LENGTH = 50;
const AXIS_FONT = "Helvetica";
const AXIS_COLOR = "#000000";
const canvasWidth = 900;
const canvasHeight = 600;
const COLOR_ARRAY = [
  "#5E4FA2",
  "#3288BD",
  "#66C2A5",
  "#FEE08B",
  "#FDAE61",
  "#F46D43",
  "#D53E4F",
  "#c9cc76",
];

const layerNames = ["umapCanvas"];
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

const PhenotypeUMAP = ({
  data,
  xParam,
  yParam,
  idParam,
  colorScale,
  onLasso,
  onLegendClick,
  disable,
  highlightIDs,

  loadingTest,
  labels = (value) => value,
}) => {
  const [canvas, setCanvas] = useState(null);
  const [{ filters, subset }, dispatch] = useData();
  const [subsetParam, setSubset] = useState(subset);

  const isCategorical =
    typeof data.filter((datum) => datum.hasOwnProperty(subsetParam))[0][
      subsetParam
    ] === "string";

  const legendFilter = isCategorical
    ? (value, datum) => datum[subsetParam] === value
    : (value, datum) =>
        datum.hasOwnProperty(subsetParam) && datum[subsetParam] >= value[0];

  const legendColors = getColorScale({
    data,
    subsetParam,
    isCategorical: isCategorical,
  });

  const legendWidth = NUM_LEGEND_WIDTH;

  const chartWidth = canvasWidth - AXIS_SPACE - legendWidth;
  const chartHeight = canvasHeight - AXIS_SPACE - PADDING - PADDING;

  const [wrapperRef] = useGL(chartWidth, chartHeight, layerNames);

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

  const [hoveredLegend, setHoveredLegend] = useState(null);
  const [clickedLegend, setClickedLegend] = useState(null);

  const legendIDs = hoveredLegend || clickedLegend;

  useEffect(() => {
    const layers = layerNames.map((layer) => {
      const canvasLayer = d3.select("#" + layer);
      canvasLayer.attr("position", "absolute");
      return canvasLayer.node();
    });
    setCanvas(layers[0]);
  }, [wrapperRef]);

  const [lassoData, drawLasso, addLassoHandler, resetLasso] = useLasso(
    data,
    xScale,
    yScale,
    xParam,
    yParam
  );

  const subsettedIDs = [];
  const canvasRef = useCanvas(
    (canvas) => {
      const context = canvas.getContext("2d");
      canvas.id = "lasso";

      drawUMAPAxis({
        context,
        xPos: AXIS_SPACE / 2,
        yPos: canvasHeight - AXIS_SPACE / 2,
        xLabel: xParam,
        yLabel: yParam,
      });

      drawLasso(context);
      const disableLasso = disable || legendIDs !== null;
      addLassoHandler(canvas, disableLasso, onLasso);
    },
    chartWidth,
    canvasHeight,
    [data, disable, subsetParam, subsettedIDs]
  );

  const getLegendData = (value) => {
    return value === null
      ? value
      : data
          .filter((datum) => legendFilter(value, datum))
          .map((datum) => datum[idParam]);
  };
  console.log(data);
  return (
    <div style={{ width: canvasWidth }}>
      <Grid
        container
        direction="row"
        style={{
          padding: 0,
          position: "relative",
          height: canvasHeight,
          width: canvasWidth,
        }}
      >
        <Grid
          item
          style={{
            float: "left",
            position: "absolute",
            left: 0,
          }}
        >
          <span style={{ marginRight: 10 }}>
            <ColorSelect
              id={"phenotypeSelect"}
              width={legendWidth}
              title={"Colored By"}
              value={subset}
              options={filters.map((datum) => datum["name"])}
              onSelect={setSubset}
            />
          </span>
          <VerticalLegend
            width={legendWidth}
            height={canvasHeight / 2}
            colorScale={legendColors}
            ticks={
              isCategorical
                ? legendColors
                    .domain()
                    .sort()
                    .map((value) => ({ value, label: value }))
                : 10
            }
            onHover={(value) => {
              const legendData = getLegendData(value);
              setHoveredLegend(legendData);
              //  onLegendHover(value);
            }}
            onClick={(value) => {
              const legendData = getLegendData(value);
              setClickedLegend(legendData);
              //  onLegendClick(value);
            }}
            disable={disable || lassoData !== null}
            reset={highlightIDs !== null}
          />
          <TtestResults
            resetLasso={resetLasso}
            idParam={idParam}
            cells={lassoData}
            height={canvasHeight / 2}
            width={NUM_LEGEND_WIDTH}
            count={lassoData ? lassoData.length : null}
          />
        </Grid>

        <Grid
          item
          style={{
            position: "relative",
            marginLeft: legendWidth,
            width: chartWidth,
          }}
          ref={wrapperRef}
        >
          {canvas && data && data.length > 0 && (
            <ReglUmap
              canvasRef={canvas}
              pointSize={4}
              data={data}
              lassoDataObj={
                lassoData && lassoData.length > 0
                  ? getLassoObj(lassoData, idParam)
                  : {}
              }
              isCategorical={isCategorical}
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
          <canvas
            ref={canvasRef}
            style={{ zIndex: 100, position: "absolute" }}
          />
        </Grid>
      </Grid>
    </div>
  );
};

export default PhenotypeUMAP;
