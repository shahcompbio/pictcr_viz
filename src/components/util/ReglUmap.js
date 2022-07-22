import React from "react";

import * as d3 from "d3";
import _ from "lodash";

import { Grid } from "@mui/material";
import { ReglFrame } from "react-regl";

import WebglUMAP from "./drawWebglUmap";
import { reglClear } from "./util";

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
const GREY_VEC4 = "[0.810, 0.786, 0.786, 1.0]";

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
      .scaleSequential(
        COLOR_ARRAY.slice(0, Math.min(subsetData.length, COLOR_ARRAY.length))
      )
      .domain([0, subsetMax])
      .nice();
  }
};

const isLassoSelected = (lassoDataObj) =>
  lassoDataObj && Object.keys(lassoDataObj).length > 0;

const ReglUmap = ({
  data,
  width = 500,
  height = 500,
  xParam,
  yParam,
  subsetParam,
  idParam = "id",
  yScale,
  xScale,
  canvasRef,
  lassoDataObj = {},
  isCategorical = true,
  pointSize = 2,
}) => {
  const subsetColors = getColorScale({
    data,
    subsetParam,
    isCategorical: isCategorical,
  });

  const hasLassoSelected = isLassoSelected(lassoDataObj);

  const newData = data.map((d) => {
    const x = xScale(d[xParam]);
    const y = yScale(d[yParam]);

    const isGrey = hasLassoSelected && !lassoDataObj.hasOwnProperty(d[idParam]);
    //get rid of this
    const color = d["color"] ? d["color"] : subsetColors(d[subsetParam]);
    //      : subsetColors(d[subsetParam]);

    return {
      ...d,
      x: x,
      y: y,
      color: color,
      pointSize: d["pointSize"],
    };
  });
  return newData.length > 0 ? (
    <ReglFrame
      canvasRef={canvasRef}
      onMount={reglClear}
      onFrame={(context, regl) => regl.clear({ color: [1, 1, 1, 0] })}
    >
      <WebglUMAP
        data={newData}
        stageWidth={width}
        stageHeight={height}
        pointSize={pointSize}
      />
    </ReglFrame>
  ) : null;
};

export default ReglUmap;
