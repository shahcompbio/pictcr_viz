import React, { useRef } from "react";
import _ from "lodash";
import * as d3 from "d3";

import ClonotypeUMAP from "./components/Umap";
import ClonotypeExpansion from "./components/ClonotypeExpansion";
import RankedOrder from "./components/RankedOrder";
import Doughnut from "./components/Doughnut";
import Sunburst from "./components/Sunburst";

import PhenotypeUMAP from "./components/WebglUMAP.js";
import ScrollBar from "./components/ScrollBar";
import ViewButtons from "./components/ViewButtons";

import { Heatmap, ProbabilityHistogram, Sankey } from "@shahlab/planetarium";

import Grid from "@mui/material/Grid";

import { styled } from "@mui/material/styles";

import { useData } from "./provider/dataContext";

import { CONSTANTS, CLONOTYPE_COLORS } from "./config";

const PHENOTYPE_COLORS = [
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
const View2 = ({ metadata, filters, view, setView }) => {
  const { clonotypeParam, subtypeParam, logProbParam, xParam, yParam } =
    CONSTANTS;

  const [
    {
      inputRefMapping,
      selectPhenotype,
      selectClone,
      selectIDs,
      activeGraph,
      selectFilters,
      subset,
      selectSubset,
    },
    dispatch,
  ] = useData();
  const inputRef = useRef([]);

  const data =
    selectFilters === null
      ? metadata
      : metadata.filter(
          (datum) => datum[selectFilters[0]] === selectFilters[1]
        );
  // Remove none
  const clonotypeCounts = _.countBy(
    data.filter((datum) => datum[clonotypeParam] !== "None"),
    clonotypeParam
  );

  const clonotypeLabels = Object.keys(clonotypeCounts)
    .sort((a, b) => clonotypeCounts[b] - clonotypeCounts[a])
    .slice(0, 10)
    .map((value, index) => ({
      value,
      label: `Clone ${value}`,
      color: CLONOTYPE_COLORS[index],
    }));

  const cloneColorScale = d3
    .scaleOrdinal()
    .domain(clonotypeLabels.map((label) => label["value"]))
    .range(
      CLONOTYPE_COLORS.slice(
        0,
        Math.min(clonotypeLabels.length, CLONOTYPE_COLORS.length)
      )
    )
    .unknown("#e8e8e8");

  const phenotypeValues = Object.keys(_.groupBy(data, subset)).sort();
  const phenotypeColorScale = d3
    .scaleOrdinal()
    .domain(phenotypeValues)
    .range(
      PHENOTYPE_COLORS.slice(
        0,
        Math.min(phenotypeValues.length, PHENOTYPE_COLORS.length)
      )
    )
    .unknown("#e8e8e8");

  const highlightData =
    selectIDs !== null
      ? data.filter((datum) => selectIDs.includes(datum["cell_id"]))
      : selectClone !== null
      ? data.filter((datum) => datum[clonotypeParam] === selectClone)
      : selectSubset !== null
      ? data.filter((datum) => datum[subset] === selectSubset)
      : null;

  const highlightIDs =
    highlightData === null
      ? null
      : highlightData.map((datum) => datum["cell_id"]);

  const probabilities = data.filter(
    (datum) =>
      datum[clonotypeParam] !== "None" || datum[logProbParam] !== "None"
  );

  const subtypeTotals = _.countBy(data, subtypeParam);

  const StyledGridItem = styled((props) => (
    <Grid
      item
      sx={{ margin: "auto", marginRight: "50px" }}
      ref={(el) => (inputRef.current[props.id] = el)}
    >
      {props.children}
    </Grid>
  ))(({ theme }) => ({
    margin: "auto",
  }));

  return (
    <span>
      <Grid
        container
        direction="row"
        justifyContent="flex-start"
        alignItems="flex-start"
        spacing={3}
        style={{
          overflowX: "scroll !important",
          flexWrap: "nowrap",
          padding: 15,
          marginBottom: 10,
          backgroundColor: "white",
          marginTop: 50,
        }}
      >
        <Grid
          item
          container
          direction="column"
          justifyContent="flex-start"
          alignItems="flex-start"
          id="filters-grid"
          style={{ width: 400 }}
        >
          <ViewButtons view={view} setView={setView} />
          <StyledGridItem id="filtersRef"></StyledGridItem>
        </Grid>
        <StyledGridItem id="sankeyRef">
          <Sankey
            width={700}
            height={400}
            data={data}
            subsetParam={subset}
            cloneParam={"clone_id"}
            timepointOrder={["Pre", "Post"]}
            timepointParam={"treatment"}
          />
        </StyledGridItem>
        <StyledGridItem id="phenotypeRef">
          <PhenotypeUMAP
            width={700}
            height={600}
            data={data}
            xParam={xParam}
            yParam={yParam}
            subsetParam={subset}
            idParam="cell_id"
            colorScale={phenotypeColorScale}
            onLasso={(data) => {
              /*  setSelectIDs(
                data === null ? null : data.map((datum) => datum["cell_id"])
              );
              setActiveGraph(data === null ? null : "phenoUMAP");*/
            }}
            onLegendClick={(value) => {
              /*  dispatch({
                type: "setSelectSubset",
                value: value,
              });
              dispatch({
                type: "setActiveGraph",
                value: value === null ? null : "phenoUMAP",
              });*/
            }}
            disable={activeGraph !== null && activeGraph !== "phenoUMAP"}
            highlightIDs={highlightIDs}
          />
        </StyledGridItem>
        <StyledGridItem id="heatmapRef">
          <Heatmap
            width={750}
            height={550}
            font={"Helvetica"}
            data={probabilities}
            column={clonotypeParam}
            row={subtypeParam}
            highlightedRow={selectPhenotype}
            highlightedColumn={selectClone}
            columnLabels={clonotypeLabels}
            rowTotal={subtypeTotals}
          />
        </StyledGridItem>
        <StyledGridItem id="clonotypeExpansionRef">
          <ClonotypeExpansion
            chartName={"BARPLOT"}
            data={probabilities}
            width={750}
            height={550}
            highlightedRow={selectPhenotype}
          />
        </StyledGridItem>
        <StyledGridItem id="probabilityHistogramRef">
          <ProbabilityHistogram
            data={probabilities}
            font={"Helvetica"}
            width={750}
            height={500}
            probParam={logProbParam}
            subgroupParam={subtypeParam}
            observationParam={clonotypeParam}
            highlightedObservation={selectClone}
            highlightedSubgroup={selectPhenotype}
          />
        </StyledGridItem>
        <StyledGridItem id="rankedOrderRef">
          <RankedOrder
            width={700}
            height={500}
            data={probabilities}
            highlight={selectPhenotype}
          />
        </StyledGridItem>
      </Grid>
      <ScrollBar
        scrollBarWidth={1600}
        height={50}
        inputRef={inputRef}
        inputRefMapping={inputRefMapping[view]}
      />
    </span>
  );
};
export default View2;
