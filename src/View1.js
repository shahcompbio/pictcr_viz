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

import Filters from "./components/Filters";
import { Heatmap, ProbabilityHistogram } from "@shahlab/planetarium";

import Grid from "@mui/material/Grid";

import { styled } from "@mui/material/styles";

import { useData } from "./provider/dataContext";

import parseClonotypeData from "./util/parseClonotypeData.js";
import parsePhenotypeData from "./util/parsePhenotypeData.js";

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
const View1 = ({ view, setView }) => {
  const { clonotypeParam, subtypeParam, logProbParam, xParam, yParam } =
    CONSTANTS;

  const [
    {
      metadata,
      filters,
      stats,
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

  const {
    clonotypeCounts,
    clonotypeLabels,
    cloneColorScale,
    clonotypeDataIsLoaded,
  } = parseClonotypeData(data, stats);

  const { phenotypeColorScale, phenotypeValues } = parsePhenotypeData(
    data,
    subset
  );

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
          <StyledGridItem id="filtersRef">
            <Filters filters={filters} />
          </StyledGridItem>
        </Grid>
        <StyledGridItem id="sunburstRef">
          <Sunburst
            data={highlightData || data}
            type={"CLONOTYPEDOUGH"}
            colors={CLONOTYPE_COLORS}
            selectedCloneColor={
              selectClone ? cloneColorScale(selectClone) : null
            }
            cloneColorScale={cloneColorScale}
            width={450}
            height={350}
            otherSubsetParam={subtypeParam}
            subsetParam={clonotypeParam}
          />
        </StyledGridItem>
        <StyledGridItem id="doughnutRef">
          <Doughnut
            data={highlightData || data}
            type={"SUBTYPEDOUGH"}
            colors={CLONOTYPE_COLORS}
            colorScale={phenotypeColorScale}
            width={450}
            height={350}
            otherSubsetParam={clonotypeParam}
            subsetParam={subset}
          />
        </StyledGridItem>
        <StyledGridItem id="cloneumapRef">
          {clonotypeDataIsLoaded && (
            <ClonotypeUMAP
              width={800}
              height={600}
              data={data}
              xParam={xParam}
              yParam={yParam}
              subsetParam={clonotypeParam}
              idParam={clonotypeParam}
              colorScale={cloneColorScale}
              labels={(value) => `Clone ${value}`}
              highlightIDs={highlightIDs}
              onLasso={(data) => {
                /*      setSelectIDs(
                  data === null ? null : data.map((datum) => datum["cell_id"])
                );
                setActiveGraph(data === null ? null : "cloneUMAP");*/
              }}
              onLegendClick={(value) => {
                /*setSelectClone(value);
                setActiveGraph(value === null ? null : "cloneUMAP");*/
              }}
              disable={activeGraph !== null && activeGraph !== "cloneUMAP"}
            />
          )}
        </StyledGridItem>
        <StyledGridItem id="phenotypeRef">
          <PhenotypeUMAP
            width={400}
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
export default View1;
