import React, { useRef } from "react";
import _ from "lodash";

import ClonotypeUMAP from "./components/Umap";
import ClonotypeExpansion from "./components/ClonotypeExpansion";
import RankedOrder from "./components/RankedOrder";
import Doughnut from "./components/Doughnut";
import Sunburst from "./components/Sunburst";
import PhenotypeUMAP from "./components/WebglUMAP.js";
import ScrollBar from "./components/ScrollBar";
import ViewButtons from "./components/ViewButtons";

import { Heatmap, ProbabilityHistogram } from "@shahlab/planetarium";

import Grid from "@mui/material/Grid";

import { styled } from "@mui/material/styles";

import { useData } from "./provider/dataContext";

import parseClonotypeData from "./util/parseClonotypeData.js";
import parsePhenotypeData from "./util/parsePhenotypeData.js";

import { CONSTANTS, CLONOTYPE_COLORS, CLONOTYPE_COLORS_2 } from "./config";

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
const parseData = (metadata, selectFilters) =>
  selectFilters === null
    ? metadata
    : metadata.filter((datum) => datum[selectFilters[0]] === selectFilters[1]);

const View1 = ({ view, setView }) => {
  const { clonotypeParam, subtypeParam, logProbParam, xParam, yParam } =
    CONSTANTS;

  const [
    {
      metadata,
      filters,
      settings,
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
  ] = useData();

  const inputRef = useRef([]);

  const data = parseData(metadata, selectFilters);
  console.log(data);
  console.log(settings);
  const statsData =
    selectFilters === null
      ? stats
      : data.reduce((final, d) => {
          const clone = d[clonotypeParam];
          final[clone] = final[clone] ? final[clone] + 1 : 1;
          return final;
        }, {});

  const {
    clonotypeCounts,
    clonotypeLabels,
    cloneColorScale,
    cloneHEXColorScale,
    clonotypeDataIsLoaded,
  } = parseClonotypeData(data, statsData, settings.filterNA);
  console.log(clonotypeLabels);
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
        <StyledGridItem id="sunburstRef">
          <Sunburst
            data={highlightData || data}
            type={"CLONOTYPEDOUGH"}
            colors={CLONOTYPE_COLORS_2}
            selectedCloneColor={
              selectClone ? cloneColorScale(selectClone) : null
            }
            cloneHEXColorScale={cloneHEXColorScale}
            cloneColorScale={cloneColorScale}
            width={450}
            height={350}
            otherSubsetParam={subtypeParam}
            subsetParam={clonotypeParam}
          />
        </StyledGridItem>
        <StyledGridItem id="cloneumapRef">
          {clonotypeDataIsLoaded && (
            <ClonotypeUMAP
              width={1000}
              height={700}
              data={data}
              xParam={xParam}
              yParam={yParam}
              subsetParam={clonotypeParam}
              idParam={clonotypeParam}
              cloneHEXColorScale={cloneHEXColorScale}
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
      </Grid>
    </span>
  );
};

export default View1;
