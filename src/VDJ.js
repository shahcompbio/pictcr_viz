import React, { useState, useEffect } from "react";
import _ from "lodash";
import * as d3 from "d3";

import ClonotypeUMAP from "./components/Umap";
import PhenotypeUMAP from "./components/subTypeUmap";
import ClonotypeExpansion from "./components/ClonotypeExpansion";
import DEGTable from "./components/DEGTable";
import RankedOrder from "./components/RankedOrder";
import Doughnut from "./components/Doughnut";
import MetaData from "./components/MetaData";
import Header from "./components/Header";
import Sunburst from "./components/Sunburst";

import Paper from "@mui/material/Paper";
import Filters from "./components/Filters";

import {
  Heatmap,
  ProbabilityHistogram,
  Layout,
  Select,
} from "@shahlab/planetarium";

import Grid from "@mui/material/Grid";

import { theme } from "./theme/theme";
import { ThemeProvider, StyledEngineProvider } from "@mui/material/styles";
import CssBaseline from "@mui/material/CssBaseline";

import { CONSTANTS, CLONOTYPE_COLORS, INFO } from "./config";

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
const DataWrapper = ({ data }) => (
  <VDJ
    metadata={data["metadata"]}
    filters={data["filters"]}
    degs={data["degs"]}
  />
);

export const VDJ = ({ metadata, degs, filters }) => {
  const { clonotypeParam, subtypeParam, logProbParam, xParam, yParam } =
    CONSTANTS;

  const [selectPhenotype, setSelectPhenotype] = useState(null);
  const [selectClone, setSelectClone] = useState(null);
  const [selectIDs, setSelectIDs] = useState(null);
  const [activeGraph, setActiveGraph] = useState(null);

  const [subset, setSubset] = useState(subtypeParam);
  const [selectSubset, setSelectSubset] = useState(null);

  const [selectFilters, setSelectFilters] = useState(null);

  const [tTestData, settTestData] = useState([]);
  const data =
    selectFilters === null
      ? metadata
      : metadata.filter(
          (datum) => datum[selectFilters[0]] === selectFilters[1]
        );

  useEffect(() => {
    if (selectIDs !== null) {
      if (selectIDs.length !== 0) {
        const param = selectIDs.join(",");
        fetch("http://localhost:5000/testing/" + param + "/", {
          credentials: "include",
        })
          .then((res) => res.json())
          .then((result) => {
            if (result.data) {
              settTestData(result.data);
            }
          });

        /*  fetch("http://localhost:5000/isLoaded/", {
          credentials: "include",
        })
          .then((res) => res.json())
          .then((result) => {
            if (!result.data) {
              console.log("loading");

              fetch(
                "http://localhost:5000/l/Users/vbojilova/Projects/pictcr_viz/src/data/hacohen_viz.h5ad/",
                {
                  credentials: "include",
                }
              )
                .then((res) => res.json())
                .then((result) => {
                  fetch("http://localhost:5000/testing/" + param + "/", {
                    credentials: "include",
                  })
                    .then((res) => res.json())
                    .then((result) => {
                      if (result.data) {
                        settTestData(result.data);
                      }
                    });
                });
            } else {
              fetch("http://localhost:5000/testing/" + param + "/", {
                credentials: "include",
              })
                .then((res) => res.json())
                .then((result) => {
                  if (result.data) {
                    settTestData(result.data);
                  }
                });
            }
            console.log(result);
          });
      } else {
        settTestData([]);
      }*/
      }
    }
  }, [selectIDs]);

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

  return (
    <StyledEngineProvider injectFirst>
      <ThemeProvider theme={theme}>
        <CssBaseline />
        <Grid
          container
          direction="column"
          justifyContent="flex-start"
          alignItems="flex-start"
          style={{ xOverflow: "none" }}
        >
          <Grid
            container
            direction="row"
            justifyContent="flex-start"
            alignItems="flex-start"
            spacing={3}
          >
            <Grid item xs={3} style={{ width: 200 }}>
              <Paper
                elevation={0}
                style={{
                  background: "none",
                  margin: 10,
                  padding: 10,
                }}
              >
                <MetaData
                  width={250}
                  data={data}
                  sample="Hacohen"
                  filters={filters}
                  highlighted={highlightData}
                  selected={selectClone || selectSubset}
                  setHighlight={() => {
                    setSelectIDs(null);
                    setSelectClone(null);
                    setSelectSubset(null);
                    setActiveGraph(null);
                    setSelectFilters(null);
                  }}
                  selectFilters={selectFilters}
                  selectedType={selectClone ? "Clone" : selectSubset}
                  setFilters={setSelectFilters}
                  totalCount={metadata.length}
                />
                <Filters
                  selected={selectFilters}
                  filters={filters}
                  setFilters={setSelectFilters}
                />
              </Paper>
            </Grid>
            <Grid item>
              <Grid
                container
                direction="column"
                justifyContent="flex-start"
                alignItems="flex-start"
              >
                <Grid item>
                  <Header />
                </Grid>
                <Grid item>
                  <Grid
                    item
                    container
                    direction="row"
                    justify="flex-start"
                    alignItems="flex-start"
                  >
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
                    <Doughnut
                      data={highlightData || data}
                      type={"SUBTYPEDOUGH"}
                      colors={CLONOTYPE_COLORS}
                      colorScale={phenotypeColorScale}
                      width={450}
                      height={350}
                      otherSubsetParam={clonotypeParam}
                      subsetParam={subset}
                      Select={
                        <span style={{ marginRight: 10 }}>
                          <Select
                            width={200}
                            title={"Color By"}
                            value={subset}
                            options={filters.map((datum) => datum["name"])}
                            onSelect={setSubset}
                          />
                        </span>
                      }
                    />
                  </Grid>
                </Grid>
              </Grid>
            </Grid>
          </Grid>
          <Grid
            item
            container
            direction="row"
            justifyContent="flex-start"
            alignItems="flex-start"
          >
            <Layout
              title={INFO["UMAP"]["title"]}
              infoText={INFO["UMAP"]["text"]}
            >
              <ClonotypeUMAP
                width={800}
                height={600}
                data={data}
                xParam={xParam}
                yParam={yParam}
                subsetParam={clonotypeParam}
                idParam="cell_id"
                colorScale={cloneColorScale}
                labels={(value) => `Clone ${value}`}
                highlightIDs={highlightIDs}
                onLasso={(data) => {
                  setSelectIDs(
                    data === null ? null : data.map((datum) => datum["cell_id"])
                  );
                  setActiveGraph(data === null ? null : "cloneUMAP");
                }}
                onLegendClick={(value) => {
                  setSelectClone(value);
                  setActiveGraph(value === null ? null : "cloneUMAP");
                }}
                disable={activeGraph !== null && activeGraph !== "cloneUMAP"}
              />
            </Layout>
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
                setSelectIDs(
                  data === null ? null : data.map((datum) => datum["cell_id"])
                );
                setActiveGraph(data === null ? null : "phenoUMAP");
              }}
              onLegendClick={(value) => {
                setSelectSubset(value);
                setActiveGraph(value === null ? null : "phenoUMAP");
              }}
              disable={activeGraph !== null && activeGraph !== "phenoUMAP"}
              highlightIDs={highlightIDs}
              Select={
                <span style={{ marginRight: 10 }}>
                  <Select
                    width={200}
                    title={"Color By"}
                    value={subset}
                    options={filters.map((datum) => datum["name"])}
                    onSelect={setSubset}
                  />
                </span>
              }
              tTestData={tTestData}
            />
          </Grid>
          <Grid
            item
            container
            direction="row"
            justify="flex-start"
            alignItems="flex-end"
          >
            <Layout
              title={INFO["HEATMAP"]["title"]}
              infoText={INFO["HEATMAP"]["text"]}
            >
              <Heatmap
                width={750}
                height={550}
                font={"MyFontLight"}
                data={probabilities}
                column={clonotypeParam}
                row={subtypeParam}
                highlightedRow={selectPhenotype}
                highlightedColumn={selectClone}
                columnLabels={clonotypeLabels}
                rowTotal={subtypeTotals}
              />
            </Layout>
            <ClonotypeExpansion
              chartName={"BARPLOT"}
              data={probabilities}
              width={750}
              height={550}
              highlightedRow={selectPhenotype}
            />
          </Grid>
          <Grid
            item
            container
            direction="row"
            justify="flex-start"
            alignItems="flex-start"
          >
            <Layout
              title={INFO["HISTOGRAM"]["title"]}
              infoText={INFO["HISTOGRAM"]["text"]}
            >
              <ProbabilityHistogram
                data={probabilities}
                width={750}
                height={500}
                probParam={logProbParam}
                subgroupParam={subtypeParam}
                observationParam={clonotypeParam}
                highlightedObservation={selectClone}
                highlightedSubgroup={selectPhenotype}
              />
            </Layout>
            <Layout
              title={INFO["RANKED"]["title"]}
              infoText={INFO["RANKED"]["text"]}
            >
              <RankedOrder
                width={800}
                height={500}
                data={probabilities}
                highlight={selectPhenotype}
              />
            </Layout>
          </Grid>
          {/* <Grid
            item
            container
            direction="row"
            justify="flex-start"
            alignItems="flex-end"
          >
            <Layout
              title={INFO["HEATMAP"]["title"]}
              infoText={INFO["HEATMAP"]["text"]}
            >
              <Heatmap
                width={750}
                height={550}
                font={"MyFontLight"}
                data={probabilities}
                column={clonotypeParam}
                row={subtypeParam}
                highlightedRow={selectPhenotype}
                highlightedColumn={selectClone}
                columnLabels={clonotypeLabels}
                rowTotal={subtypeTotals}
              />
            </Layout>
            <ClonotypeExpansion
              chartName={"BARPLOT"}
              data={probabilities}
              width={750}
              height={550}
              highlightedRow={selectPhenotype}
            />
          </Grid>
          <Grid
            item
            container
            direction="row"
            justify="flex-start"
            alignItems="flex-start"
          >
            <Layout
              title={INFO["HISTOGRAM"]["title"]}
              infoText={INFO["HISTOGRAM"]["text"]}
            >
              <ProbabilityHistogram
                data={probabilities}
                width={750}
                height={500}
                probParam={logProbParam}
                subgroupParam={subtypeParam}
                observationParam={clonotypeParam}
                highlightedObservation={selectClone}
                highlightedSubgroup={selectPhenotype}
              />
            </Layout>
            <DEGTable
              chartName={"TABLE"}
              data={degs}
              selectedSubtype={selectPhenotype || selectClone}
              chartDim={{
                height: 500,
                width: 750,
              }}
            />
          </Grid>
          <Layout
            title={INFO["RANKED"]["title"]}
            infoText={INFO["RANKED"]["text"]}
          >
            <RankedOrder
              width={800}
              height={500}
              data={probabilities}
              highlight={selectPhenotype}
            />
          </Layout> */}
        </Grid>
      </ThemeProvider>
    </StyledEngineProvider>
  );
};

const PhenotypeSelect = ({ options, subsetParam, onSelect }) => (
  <Select
    width={200}
    options={options}
    value={subsetParam}
    title={"Color"}
    onSelect={onSelect}
  />
);

export default DataWrapper;
