import React, { useState } from "react";
import _ from "lodash";
import * as d3 from "d3";

import ClonotypeUMAP from "./components/Umap";
import ClonotypeExpansion from "./components/ClonotypeExpansion";
import DEGTable from "./components/DEGTable";
import RankedOrder from "./components/RankedOrder";
import Doughnut from "./components/Doughnut";
import MetaData from "./components/MetaData";
import Header from "./components/Header";

import {
  Heatmap,
  ProbabilityHistogram,
  Layout,
  UMAP,
} from "@shahlab/planetarium";

import { makeStyles } from "@material-ui/core/styles";
import Popper from "@material-ui/core/Popper";
import Typography from "@material-ui/core/Typography";
import Grid from "@material-ui/core/Grid";
import Button from "@material-ui/core/Button";
import Card from "@material-ui/core/Card";
import CardContent from "@material-ui/core/CardContent";
import CardHeader from "@material-ui/core/CardHeader";

import { theme } from "./theme/theme";
import { MuiThemeProvider } from "@material-ui/core/styles";
import CssBaseline from "@material-ui/core/CssBaseline";

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
    probabilities={data["probabilities"]}
    degs={data["degs"]}
  />
);

export const VDJ = ({ metadata, degs }) => {
  const [selectPhenotype, setSelectPhenotype] = useState(null);
  const [selectClone, setSelectClone] = useState(null);
  const [selectIDs, setSelectIDs] = useState(null);
  const [activeGraph, setActiveGraph] = useState(null);

  const { clonotypeParam, subtypeParam, logProbParam, xParam, yParam } =
    CONSTANTS;

  // Remove none
  const clonotypeCounts = _.countBy(
    metadata.filter((datum) => datum[clonotypeParam] !== "None"),
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

  const phenotypeValues = Object.keys(_.groupBy(metadata, subtypeParam)).sort();
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
      ? metadata.filter((datum) => selectIDs.includes(datum["cell_id"]))
      : selectClone !== null
      ? metadata.filter((datum) => datum[clonotypeParam] === selectClone)
      : selectPhenotype !== null
      ? metadata.filter((datum) => datum[subtypeParam] === selectPhenotype)
      : null;

  const highlightIDs =
    highlightData === null
      ? null
      : highlightData.map((datum) => datum["cell_id"]);

  const probabilities = metadata.filter(
    (datum) =>
      datum[clonotypeParam] !== "None" || datum[logProbParam] !== "None"
  );

  const subtypeTotals = _.countBy(metadata, subtypeParam);

  return (
    <MuiThemeProvider theme={theme}>
      <CssBaseline />
      <Header />
      <Grid
        container
        direction="column"
        justify="flex-start"
        alignItems="flex-start"
      >
        <Grid
          item
          container
          direction="row"
          justify="flex-start"
          alignItems="flex-start"
        >
          <MetaData
            width={250}
            data={metadata}
            sample="SAMPLE-TITLE-NDVL"
            highlighted={highlightData}
            selected={selectClone || selectPhenotype}
            setHighlight={() => {
              setSelectIDs(null);
              setSelectClone(null);
              setSelectPhenotype(null);
              setActiveGraph(null);
            }}
            selectedType={selectClone ? "Clone" : "Phenotype"}
          />
          <Doughnut
            data={highlightData || metadata}
            type={"CLONOTYPEDOUGH"}
            colors={CLONOTYPE_COLORS}
            width={450}
            height={350}
            otherSubsetParam={subtypeParam}
            subsetParam={clonotypeParam}
          />
          <Doughnut
            data={highlightData || metadata}
            type={"SUBTYPEDOUGH"}
            colors={CLONOTYPE_COLORS}
            width={450}
            height={350}
            otherSubsetParam={clonotypeParam}
            subsetParam={subtypeParam}
          />
        </Grid>
        <Grid
          item
          container
          direction="row"
          justify="flex-start"
          alignItems="flex-start"
        >
          <Layout title={INFO["UMAP"]["title"]} infoText={INFO["UMAP"]["text"]}>
            <ClonotypeUMAP
              width={800}
              height={600}
              data={metadata}
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
          <Layout
            title={INFO["SUBTYPEUMAP"]["title"]}
            infoText={INFO["SUBTYPEUMAP"]["text"]}
          >
            <UMAP
              width={700}
              height={600}
              data={metadata}
              xParam={xParam}
              yParam={yParam}
              subsetParam={subtypeParam}
              idParam="cell_id"
              colorScale={phenotypeColorScale}
              onLasso={(data) => {
                setSelectIDs(
                  data === null ? null : data.map((datum) => datum["cell_id"])
                );
                setActiveGraph(data === null ? null : "phenoUMAP");
              }}
              onLegendClick={(value) => {
                setSelectPhenotype(value);
                setActiveGraph(value === null ? null : "phenoUMAP");
              }}
              disable={activeGraph !== null && activeGraph !== "phenoUMAP"}
              highlightIDs={highlightIDs}
            />
          </Layout>
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
          <RankedOrder width={800} height={500} data={probabilities} />
        </Layout>
      </Grid>
    </MuiThemeProvider>
  );
};
const useStyles = makeStyles({
  root: {
    minWidth: 275,
  },
  header: { padding: 10, paddingBottom: 0 },
  body: {
    padding: 5,
    paddingLeft: 20,
  },
  button: { margin: 5, float: "right" },
  poppr: {
    width: 150,
    float: "right",
    right: "100px",
    top: "10px",
    left: "auto",
    margin: 10,
  },
});
const Popup = ({ selected, setSelected, type }) => {
  const classes = useStyles();
  return (
    <Popper
      open={true}
      placement={"bottom"}
      transition
      className={classes.popper}
    >
      <Card className={classes.root} variant="outlined">
        <CardHeader className={classes.header} title={"Selected " + type} />
        <CardContent className={classes.body}>
          <Typography variant="body">{selected}</Typography>
        </CardContent>
        <Button
          color="primary"
          size="small"
          variant="outlined"
          className={classes.button}
          onClick={setSelected}
        >
          Clear
        </Button>
      </Card>
    </Popper>
  );
};
export default DataWrapper;
