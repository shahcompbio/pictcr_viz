import React, { useState } from "react";

import Card from "@mui/material/Card";
import CardContent from "@mui/material/CardContent";
import Grid from "@mui/material/Grid";
import Button from "@mui/material/Button";
import Typography from "@mui/material/Typography";
import Paper from "@mui/material/Paper";
import List from "@mui/material/List";
import ListItem from "@mui/material/ListItem";
import ListItemText from "@mui/material/ListItemText";
import Divider from "@mui/material/Divider";

import Accordion from "@mui/material/Accordion";
import AccordionDetails from "@mui/material/AccordionDetails";
import AccordionSummary from "@mui/material/AccordionSummary";

import ExpandMoreIcon from "@mui/icons-material/ExpandMore";

import makeStyles from "@mui/styles/makeStyles";

const useStyles = makeStyles({
  accordian: {
    //  background: "#f5f5f5",
    backgroundColor: "#e9eaed",
    borderRadius: 5,
    marginBottom: 5,
    "&.MuiAccordion-root:before": {
      top: 0,
    },
  },
  accordianList: {
    width: "100%",
    paddingTop: 0,
  },
  selectedAccordianItem: {
    backgroundColor: "#e7e7f1",
    fontSize: 12,
    paddingTop: 0,
    paddingBottom: 0,
  },
  accordianItem: {
    fontSize: 12,
    paddingTop: 0,
    paddingBottom: 0,
  },
  accordianDetails: {
    padding: 0,
    paddingLeft: 10,
    backgroundColor: "#f5f5f5",
    //marginLeft: 15,
  },
  root: {
    backgroundColor: "#f5f5f5",
    //  backgroundColor: "white",
    padding: 15,
    width: "100%",
    marginTop: 0,
    marginLeft: 15,
    marginBottom: 0,
  },
  rootWithMarginTop: {
    //  backgroundColor: "#f5f5f5",
    backgroundColor: "none",
    //  border: "#ababab",
    //  borderStyle: "solid",
    //  borderRadius: 5,
    //  borderWidth: 1,

    padding: 15,
    minWidth: 275,
    marginTop: 15,
    marginLeft: 15,
    marginBottom: 0,
  },
  clearButton: {
    float: "right",
    boxShadow: "none !important",
    backgroundColor: "#f7f8fb",
    marginRight: 7,
    //  border: "solid 1px",
  },
  clearSelectionWrapper: {
    paddingTop: 10,
    paddingBottom: 10,
    paddingLeft: 0,
  },
  heading: {
    marginLeft: 3,
  },
  button: {
    fontSize: 18,
    boxShadow: "none !important",
    backgroundColor: "#f7f8fb",
    border: "solid 1px",
  },
  content: {
    width: "100%",
    padding: 0,
    paddingBottom: "0 !important",
  },
  selectionContent: {
    //width: 200,
    padding: 10,
    paddingBottom: "10 !important",
    background: "#f5f5f5",
    border: "#ababab",
    borderStyle: "solid",
    borderRadius: 5,
    borderWidth: 1,
  },
  key: {
    fontWeight: "light",
    fontSize: 20,
    color: "#b3b3b3",
  },
  value: {
    fontSize: 20,
    color: "#81a6f9",
  },
  overallCells: {
    fontSize: 18,
    fontWeight: "bold",
    fontFamily: "Noto Sans",

    color: "#81a6f9",
    paddingTop: 5,
  },
  selectedClone: {
    fontSize: 28,
    fontFamily: "Noto Sans",
    color: "black",
  },
  selectedCells: {
    fontWeight: "bold",
    fontSize: 18,
    fontFamily: "Noto Sans",
    color: "black",
  },
  overallCellsValue: {
    fontSize: 25,
    fontFamily: "Noto Sans",
    color: "black",
  },
  popOver: { width: 350, maxWidth: 500 },
  bullet: {
    display: "inline-block",
    margin: "0 2px",
    transform: "scale(0.8)",
  },
  tabPanel: {
    height: 300,
  },
  tabTitle: {
    minWidth: "75px !important",
    textTransform: "none",
  },
  light: {
    fontFamily: "Noto Sans",
    fontWeight: "light",
  },
  title: {
    fontSize: 14,
    fill: "black",
    fontFamily: "Noto Sans",
    fontWeight: "bold",
  },
  pos: {
    marginBottom: 12,
  },
});
const filterMapping = {
  response: "Response",
  patient: "Patient",
  treatment: "Treatment",
  subtype: "Subtype",
  clone_id: "Clone ID",
};

const MetaData = ({
  width,
  height,
  data,
  sample,
  selected,
  highlighted,
  selectedType,
  setHighlight,
  filters,
  setFilters,
  selectFilters,
  totalCount,
}) => {
  const classes = useStyles();

  return (
    <Grid
      container
      direction="column"
      justifyContent="flex-start"
      alignItems="stretch"
    >
      <Grid item>
        <Header classes={classes} sample={sample} totalCount={totalCount} />
        <Divider variant="middle" style={{ marginLeft: 30 }} />

        <Card className={classes.root} elevation={0}>
          <CardContent className={classes.content}>
            <Grid
              container
              direction="row"
              justifyContent="space-around"
              alignItems="flex-end"
              style={{ marginBottom: 0 }}
            >
              <Grid item xs={9}>
                {highlighted ? (
                  <Typography
                    varient="h4"
                    style={{
                      color: "black",
                    }}
                  >
                    {highlighted.length} selected
                  </Typography>
                ) : selected ? (
                  <Typography
                    varient="h4"
                    style={{
                      color: "black",
                    }}
                  >
                    {selected.length} selected
                  </Typography>
                ) : (
                  <Typography
                    varient="h4"
                    style={{
                      color: "grey",
                    }}
                  >
                    Selection
                  </Typography>
                )}
              </Grid>
              <Grid item>
                <Button
                  variant="contained"
                  disabled={!selected && !highlighted}
                  className={classes.clearButton}
                  onClick={() => setHighlight()}
                  disableElevation
                  style={{ backgroundColor: "#ffffff" }}
                >
                  Clear
                </Button>
              </Grid>
            </Grid>
          </CardContent>
        </Card>
        <Divider variant="middle" style={{ marginLeft: 30 }} />
      </Grid>
    </Grid>
  );
};

const Header = ({ classes, sample, totalCount }) => (
  <div className={classes.rootWithMarginTop}>
    <Grid
      container
      direction="row"
      justifyContent="flex-start"
      alignItems="center"
    >
      <Typography>Sample:</Typography>
      <Typography>{sample}</Typography>
    </Grid>
    <Grid
      container
      direction="row"
      justifyContent="flex-start"
      alignItems="stretch"
    >
      <Typography>Data Points:</Typography>
      <Typography>{totalCount}</Typography>
    </Grid>
  </div>
);

export default MetaData;
