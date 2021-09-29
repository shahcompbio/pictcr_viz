import React, { useState } from "react";

import Card from "@material-ui/core/Card";
import CardContent from "@material-ui/core/CardContent";
import Grid from "@material-ui/core/Grid";
import Button from "@material-ui/core/Button";
import Typography from "@material-ui/core/Typography";
import Paper from "@material-ui/core/Paper";
import List from "@material-ui/core/List";
import ListItem from "@material-ui/core/ListItem";
import ListItemText from "@material-ui/core/ListItemText";
import Divider from "@material-ui/core/Divider";

import Accordion from "@material-ui/core/Accordion";
import AccordionDetails from "@material-ui/core/AccordionDetails";
import AccordionSummary from "@material-ui/core/AccordionSummary";

import ExpandMoreIcon from "@material-ui/icons/ExpandMore";

import { makeStyles } from "@material-ui/core/styles";

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
    border: "solid 1px",
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
    background: "white",
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
      justify="flex-start"
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
    <Grid container direction="row" justify="flex-start" alignItems="stretch">
      <Typography>Data Points:</Typography>
      <Typography>{totalCount}</Typography>
    </Grid>
  </div>
);

export default MetaData;
