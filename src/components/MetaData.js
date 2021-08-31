import React from "react";
import * as d3 from "d3";
import _ from "lodash";

import { useCanvas, PieChart, drawCanvasAxis } from "@shahlab/planetarium";

import Card from "@material-ui/core/Card";
import CardActions from "@material-ui/core/CardActions";
import CardContent from "@material-ui/core/CardContent";
import Grid from "@material-ui/core/Grid";
import Button from "@material-ui/core/Button";
import Typography from "@material-ui/core/Typography";
import Paper from "@material-ui/core/Paper";

import { Layout } from "@shahlab/planetarium";

import { makeStyles } from "@material-ui/core/styles";

const useStyles = makeStyles({
  root: {
    minWidth: 275,
    marginTop: 15,
    marginLeft: 15,
    paddingBottom: 0,
  },
  content: {
    paddingBottom: 0,
  },
  key: {
    fontWeight: "light",
    fontSize: 15,
    color: "#b3b3b3",
  },
  value: {
    fontSize: 20,
    color: "#81a6f9",
  },
  overallCells: {
    fontSize: 25,
    fontFamily: "MyFontLight",
    color: "#81a6f9",
    marginLeft: 15,
    paddingTop: 5,
  },
  selectedCells: {
    fontSize: 30,
    fontFamily: "MyFontRegular",
    color: "black",
  },
  overallCellsValue: {
    fontSize: 30,
    fontFamily: "MyFontRegular",
    color: "black",
  },
  bullet: {
    display: "inline-block",
    margin: "0 2px",
    transform: "scale(0.8)",
  },
  light: {
    fontFamily: "MyFontLight",
  },
  title: {
    fontSize: 14,
    fill: "black",
    fontFamily: "MyFontBold",
  },
  pos: {
    marginBottom: 12,
  },
});
const COLOR_ARRAY = [
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
const PADDING = 50;
const format = d3.format(".3f");

const MetaData = ({
  width,
  height,
  data,
  sample,
  selected,
  highlighted,
  selectedType,
  setHighlight,
}) => {
  const classes = useStyles();
  return (
    <Paper
      elevation={0}
      style={{
        background: "none",
        margin: 10,
        padding: 10,
      }}
    >
      <Grid
        container
        direction="column"
        justify="flex-start"
        alignItems="stretch"
      >
        <Grid item>
          <Card className={classes.root}>
            <CardContent>
              <Typography className={classes.key}>Sample:</Typography>
              <Typography className={classes.value}>{sample}</Typography>
            </CardContent>
          </Card>
          <Card className={classes.root}>
            <CardContent className={classes.content}>
              {highlighted && (
                <span>
                  <Typography className={classes.key}>Data Points:</Typography>
                  <Grid
                    container
                    direction="row"
                    justify="flex-start"
                    alignItems="stretch"
                  >
                    <Typography className={classes.selectedCells}>
                      {highlighted.length}
                    </Typography>
                    <Typography className={classes.overallCells}>
                      / {data.length} selected
                    </Typography>
                  </Grid>
                </span>
              )}
              {selected && (
                <span>
                  <Typography className={classes.key}>Data Points:</Typography>
                  <Grid
                    container
                    direction="row"
                    justify="flex-start"
                    alignItems="stretch"
                  >
                    <Typography variant="h4">{selectedType}</Typography>
                    <Typography className={classes.selectedCells}>
                      {selected.length}
                    </Typography>
                    <Typography className={classes.overallCells}>
                      / {data.length} selected
                    </Typography>
                  </Grid>
                </span>
              )}
              {!highlighted && !selected && (
                <Grid
                  container
                  direction="row"
                  justify="flex-start"
                  alignItems="stretch"
                >
                  <Typography className={classes.overallCellsValue}>
                    {data.length}
                  </Typography>
                  <Typography className={classes.overallCells}>
                    data points
                  </Typography>
                </Grid>
              )}
            </CardContent>

            {(selected || highlighted) && (
              <CardActions>
                <Button
                  size="small"
                  color="secondary"
                  onClick={() => setHighlight()}
                >
                  Clear
                </Button>
              </CardActions>
            )}
          </Card>
        </Grid>
      </Grid>
    </Paper>
  );
};

export default MetaData;
