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
    margin: 15,
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

const MetaData = ({ width, height, data, sample, selected }) => {
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
              <Typography
                className={classes.title}
                color="textSecondary"
                gutterBottom
              >
                Sample:{sample}
              </Typography>
            </CardContent>
          </Card>
          <Card className={classes.root}>
            <CardContent>
              {selected ? (
                <Typography className={classes.light}>
                  {selected.length} / {data.length}
                </Typography>
              ) : (
                <Typography className={classes.title}>
                  {data.length} / 300000 cells
                </Typography>
              )}
            </CardContent>
          </Card>
        </Grid>
      </Grid>
    </Paper>
  );
};

export default MetaData;
