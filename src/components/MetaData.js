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

import { Layout } from "@shahlab/planetarium";

import { makeStyles } from "@material-ui/core/styles";

const useStyles = makeStyles({
  root: {
    minWidth: 275,
  },
  bullet: {
    display: "inline-block",
    margin: "0 2px",
    transform: "scale(0.8)",
  },
  title: {
    fontSize: 14,
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

const MetaData = ({ width, height, data }) => {
  const classes = useStyles();
  return (
    <Layout>
      <Card className={classes.root}>
        <CardContent>
          <Typography
            className={classes.title}
            color="textSecondary"
            gutterBottom
          >
            Sample Name
          </Typography>
        </CardContent>
      </Card>
      <Card className={classes.root}>
        <CardContent>
          <Typography
            className={classes.title}
            color="textSecondary"
            gutterBottom
          >
            Selected Cells
          </Typography>
        </CardContent>
        <CardActions>
          <Button size="small">Clear</Button>
        </CardActions>
      </Card>
      <Card className={classes.root}>here</Card>
    </Layout>
  );
};

export default MetaData;
