import React, { useState } from "react";
import * as d3 from "d3";
import _ from "lodash";

import { useCanvas, PieChart, drawCanvasAxis } from "@shahlab/planetarium";

import InboxIcon from "@material-ui/icons/Inbox";
import Card from "@material-ui/core/Card";
import CardActions from "@material-ui/core/CardActions";
import CardContent from "@material-ui/core/CardContent";
import Grid from "@material-ui/core/Grid";
import Button from "@material-ui/core/Button";
import Typography from "@material-ui/core/Typography";
import Paper from "@material-ui/core/Paper";
import Popover from "@material-ui/core/Popover";
import Tabs from "@material-ui/core/Tabs";
import Tab from "@material-ui/core/Tab";
import Box from "@material-ui/core/Box";

import List from "@material-ui/core/List";
import ListItem from "@material-ui/core/ListItem";
import ListItemIcon from "@material-ui/core/ListItemIcon";
import ListItemText from "@material-ui/core/ListItemText";

import { Layout } from "@shahlab/planetarium";

import { makeStyles } from "@material-ui/core/styles";

const useStyles = makeStyles({
  root: {
    backgroundColor: "#f5f5f5",
    minWidth: 275,
    marginTop: 15,
    marginLeft: 15,
    marginBottom: 0,
    paddingBottom: 0,
  },
  clearButton: {
    boxShadow: "none !important",
    backgroundColor: "#f7f8fb",
    border: "solid 1px",
  },
  clearSelectionWrapper: {
    paddingTop: 10,
    paddingBottom: 10,
    paddingLeft: 0,
  },
  button: {
    boxShadow: "none !important",
    backgroundColor: "#f7f8fb",
    border: "solid 1px",
  },
  content: {
    width: 200,
    padding: 0,
    paddingBottom: "0 !important",
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
    fontSize: 18,
    fontWeight: "bold",
    fontFamily: "MyFontLight",
    color: "#81a6f9",
    marginLeft: 15,
    paddingTop: 5,
  },
  selectedClone: {
    fontSize: 28,
    fontFamily: "MyFontRegular",
    color: "black",
  },
  selectedCells: {
    fontSize: 25,
    fontFamily: "MyFontRegular",
    color: "black",
  },
  overallCellsValue: {
    fontSize: 25,
    fontFamily: "MyFontRegular",
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
const filterMapping = {
  response: "Response",
  patient: "Patient",
  treatment: "Treatment",
  subtype: "Subtype",
  clone_id: "Clone ID",
};
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
  filters,
}) => {
  console.log(selected);
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
          <Card className={classes.root} elevation={0}>
            <CardContent className={classes.content}>
              <Typography className={classes.key}>Sample:</Typography>
              <Typography className={classes.value}>{sample}</Typography>
            </CardContent>
          </Card>
          <Card className={classes.root} elevation={0}>
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
                  <Typography className={classes.key}>
                    {selectedType}:
                  </Typography>
                  <Grid
                    container
                    direction="row"
                    justify="flex-start"
                    alignItems="stretch"
                  >
                    <Typography variant="h4">{selected} </Typography>
                  </Grid>
                </span>
              )}
              {!highlighted && !selected && (
                <span>
                  <Typography className={classes.key}>Data Points:</Typography>
                  <Grid
                    container
                    direction="row"
                    justify="flex-start"
                    alignItems="stretch"
                  >
                    <Typography className={classes.overallCellsValue}>
                      {data.length}
                    </Typography>
                  </Grid>
                </span>
              )}
            </CardContent>

            {(selected || highlighted) && (
              <CardActions className={classes.clearSelectionWrapper}>
                <Button
                  variant="contained"
                  className={classes.clearButton}
                  onClick={() => setHighlight()}
                >
                  Clear Selection
                </Button>
              </CardActions>
            )}
          </Card>
          <Card className={classes.root} elevation={0}>
            <CardContent className={classes.content}>
              <Filters classes={classes} filters={filters} />
            </CardContent>
          </Card>
        </Grid>
      </Grid>
    </Paper>
  );
};

const Filters = ({ filters, classes }) => {
  const [anchorEl, setAnchorEl] = useState(null);
  const [tabIndex, setTabIndex] = useState(0);

  function a11yProps(index) {
    return {
      id: `scrollable-auto-tab-${index}`,
      "aria-controls": `scrollable-auto-tabpanel-${index}`,
    };
  }
  const [selectedIndex, setSelectedIndex] = useState(1);

  const handleListItemClick = (event, index) => {
    setSelectedIndex(index);
  };

  const open = Boolean(anchorEl);
  const id = open ? "simple-popover" : undefined;
  return (
    <div>
      <Button
        aria-describedby={id}
        variant="contained"
        onClick={(event) => setAnchorEl(event.currentTarget)}
        className={classes.button}
      >
        Filter Data
      </Button>
      <Popover
        className={classes.popOver}
        id={id}
        open={open}
        anchorEl={anchorEl}
        onClose={() => setAnchorEl(null)}
        anchorOrigin={{
          vertical: "bottom",
          horizontal: "left",
        }}
        transformOrigin={{
          vertical: "top",
          horizontal: "left",
        }}
      >
        <Typography className={classes.typography}>
          <Tabs
            value={tabIndex}
            onChange={(event, newValue) => setTabIndex(newValue)}
            indicatorColor="primary"
            textColor="primary"
            variant="scrollable"
            scrollButtons="auto"
            aria-label="scrollable auto tabs"
          >
            {filters.map((filter, index) => (
              <Tab
                label={filterMapping[filter["name"]]}
                {...a11yProps(index)}
                className={classes.tabTitle}
              />
            ))}
          </Tabs>
          {filters.map((filter, index) => {
            return (
              <TabPanel
                className={classes.tabPanel}
                value={tabIndex}
                index={index}
                key={filter["name"] + "TabPanel"}
              >
                <List component="nav" aria-label="main mailbox folders">
                  {filter["values"].map((value, index) => {
                    const labelId = `checkbox-list-label-${value}`;
                    return (
                      <ListItem
                        button
                        selected={selectedIndex === index}
                        onClick={(event) => handleListItemClick(event, 0)}
                      >
                        <ListItemText primary={value} />
                      </ListItem>
                    );
                  })}
                </List>
              </TabPanel>
            );
          })}
        </Typography>
      </Popover>
    </div>
  );
};
const TabPanel = ({ children, value, index, ...other }) => {
  return (
    <div
      role="tabpanel"
      hidden={value !== index}
      id={`scrollable-auto-tabpanel-${index}`}
      aria-labelledby={`scrollable-auto-tab-${index}`}
      {...other}
    >
      {value === index && (
        <Box p={3}>
          <Typography>{children}</Typography>
        </Box>
      )}
    </div>
  );
};
export default MetaData;
