import React, { useState } from "react";

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
import ListItemText from "@material-ui/core/ListItemText";
import Drawer from "@material-ui/core/Drawer";

import Accordion from "@material-ui/core/Accordion";
import AccordionDetails from "@material-ui/core/AccordionDetails";
import AccordionSummary from "@material-ui/core/AccordionSummary";

import ExpandMoreIcon from "@material-ui/icons/ExpandMore";

import { makeStyles } from "@material-ui/core/styles";

const useStyles = makeStyles({
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
    backgroundColor: "#f5f5f5",
    backgroundColor: "white",
    border: "#ababab",
    borderStyle: "solid",
    borderRadius: 5,
    borderWidth: 1,

    padding: 15,
    minWidth: 275,
    marginTop: 15,
    marginLeft: 30,
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
  button: {
    fontSize: 18,
    boxShadow: "none !important",
    backgroundColor: "#f7f8fb",
    border: "solid 1px",
  },
  content: {
    width: 200,
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

  console.log(selectFilters);
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
          <Header classes={classes} sample={sample} totalCount={totalCount} />
          <Card className={classes.root} elevation={0}>
            <CardContent className={classes.content}>
              <Filters
                classes={classes}
                filters={filters}
                setFilters={setFilters}
                setHighlight={setHighlight}
              />
            </CardContent>
          </Card>
          {(highlighted || selected || selectFilters) && (
            <Card className={classes.root} elevation={0}>
              <CardContent className={classes.selectionContent}>
                {highlighted && (
                  <HighlightedContent
                    dataCount={totalCount}
                    highlightedCount={highlighted.length}
                    classes={classes}
                  />
                )}
                {selected && (
                  <SelectedContent
                    classes={classes}
                    selectedType={selectedType}
                    selected={selected}
                  />
                )}
                {selectFilters && (
                  <SelectedFilterContent
                    dataCount={totalCount}
                    classes={classes}
                    highlightedCount={data.length}
                    filter={selectFilters}
                  />
                )}
              </CardContent>

              {(selected || highlighted || selectFilters) && (
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
          )}
        </Grid>
      </Grid>
    </Paper>
  );
};

const Filters = ({ filters, classes, setFilters, setHighlight }) => {
  const [isDrawerOpen, setDrawerOpen] = useState(false);

  return (
    <div style={{ height: 400, overflowY: "scroll" }}>
      <Grid
        container
        direction="row"
        justifyContent="space-around"
        alignItems="flex-end"
        style={{ marginBottom: 10 }}
      >
        <Grid item xs={9}>
          <Typography varient="h4">Filters</Typography>
        </Grid>
        <Grid item>
          <Button
            variant="contained"
            className={classes.clearButton}
            onClick={() => setHighlight()}
          >
            Clear
          </Button>
        </Grid>
      </Grid>
      <DrawerContent
        filters={filters}
        classes={classes}
        setFilters={setFilters}
        setDrawerOpen={setDrawerOpen}
      />
    </div>
  );
};
const DrawerContent = ({ filters, classes, setFilters, setDrawerOpen }) => {
  const [expanded, setExpanded] = useState({});

  const handleChange = (panel) => (event, isExpanded) => {
    var newExpanded = expanded;
    newExpanded[panel] = isExpanded;
    setExpanded(newExpanded);
  };

  return [
    ...filters.map((filter, index) => (
      <Accordion
        key={"panel-" + index}
        expanded={expanded["panel" + index]}
        onChange={handleChange("panel" + index)}
      >
        <AccordionSummary
          expandIcon={<ExpandMoreIcon />}
          aria-controls={"panel" + index + "bh-content"}
          id={"panel" + index + "bh-header"}
          key={"panel-summary-" + index}
        >
          <Typography className={classes.heading} key={"panel-title-" + index}>
            {filterMapping[filter["name"]]}
          </Typography>
        </AccordionSummary>
        <AccordionDetails key={"panel-details-" + index}>
          <List component="nav" key={"panel-list-values-" + index}>
            {filter["values"].map((value, i) => (
              <ListItem
                key={"panel-item-" + value}
                style={{ fontSize: 14 }}
                button
                //  selected={index === selectedTabIndex && selectedIndex === i}
                onClick={(event) => {
                  //  setSelectedTabIndex(tabIndex);
                  //  setSelectedIndex(i);
                  setFilters([filter["name"], value]);
                  handleChange("panel" + index);
                  setDrawerOpen(false);
                  // handleListItemClick(event, index, value);
                }}
              >
                <ListItemText
                  primary={value}
                  style={{ fontSize: 14 }}
                  key={"panel-item-text-" + value}
                />
              </ListItem>
            ))}
          </List>
        </AccordionDetails>
      </Accordion>
    )),
  ];
};
/*      <Popover
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
    >*/
/*const FilterDrawerContent = ({ setTabIndex, tabIndex, filters }) => (
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
    {filters.map((filter, tindex) => {
      return (
        <TabPanel
          className={classes.tabPanel}
          value={tabIndex}
          index={tindex}
          key={filter["name"] + "TabPanel"}
        >
          <List component="nav" aria-label="main mailbox folders">
            {filter["values"].map((value, index) => (
              <ListItem
                button
                selected={
                  tabIndex === selectedTabIndex && selectedIndex === index
                }
                onClick={(event) => {
                  setSelectedTabIndex(tabIndex);
                  setSelectedIndex(index);
                  setFilters([filter["name"], value]);
                  // handleListItemClick(event, index, value);
                }}
              >
                <ListItemText primary={value} />
              </ListItem>
            ))}
          </List>
        </TabPanel>
      );
    })}
  </Typography>
);*/

const Header = ({ classes, sample, totalCount }) => (
  <Card className={classes.rootWithMarginTop} elevation={0}>
    <CardContent className={classes.content}>
      <Grid
        container
        direction="row"
        justifyContent="flex-start"
        alignItems="center"
      >
        <Typography className={classes.key}>Sample:</Typography>
        <Typography className={classes.value}>{sample}</Typography>
      </Grid>
      <Grid container direction="row" justify="flex-start" alignItems="stretch">
        <Typography className={classes.key}>Data Points:</Typography>
        <Typography className={classes.value}>{totalCount}</Typography>
      </Grid>
    </CardContent>
  </Card>
);
const SelectedFilterContent = ({
  classes,
  dataCount,
  highlightedCount,
  filter,
}) => (
  <span>
    <Grid
      container
      direction="column"
      justify="flex-start"
      alignItems="stretch"
    >
      <Grid>
        <Typography variant="h7" className={classes.key}>
          Filtered by:{" "}
          {filter[0][0].toUpperCase() +
            filter[0].substring(1, filter[0].length)}
        </Typography>
      </Grid>
      <Grid>
        <Typography variant="h4">{filter[1]} </Typography>
      </Grid>
      <Grid container direction="row" justify="flex-start" alignItems="stretch">
        <Typography className={classes.selectedCells}>
          {highlightedCount}
        </Typography>
        <Typography className={classes.overallCells}>
          /{dataCount} selected
        </Typography>
      </Grid>
    </Grid>
  </span>
);

const NoSelectionContent = ({ classes, dataCount }) => (
  <span>
    <Typography className={classes.key}>Data Points:</Typography>
    <Grid container direction="row" justify="flex-start" alignItems="stretch">
      <Typography className={classes.overallCellsValue}>{dataCount}</Typography>
    </Grid>
  </span>
);
const SelectedContent = ({ classes, selectedType, selected }) => (
  <span>
    <Typography className={classes.key}>{selectedType}:</Typography>
    <Grid container direction="row" justify="flex-start" alignItems="stretch">
      <Typography variant="h4">{selected} </Typography>
    </Grid>
  </span>
);

const HighlightedContent = ({ classes, highlightedCount, dataCount }) => (
  <span>
    <Typography className={classes.key}>Highlighted:</Typography>
    <Grid container direction="row" justify="flex-start" alignItems="stretch">
      <Typography className={classes.selectedCells}>
        {highlightedCount}
      </Typography>
      <Typography className={classes.overallCells}>
        /{dataCount} selected
      </Typography>
    </Grid>
  </span>
);

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
