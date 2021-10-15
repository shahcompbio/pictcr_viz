import React, { useState } from "react";
import Card from "@mui/material/Card";
import CardContent from "@mui/material/CardContent";
import Grid from "@mui/material/Grid";
import Button from "@mui/material/Button";
import Typography from "@mui/material/Typography";
import Box from "@mui/material/Box";
import List from "@mui/material/List";
import ListSubheader from "@mui/material/ListSubheader";
import ListItem from "@mui/material/ListItem";
import ListItemButton from "@mui/material/ListItemButton";
import ListItemText from "@mui/material/ListItemText";
import Collapse from "@mui/material/Collapse";
import Accordion from "@mui/material/Accordion";
import AccordionDetails from "@mui/material/AccordionDetails";
import AccordionSummary from "@mui/material/AccordionSummary";

import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import ExpandLess from "@mui/icons-material/ExpandLess";
import ExpandMore from "@mui/icons-material/ExpandMore";

import makeStyles from "@mui/styles/makeStyles";

// const useStyles = makeStyles({
//   accordian: {
//     //  background: "#f5f5f5",
//     backgroundColor: "#e9eaed",
//     borderRadius: 5,
//     marginBottom: 5,
//     "&.MuiAccordion-root:before": {
//       top: 0,
//     },
//   },
//   accordianList: {
//     width: "100%",
//     paddingTop: 0,
//   },
//   selectedAccordianItem: {
//     backgroundColor: "#e7e7f1",
//     fontSize: 12,
//     paddingTop: 0,
//     paddingBottom: 0,
//   },
//   accordianItem: {
//     fontSize: 12,
//     paddingTop: 0,
//     paddingBottom: 0,
//   },
//   accordianDetails: {
//     padding: 0,
//     paddingLeft: 10,
//     backgroundColor: "#f5f5f5",
//     //marginLeft: 15,
//   },
//   root: {
//     backgroundColor: "#f5f5f5",
//     //  backgroundColor: "white",
//     padding: 15,
//     width: "100%",
//     marginTop: 0,
//     marginLeft: 15,
//     marginBottom: 0,
//   },
//   rootWithMarginTop: {
//     //  backgroundColor: "#f5f5f5",
//     backgroundColor: "none",
//     //  border: "#ababab",
//     //  borderStyle: "solid",
//     //  borderRadius: 5,
//     //  borderWidth: 1,

//     padding: 15,
//     minWidth: 275,
//     marginTop: 15,
//     marginLeft: 15,
//     marginBottom: 0,
//   },
//   clearButton: {
//     float: "right",
//     boxShadow: "none !important",
//     backgroundColor: "#f7f8fb",
//     border: "solid 1px",
//   },
//   clearSelectionWrapper: {
//     paddingTop: 10,
//     paddingBottom: 10,
//     paddingLeft: 0,
//   },
//   heading: {
//     marginLeft: 3,
//   },
//   button: {
//     fontSize: 18,
//     boxShadow: "none !important",
//     backgroundColor: "#f7f8fb",
//     border: "solid 1px",
//   },
//   content: {
//     width: "100%",
//     padding: 0,
//     paddingBottom: "0 !important",
//   },
//   selectionContent: {
//     //width: 200,
//     padding: 10,
//     paddingBottom: "10 !important",
//     background: "white",
//     border: "#ababab",
//     borderStyle: "solid",
//     borderRadius: 5,
//     borderWidth: 1,
//   },
//   key: {
//     fontWeight: "light",
//     fontSize: 20,
//     color: "#b3b3b3",
//   },
//   value: {
//     fontSize: 20,
//     color: "#81a6f9",
//   },
//   overallCells: {
//     fontSize: 18,
//     fontWeight: "bold",
//     fontFamily: "Noto Sans",

//     color: "#81a6f9",
//     paddingTop: 5,
//   },
//   selectedClone: {
//     fontSize: 28,
//     fontFamily: "Noto Sans",
//     color: "black",
//   },
//   selectedCells: {
//     fontWeight: "bold",
//     fontSize: 18,
//     fontFamily: "Noto Sans",
//     color: "black",
//   },
//   overallCellsValue: {
//     fontSize: 25,
//     fontFamily: "Noto Sans",
//     color: "black",
//   },
//   popOver: { width: 350, maxWidth: 500 },
//   bullet: {
//     display: "inline-block",
//     margin: "0 2px",
//     transform: "scale(0.8)",
//   },
//   tabPanel: {
//     height: 300,
//   },
//   tabTitle: {
//     minWidth: "75px !important",
//     textTransform: "none",
//   },
//   light: {
//     fontFamily: "Noto Sans",
//     fontWeight: "light",
//   },
//   title: {
//     fontSize: 14,
//     fill: "black",
//     fontFamily: "Noto Sans",
//     fontWeight: "bold",
//   },
//   pos: {
//     marginBottom: 12,
//   },
// });
const filterMapping = {
  response: "Response",
  patient: "Patient",
  treatment: "Treatment",
  subtype: "Subtype",
  clone_id: "Clone ID",
  timepoint: "Timepoint",
};
const Filters = ({ filters, setFilters, selected }) => {
  //   const classes = useStyles();
  return (
    <Card elevation={0}>
      <CardContent>
        <div style={{ height: 400, overflowY: "scroll", overflowX: "clip" }}>
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
                disabled={!selected}
                // className={classes.clearButton}
                onClick={() => setFilters(null)}
              >
                Clear
              </Button>
            </Grid>
          </Grid>
          <DrawerContent
            filters={filters}
            // classes={classes}
            setFilters={setFilters}
            selected={selected}
          />
        </div>
      </CardContent>
    </Card>
  );
};

const DrawerContent = ({ filters, classes, setFilters, selected }) => {
  const [expanded, setExpanded] = useState({});

  const handleChange = (panel) => (event, isExpanded) => {
    var newExpanded = expanded;
    newExpanded[panel] = isExpanded;
    setExpanded(newExpanded);
  };
  return (
    <div style={{ position: "relative" }}>
      <div
        style={{
          position: "absolute",
          marginLeft: 11,
          borderLeft: "3px solid rgb(211 211 211)",
          height: "100%",
          zIndex: 10,
        }}
      />
      <div
        style={{ position: "relative", borderRight: "solid 10px transparent" }}
      >
        {filters
          .filter((filter) => filter["name"] !== "clone_id")
          .map((filter, index) => (
            <Accordion
              elevation={0}
              // classes={{
              //   root: classes.accordian,
              // }}
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
                <svg
                  height="20"
                  width="10"
                  style={{
                    zIndex: 20,
                  }}
                >
                  <circle
                    cx="5"
                    cy="12"
                    r="5"
                    stroke="black"
                    stroke-width="1"
                    style={{
                      fill:
                        selected && selected[0] === filter["name"]
                          ? "green"
                          : "rgb(211 211 211)",
                      //fill: "rgb(211 211 211)",
                      stroke: "rgb(211 211 211)",
                    }}
                  />
                </svg>
                <Typography
                  // className={classes.heading}
                  key={"panel-title-" + index}
                >
                  {filterMapping[filter["name"]]}
                </Typography>
              </AccordionSummary>
              <AccordionDetails
                key={"panel-details-" + index}
                //   className={classes.accordianDetails}
              >
                <List
                  // className={classes.accordianList}
                  component="nav"
                  key={"panel-list-values-" + index}
                >
                  {filter["values"].map((value, i) => (
                    <ListItem
                      key={"panel-item-" + value}
                      // className={
                      //   selected && selected[1] === value
                      //     ? classes.selectedAccordianItem
                      //     : classes.accordianItem
                      // }
                      button
                      onClick={(event) => {
                        setFilters([filter["name"], value]);
                      }}
                    >
                      <ListItemText
                        primary={value}
                        style={{ fontSize: 12 }}
                        key={"panel-item-text-" + value}
                      />
                    </ListItem>
                  ))}
                </List>
              </AccordionDetails>
            </Accordion>
          ))}
      </div>
    </div>
  );
};

const Filters2 = ({ filters, selected, setFilters }) => (
  <Box
    sx={{
      width: "100%",
      maxWidth: 360,
      bgcolor: "background.paper",
      maxHeight: 400,
      overflowY: "scroll",
      overflowX: "clip",
      ml: 4,
    }}
  >
    <List
      sx={{ bgcolor: "background.paper" }}
      component="div"
      aria-labelledby="nested-list-subheader"
      subheader={
        <ListSubheader id="nested-list-subheader">
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
                disabled={!selected}
                onClick={() => setFilters(null)}
              >
                Clear
              </Button>
            </Grid>
          </Grid>
        </ListSubheader>
      }
    >
      {filters.map((filter, index) => (
        <FilterDropdown
          key={filter["name"]}
          title={filter["name"]}
          values={filter["values"]}
          onValueClick={setFilters}
          selected={selected}
          top={index !== 0}
          bottom={index !== filters.length - 1}
        />
      ))}
    </List>
  </Box>
);

const FilterDropdown = ({
  title,
  values,
  onValueClick,
  selected,
  top = true,
  bottom = true,
}) => {
  const [open, setOpen] = useState(false);

  const handleClick = () => {
    setOpen(!open);
  };

  return [
    <ListItemButton
      onClick={handleClick}
      style={{
        display: "flex",
        paddingTop: 0,
        paddingBottom: 0,
      }}
    >
      <svg
        height="30"
        width="20"
        style={{
          zIndex: 20,
        }}
      >
        <rect
          x="8"
          y={top ? "0" : "10"}
          width="3"
          height={bottom || open ? 30 : 20}
          style={{
            fill: "rgb(211 211 211)",
            stroke: "rgb(211 211 211)",
          }}
        />
        <circle
          cx="10"
          cy="15"
          r="5"
          stroke="black"
          stroke-width="1"
          style={{
            fill:
              selected && selected[0] === title ? "green" : "rgb(211 211 211)",
            stroke: "rgb(211 211 211)",
          }}
        />
      </svg>
      <ListItemText primary={filterMapping[title]} sx={{ pl: 2, m: 0 }} />
      {open ? <ExpandLess /> : <ExpandMore />}
    </ListItemButton>,
    <Collapse in={open} timeout="auto" unmountOnExit>
      <List component="div" disablePadding>
        {values.map((value) => (
          <ListItemButton
            style={{ display: "flex", paddingTop: 0, paddingBottom: 0 }}
            key={`${title}-${value}`}
            onClick={() => {
              onValueClick([title, value]);
            }}
            selected={
              selected && selected[0] === title && selected[1] === value
            }
          >
            <svg
              height="35"
              width="20"
              style={{
                zIndex: 20,
              }}
            >
              <rect
                x="8"
                y="0"
                width="3"
                height="35"
                style={{
                  fill: "rgb(211 211 211)",
                  stroke: "rgb(211 211 211)",
                }}
              />
            </svg>
            <ListItemText primary={value} sx={{ pl: 4 }} />
          </ListItemButton>
        ))}
      </List>
    </Collapse>,
  ];
};

export default Filters2;
