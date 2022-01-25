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

import { createUseStyles } from "react-jss";

const greyColor = "rgb(211 211 211)";
const darkGrey = "rgb(153 153 153)";
const fillGreen = "#47e547";
const green = "#5fd538";

const filterMapping = {
  response: "Response",
  patient: "Patient",
  treatment: "Treatment",
  subtype: "Subtype",
  clone_id: "Clone ID",
  timepoint: "Timepoint",
};

const Filters = ({ filters, selected, setFilters }) => (
  <Box
    sx={{
      width: "100%",
      //maxWidth: 360,
      backgroundColor: "#f5f5f5",
      maxHeight: 400,
      overflowY: "scroll",
      overflowX: "clip",
      //ml: 4,
      marginLeft: 4,
      paddingRight: 2,
    }}
  >
    <List
      style={{ backgroundColor: "#f5f5f5" }}
      component="div"
      aria-labelledby="subheader"
      subheader={
        <ListSubheader
          id="subheader"
          style={{ backgroundColor: "#f5f5f5", position: "relative" }}
        >
          <Grid
            container
            direction="row"
            justifyContent="flex-start"
            alignItems="flex-start"
            wrap="nowrap"
            spacing={2}
            style={{ marginBottom: 10, backgroundColor: "#f5f5f5" }}
          >
            <Grid item xs={8}>
              <Typography
                varient="h3"
                style={{ margin: 5, marginBottom: -12, paddingBottom: 0 }}
              >
                Filter Data
              </Typography>
            </Grid>
            <Grid />
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
const useStyles = createUseStyles({
  root: {
    marginTop: "-2.5px !important",
    marginBottom: "-2.5px !important",
    "& span": {
      marginTop: "15px",
      marginLeft: "3px",
    },
  },
});

const FilterDropdown = ({
  title,
  values,
  onValueClick,
  selected,
  top = true,
  bottom = true,
}) => {
  const classes = useStyles();
  const [open, setOpen] = useState(false);

  const handleClick = () => {
    setOpen(!open);
  };
  const isSelected = selected && selected[0] === title;
  const selectedValueIndex = isSelected
    ? values
        .map((value, index) => ({ v: value, index: index }))
        .filter((v) => v["v"] === selected[1])[0]["index"]
    : -1;

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
        height="40"
        width="19.5"
        style={{
          zIndex: 20,
        }}
      >
        <rect
          x="8"
          y={top ? "0" : "10"}
          width="3"
          height={bottom || open ? 40 : 20}
          style={{
            fill: isSelected ? fillGreen : greyColor,
            stroke: isSelected ? fillGreen : greyColor,
          }}
        />
        <circle
          cx="10"
          cy="15"
          r="6"
          stroke="black"
          stroke-width="1"
          style={{
            fill: isSelected ? fillGreen : greyColor,
            stroke: isSelected ? fillGreen : greyColor,
          }}
        />
      </svg>
      <ListItemText
        primary={filterMapping[title]}
        sx={{ pl: 2, m: 0 }}
        style={{ margin: 5 }}
      />
      {open ? <ExpandLess /> : <ExpandMore />}
    </ListItemButton>,
    <Collapse in={open} timeout="auto" unmountOnExit>
      <List component="div" disablePadding>
        {values.map((value, i) => {
          const isFirstItem = i == 0;
          const isLastItem = i === values.length - 1;
          const isLitUp = selectedValueIndex !== -1 && i <= selectedValueIndex;
          const isSelectedItem = isSelected && selectedValueIndex === i;
          const isMiddleItem = !isFirstItem && !isLastItem;
          return (
            <ListItemButton
              sx={{
                "&.MuiListItemText-root:hover": {
                  bgcolor: "none",
                },
              }}
              style={{ display: "flex", paddingTop: 0, paddingBottom: 0 }}
              key={`${title}-${value}`}
              onClick={() => {
                onValueClick([title, value]);
              }}
              selected={isSelected && selected[1] === value}
            >
              <svg
                height="35"
                width="37"
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
                    fill: greyColor,
                    stroke: greyColor,
                  }}
                />
                //nothing is selected
                {isMiddleItem &&
                  !isSelected && [
                    <rect
                      x={"31"}
                      y={"0"}
                      width="3"
                      height="100%"
                      style={{
                        fill: greyColor,
                        stroke: greyColor,
                      }}
                    />,
                    <circle
                      cx="32"
                      cy="25"
                      r="4"
                      stroke="black"
                      stroke-width="1"
                      style={{
                        stroke: greyColor,
                        fill: darkGrey,
                      }}
                    />,
                  ]}
                //you are below the selected item
                {isMiddleItem &&
                  !isLitUp &&
                  isSelected &&
                  !isSelectedItem && [
                    <rect
                      x={"31"}
                      y={"0"}
                      width="2.5"
                      height="100%"
                      style={{
                        fill: greyColor,
                        stroke: greyColor,
                      }}
                    />,
                    <circle
                      cx="32"
                      cy="25"
                      r="4"
                      stroke="black"
                      stroke-width="1"
                      style={{
                        stroke: greyColor,
                        fill: greyColor,
                      }}
                    />,
                  ]}
                {isMiddleItem &&
                  isSelectedItem && [
                    <rect
                      x={"31"}
                      y={"0"}
                      width="2.5"
                      height="25"
                      style={{
                        fill: fillGreen,
                        stroke: fillGreen,
                      }}
                    />,
                    <rect
                      x={"31"}
                      y={"25"}
                      width="2.5"
                      height="100%"
                      style={{
                        fill: greyColor,
                        stroke: greyColor,
                      }}
                    />,
                    <circle
                      cx="32"
                      cy="25"
                      r="4"
                      stroke="green"
                      stroke-width="2"
                      style={{
                        stroke: fillGreen,
                        fill: fillGreen,
                      }}
                    />,
                  ]}
                //you are a green rect above a selected item
                {isMiddleItem && !isSelectedItem && isSelected && isLitUp && (
                  <rect
                    x={"31"}
                    y={"0"}
                    width="2.5"
                    height="100%"
                    style={{
                      fill: fillGreen,
                      stroke: fillGreen,
                    }}
                  />
                )}
                {isFirstItem && [
                  <path
                    d="M10,0 C9,24 25,8 32.5,25"
                    stroke={isSelectedItem || isLitUp ? fillGreen : greyColor}
                    strokeWidth="3"
                    fill="transparent"
                  />,
                  <rect
                    x={"31"}
                    y={"25"}
                    width="2.5"
                    height="20"
                    style={{
                      fill: isLitUp && !isSelectedItem ? fillGreen : greyColor,
                      stroke:
                        isLitUp && !isSelectedItem ? fillGreen : greyColor,
                    }}
                  />,
                  isLitUp && !isSelectedItem ? (
                    <div />
                  ) : (
                    <circle
                      cx="32"
                      cy="25"
                      r="4"
                      stroke="black"
                      stroke-width="1"
                      style={{
                        stroke: isSelectedItem ? fillGreen : greyColor,
                        fill: isSelectedItem ? fillGreen : darkGrey,
                      }}
                    />
                  ),
                ]}
                {isLastItem && [
                  <rect
                    x={"31"}
                    y={"0"}
                    width="2.5"
                    height="21"
                    style={{
                      fill: isLitUp ? fillGreen : greyColor,
                      stroke: isLitUp ? fillGreen : greyColor,
                    }}
                  />,
                  <path
                    d="M8,40 C12,21 29,39 33,22"
                    stroke={greyColor}
                    strokeWidth="3"
                    fill="transparent"
                  />,
                  <circle
                    cx="32"
                    cy="25"
                    r="4"
                    stroke="black"
                    stroke-width="1"
                    style={{
                      stroke: isSelectedItem ? fillGreen : greyColor,
                      fill: isSelectedItem ? fillGreen : darkGrey,
                    }}
                  />,
                ]}
              </svg>
              <ListItemText
                primary={value}
                className={classes.root}
                sx={{ pl: 4 }}
                style={{
                  fontSize: "12px !important",
                  paddingTop: 0,
                  paddingLeft: 0,
                }}
              />
            </ListItemButton>
          );
        })}
      </List>
    </Collapse>,
  ];
};

export default Filters;
