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
      maxWidth: 360,
      backgroundColor: "#f5f5f5",
      maxHeight: 400,
      overflowY: "scroll",
      overflowX: "clip",
      ml: 4,
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
            justifyContent="space-around"
            alignItems="flex-end"
            style={{ marginBottom: 10, backgroundColor: "#f5f5f5" }}
          >
            <Grid item xs={9}>
              <Typography varient="h4">Filters</Typography>
            </Grid>
            <Grid item>
              <Button
                disableElevation
                variant="contained"
                disabled={!selected}
                style={{ backgroundColor: "#ffffff" }}
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
        height="30"
        width="19.5"
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
            fill: isSelected ? green : greyColor,
            stroke: isSelected ? green : greyColor,
          }}
        />
        <circle
          cx="10"
          cy="15"
          r="6"
          stroke="black"
          stroke-width="1"
          style={{
            fill: isSelected ? green : greyColor,
            stroke: isSelected ? green : greyColor,
          }}
        />
      </svg>
      <ListItemText primary={filterMapping[title]} sx={{ pl: 2, m: 0 }} />
      {open ? <ExpandLess /> : <ExpandMore />}
    </ListItemButton>,
    <Collapse in={open} timeout="auto" unmountOnExit>
      <List component="div" disablePadding>
        {values.map((value, i) => {
          const isFirstItem = i == 0;
          const isLastItem = i === values.length - 1;
          const isLitUp = selectedValueIndex !== -1 && i <= selectedValueIndex;
          const isSelectedItem = isSelected && selectedValueIndex === i;
          return (
            <ListItemButton
              style={{ display: "flex", paddingTop: 0, paddingBottom: 0 }}
              key={`${title}-${value}`}
              onClick={() => {
                onValueClick([title, value]);
              }}
              selected={isSelected && selected[1] === value}
            >
              <svg
                height="35"
                width="32"
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
                <rect
                  x={isFirstItem ? "5" : isLastItem ? "30" : "27"}
                  y={isFirstItem ? "4" : isLastItem ? "-10" : "0"}
                  width="3"
                  rx={isFirstItem ? "10" : isLastItem ? "10" : "0"}
                  height={isFirstItem ? "30" : isLastItem ? "25" : "35"}
                  style={{
                    transform: isFirstItem
                      ? "rotate(316deg)"
                      : isLastItem
                      ? "rotate(49deg)"
                      : 0,
                    fill: isLitUp
                      ? isLastItem
                        ? greyColor
                        : green
                      : greyColor,
                    stroke: isLitUp
                      ? isLastItem
                        ? greyColor
                        : green
                      : greyColor,
                  }}
                />
                {!isFirstItem && !isLastItem && !isLitUp ? (
                  <circle
                    cx="28"
                    cy="15"
                    r="3"
                    stroke="black"
                    stroke-width="1"
                    style={{
                      stroke: darkGrey,
                      fill: isSelectedItem ? fillGreen : darkGrey,
                    }}
                  />
                ) : !isFirstItem && !isLastItem && isSelectedItem ? (
                  <circle
                    cx="28"
                    cy="15"
                    r="3"
                    stroke="green"
                    stroke-width="2"
                    style={{
                      stroke: darkGrey,
                      fill: isSelectedItem ? fillGreen : darkGrey,
                    }}
                  />
                ) : (
                  <div />
                )}
                {isFirstItem && [
                  <rect
                    x={"27"}
                    y={"18"}
                    width="3"
                    height="20"
                    style={{
                      fill: isLitUp ? green : greyColor,
                      stroke: isLitUp ? green : greyColor,
                    }}
                  />,
                  isLitUp && !isSelectedItem ? (
                    <div />
                  ) : (
                    <circle
                      cx="28.5"
                      cy="18"
                      r="3"
                      stroke="black"
                      stroke-width="1"
                      style={{
                        stroke: greyColor,
                        fill: isSelectedItem ? fillGreen : darkGrey,
                      }}
                    />
                  ),
                ]}
                {isLastItem && [
                  <rect
                    x={"27"}
                    y={"0"}
                    width="3"
                    height="17"
                    style={{
                      fill: isLitUp ? green : greyColor,
                      stroke: isLitUp ? green : greyColor,
                    }}
                  />,
                  <circle
                    cx="28.5"
                    cy="15"
                    r="3"
                    stroke="black"
                    stroke-width="1"
                    style={{
                      stroke: greyColor,
                      fill: isSelectedItem ? fillGreen : darkGrey,
                    }}
                  />,
                ]}
              </svg>
              <ListItemText
                primary={value}
                sx={{ pl: 4 }}
                style={{ fontSize: "12px !important" }}
              />
            </ListItemButton>
          );
        })}
      </List>
    </Collapse>,
  ];
};

export default Filters;
