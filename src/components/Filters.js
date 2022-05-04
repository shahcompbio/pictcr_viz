import React, { useState, memo } from "react";
import CardContent from "@mui/material/CardContent";
import Grid from "@mui/material/Grid";
import Typography from "@mui/material/Typography";
import Box from "@mui/material/Box";
import List from "@mui/material/List";
import ListSubheader from "@mui/material/ListSubheader";
import ListItemButton from "@mui/material/ListItemButton";
import ListItemText from "@mui/material/ListItemText";
import MuiAccordion from "@mui/material/Accordion";
import MuiAccordionSummary from "@mui/material/AccordionSummary";
import MuiAccordionDetails from "@mui/material/AccordionDetails";
import ArrowForwardIosSharpIcon from "@mui/icons-material/ArrowForwardIosSharp";
import ExpandLess from "@mui/icons-material/ExpandLess";
import Checkbox from "@mui/material/Checkbox";
import Paper from "@mui/material/Paper";

import MenuList from "@mui/material/MenuList";
import MenuItem from "@mui/material/MenuItem";
import ListItemIcon from "@mui/material/ListItemIcon";
import LayersClearIcon from "@mui/icons-material/LayersClear";
import Chip from "@mui/material/Chip";

import { styled } from "@mui/material/styles";
import { createUseStyles } from "react-jss";
import { useData } from "../provider/dataContext";

const filterMapping = {
  response: "Response",
  patient: "Patient",
  treatment: "Treatment",
  subtype: "Subtype",
  clone_id: "Clone ID",
  timepoint: "Timepoint",
};

const Filters = ({ hasSelection, setHighlight }) => {
  const [{ selectFilters, filters }, dispatch] = useData();
  const [expanded, setExpanded] = useState(
    filters.reduce((final, curr) => {
      final[curr.name] =
        selectFilters && selectFilters[0] === curr.name ? true : false;
      return final;
    }, {})
  );
  return (
    <Box
      sx={{
        width: "100%",
        //maxWidth: 360,
        //backgroundColor: "#f5f5f5",
        maxHeight: 500,
        overflowY: "scroll",
        overflowX: "clip",
        //ml: 4,
        marginLeft: 4,
        paddingRight: 2,
      }}
    >
      <List
        style={{ backgroundColor: "white" }}
        component="div"
        aria-labelledby="subheader"
        subheader={
          <ListSubheader
            id="subheader"
            sx={{ position: "relative", paddingLeft: 0 }}
          >
            <Grid
              container
              direction="row"
              justifyContent="flex-start"
              alignItems="flex-start"
              wrap="nowrap"
              spacing={2}
              style={{ marginBottom: 10 }}
            >
              <Grid item xs={8}>
                <Typography
                  variant="h5"
                  style={{
                    marginBottom: -12,
                    paddingBottom: 0,
                  }}
                >
                  Filter Data
                </Typography>
              </Grid>
              <Grid />
            </Grid>
            <Chips />
          </ListSubheader>
        }
      >
        {filters.map((filter, index) => (
          <FilterDropdown
            setExpand={(title) => {
              setExpanded({ ...expanded, [title]: !expanded[title] });
            }}
            expanded={expanded}
            key={filter["name"]}
            title={filter["name"]}
            values={filter["values"]}
            onValueClick={(value) => {
              dispatch({ type: "setSelectedFilters", value: value });
            }}
            selected={selectFilters}
            top={index !== 0}
            bottom={index !== filters.length - 1}
          />
        ))}
      </List>
      <ClearBox disabled={!hasSelection} onCLick={setHighlight} />
    </Box>
  );
};
export const ClearBox = ({ disabled, onClick, clearText = "Clear" }) => (
  <Paper
    sx={{
      position: "sticky",
      width: "100%",
      border: `1px solid grey`,
    }}
  >
    <MenuList style={{ padding: 0 }}>
      <MenuItem
        disabled={disabled}
        onClick={() => {
          onClick();
        }}
        style={{ height: "100%" }}
      >
        <ListItemIcon>
          <LayersClearIcon fontSize="small" />
        </ListItemIcon>
        <ListItemText>{clearText}</ListItemText>
        <Typography variant="body2" color="text.secondary">
          ⌘X
        </Typography>
      </MenuItem>
    </MenuList>
  </Paper>
);

const Accordion = styled((props) => (
  <MuiAccordion disableGutters elevation={0} square {...props} />
))(({ theme }) => ({
  //backgroundColor: "white !important",
  border: `1px solid ${theme.palette.divider}`,
  marginBottom: theme.spacing(1),
  borderRadius: "5px",
  "&:not(:last-child)": {
    borderBottom: 0,
  },
  "&:before": {
    display: "none",
  },
}));

const AccordionSummary = styled((props) => (
  <MuiAccordionSummary
    expandIcon={<ArrowForwardIosSharpIcon sx={{ fontSize: "0.9rem" }} />}
    {...props}
  />
))(({ theme }) => ({
  backgroundColor: "white",
  borderRadius: "5px",
  //marginBottom: theme.spacing(1),
  "& .MuiAccordionSummary-expandIconWrapper.Mui-expanded": {
    transform: "rotate(90deg)",
    marginRight: theme.spacing(1),
  },
  "& .MuiAccordionSummary-content": {
    marginLeft: theme.spacing(1),
  },
}));
const Chips = () => {
  const [{ selectFilters, filters }, dispatch] = useData();
  return selectFilters !== null ? (
    <Chip
      label={
        <span>
          <span style={{ fontWeight: "bold" }}>
            {filterMapping[selectFilters[0]]}:{" "}
          </span>
          <span>{selectFilters[1]}</span>
        </span>
      }
      sx={{ paddingLeft: "-16px" }}
      variant="outlined"
      onClick={() => {}}
      onDelete={() => {}}
    />
  ) : null;
};
const AccordionDetails = styled(MuiAccordionDetails)(({ theme }) => ({
  padding: theme.spacing(2),
  borderTop: "1px solid rgba(0, 0, 0, .125)",
  backgroundColor: "#f5f5f5",
}));

const FilterDropdown = ({
  expanded,
  setExpand,
  title,
  values,
  onValueClick,
  selected,
  top = true,
  bottom = true,
}) => {
  const isSelected = selected && selected[0] === title;

  return (
    <div>
      <Accordion expanded={expanded[title]} onChange={() => setExpand(title)}>
        <AccordionSummary aria-controls="panel1d-content" id="panel1d-header">
          <Typography>{filterMapping[title]}</Typography>
        </AccordionSummary>
        <AccordionDetails>
          {values.map((value, i) => (
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
              //selected={isSelected && selected[1] === value}
            >
              <ListItemIcon>
                <Checkbox
                  edge="start"
                  sx={{
                    color: "#7d867d",

                    "&.Mui-checked": {
                      color: "#3D70B2",
                    },
                  }}
                  checked={isSelected && selected[1] === value}
                  tabIndex={-1}
                  disableRipple
                  inputProps={{ "aria-labelledby": value + "label" }}
                />
              </ListItemIcon>
              <ListItemText
                primary={value}
                sx={{ pl: 4 }}
                style={{
                  fontSize: "12px !important",
                  paddingTop: 0,
                  paddingLeft: 0,
                }}
              />
            </ListItemButton>
          ))}
        </AccordionDetails>
      </Accordion>
    </div>
  );
};

export default memo(Filters);
