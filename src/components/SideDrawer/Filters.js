import React, { useState, memo } from "react";
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
import Checkbox from "@mui/material/Checkbox";
import Paper from "@mui/material/Paper";
import MenuList from "@mui/material/MenuList";
import MenuItem from "@mui/material/MenuItem";
import ListItemIcon from "@mui/material/ListItemIcon";
import LayersClearIcon from "@mui/icons-material/LayersClear";
import Chip from "@mui/material/Chip";

import { styled } from "@mui/material/styles";
import { useData } from "../../provider/dataContext";

const filterMapping = {
  response: "Response",
  patient: "Patient",
  treatment: "Treatment",
  subtype: "Subtype",
  clone_id: "Clone ID",
  IR_VDJ_1_junction_aa: "Clone ID",
  cell_type: "Cell Type",
  timepoint: "Timepoint",
};

const Filters = () => {
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
        maxHeight: 500,
        overflowY: "scroll",
        overflowX: "clip",
      }}
    >
      <List
        style={{ backgroundColor: "whitesmoke" }}
        component="div"
        aria-labelledby="subheader"
        subheader={
          selectFilters !== null && (
            <ListSubheader
              id="subheader"
              sx={{ position: "relative", paddingLeft: 0 }}
            >
              <Chips />
            </ListSubheader>
          )
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
              dispatch({ type: "setSelectFilters", value: value });
            }}
            selected={selectFilters}
            top={index !== 0}
            bottom={index !== filters.length - 1}
          />
        ))}
      </List>
      <ClearBox
        disabled={!(selectFilters !== null)}
        onClick={() => {
          dispatch({ type: "setSelectFilters", value: null });
        }}
      />
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

export const Accordion = styled((props) => (
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

export const AccordionSummary = styled((props) => (
  <MuiAccordionSummary
    expandIcon={<ArrowForwardIosSharpIcon sx={{ fontSize: "0.9rem" }} />}
    {...props}
  />
))(({ theme }) => ({
  backgroundColor: "white",
  borderRadius: "5px",
  marginTop: "0px",
  marginBottom: "0px",
  //marginBottom: theme.spacing(1),
  "& .MuiAccordionSummary-expandIconWrapper.Mui-expanded": {
    transform: "rotate(90deg)",
    marginRight: theme.spacing(1),
  },
  "& .MuiAccordionSummary-content": {
    marginLeft: theme.spacing(1),
    marginTop: "0px",
    marginBottom: "0px",
  },
}));
const Chips = () => {
  const [{ selectFilters }, dispatch] = useData();

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
      onDelete={() => {
        dispatch({ type: "setSelectFilters", value: null });
      }}
    />
  ) : null;
};
export const AccordionDetails = styled(MuiAccordionDetails, {
  shouldForwardProp: (prop) => true,
})(({ theme, background }) => ({
  padding: "0px 0px 8px 8px",
  borderTop: "1px solid rgba(0, 0, 0, .125)",
  backgroundColor: background ? background : "#f5f5f5",
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
        <AccordionDetails background={"white"}>
          {values.map((value, i) => (
            <ListItemButton
              sx={{
                pl: 0,
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
              <ListItemIcon sx={{ minWidth: "10px" }}>
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
