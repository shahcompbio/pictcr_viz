import FilterListIcon from "@mui/icons-material/FilterList";
import { AccordionSummary, Accordion, AccordionDetails } from "./Filters";
import ListItem from "@mui/material/ListItem";
import ListItemButton from "@mui/material/ListItemButton";
import ListItemIcon from "@mui/material/ListItemIcon";
import ListItemText from "@mui/material/ListItemText";
import KeyboardArrowDownIcon from "@mui/icons-material/KeyboardArrowDown";
import IconButton from "@mui/material/IconButton";

import Filters from "./Filters";
const type = "filter";
const OpenMenuFilterItem = ({ open, setExpand, expanded }) => (
  <Accordion expanded={open ? expanded[type] : false}>
    <AccordionSummary
      aria-controls={type + "-content"}
      id={type + "-summary"}
      style={{ borderRadius: "0px", border: "0px white" }}
    >
      <ListItemButton
        sx={{
          pl: 0,
          minHeight: 48,
          justifyContent: open ? "initial" : "center",
        }}
        onClick={() => setExpand(null)}
      >
        <ListItemIcon
          sx={{
            pl: 0,
            minWidth: 0,
            mr: open ? 3 : "auto",
            justifyContent: "center",
          }}
        >
          <FilterListIcon onClick={() => setExpand(null, type)} />
        </ListItemIcon>
        <ListItemText primary={"Filter Data"} sx={{ opacity: open ? 1 : 0 }} />
      </ListItemButton>
    </AccordionSummary>
    <AccordionDetails>
      <Filters />
    </AccordionDetails>
  </Accordion>
);
export default OpenMenuFilterItem;
