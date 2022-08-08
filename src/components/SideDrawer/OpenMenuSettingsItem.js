import SettingsIcon from "@mui/icons-material/Settings";
import { AccordionSummary, Accordion, AccordionDetails } from "./Filters";
import ListItem from "@mui/material/ListItem";
import ListItemButton from "@mui/material/ListItemButton";
import ListItemIcon from "@mui/material/ListItemIcon";
import ListItemText from "@mui/material/ListItemText";
import FormControlLabel from "@mui/material/FormControlLabel";
import Switch from "@mui/material/Switch";

import { useData } from "../../provider/dataContext";
const type = "settings";

const OpenMenuSettingsItem = ({ open, setExpand, expanded }) => {
  const [{ settings }, dispatch] = useData();
  console.log(settings);
  return (
    <Accordion expanded={open ? expanded[type] : false}>
      <AccordionSummary aria-controls="settings-panel" id={type + "-panel"}>
        <ListItemButton
          sx={{
            pl: 0,
            minHeight: 48,
            justifyContent: open ? "initial" : "center",
          }}
          onClick={() => setExpand(null, type)}
        >
          <ListItemIcon
            sx={{
              pl: 0,
              minWidth: 0,
              mr: open ? 3 : "auto",
              justifyContent: "center",
            }}
          >
            <SettingsIcon
              key={type + "-icon"}
              onClick={() => setExpand(null, type)}
            />
          </ListItemIcon>
          <ListItemText
            key={type + "-title"}
            primary={"Settings"}
            sx={{ opacity: open ? 1 : 0 }}
          />
        </ListItemButton>
      </AccordionSummary>
      <AccordionDetails>
        <FormControlLabel
          control={
            <Switch
              label="Filter NA"
              checked={settings.filterNA}
              onChange={() => {
                dispatch({
                  type: "setFilterNA",
                  value: { filterNA: !settings.filterNA },
                });
              }}
            />
          }
          label="Filter NA"
        />
      </AccordionDetails>
    </Accordion>
  );
};
export default OpenMenuSettingsItem;
