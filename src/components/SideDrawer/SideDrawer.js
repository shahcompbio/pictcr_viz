import { useState, useEffect } from "react";
import { styled } from "@mui/material/styles";
import Box from "@mui/material/Box";
import MuiDrawer from "@mui/material/Drawer";
import List from "@mui/material/List";
import Divider from "@mui/material/Divider";
import IconButton from "@mui/material/IconButton";
import ChevronLeftIcon from "@mui/icons-material/ChevronLeft";
import ChevronRightIcon from "@mui/icons-material/ChevronRight";
import ListItem from "@mui/material/ListItem";
import ListItemButton from "@mui/material/ListItemButton";
import ListItemIcon from "@mui/material/ListItemIcon";
import ListItemText from "@mui/material/ListItemText";
import ListItemWrapper from "./ListItemWrapper";
import ViewCompactIcon from "@mui/icons-material/ViewCompact";
import HelpIcon from "@mui/icons-material/Help";
import SettingsIcon from "@mui/icons-material/Settings";
import FilterListIcon from "@mui/icons-material/FilterList";

import OpenMenuSettingsItem from "./OpenMenuSettingsItem";
import OpenMenuFilterItem from "./OpenMenuFilterItem";

const drawerWidth = 240;

const openedMixin = (theme) => ({
  width: drawerWidth,
  transition: theme.transitions.create("width", {
    easing: theme.transitions.easing.sharp,
    duration: theme.transitions.duration.enteringScreen,
  }),
  overflowX: "hidden",
});

const closedMixin = (theme) => ({
  transition: theme.transitions.create("width", {
    easing: theme.transitions.easing.sharp,
    duration: theme.transitions.duration.leavingScreen,
  }),
  overflowX: "hidden",
  width: `calc(${theme.spacing(15)} + 10px)`,
  [theme.breakpoints.up("sm")]: {
    width: `calc(${theme.spacing(8)} + 10px)`,
  },
});

const DrawerHeader = styled("div")(({ theme }) => ({
  display: "flex",
  alignItems: "center",
  justifyContent: "flex-end",
  padding: theme.spacing(0, 1),
  // necessary for content to be below app bar
  ...theme.mixins.toolbar,
}));

const Drawer = styled(MuiDrawer, {
  shouldForwardProp: (prop) => prop !== "open",
})(({ theme, open }) => ({
  width: drawerWidth,
  flexShrink: 0,
  whiteSpace: "nowrap",
  boxSizing: "border-box",
  ...(open && {
    ...openedMixin(theme),
    "& .MuiDrawer-paper": openedMixin(theme),
  }),
  ...(!open && {
    ...closedMixin(theme),
    "& .MuiDrawer-paper": closedMixin(theme),
  }),
}));

const SideDrawer = ({ view, setView }) => {
  const [open, setOpen] = useState(false);
  const [expanded, setExpand] = useState({ filter: false, view: false });
  useEffect(() => {
    if (open) {
      setExpanded(true, "filter");
    } else {
      setExpanded(false, "filter");
    }
  }, [open]);

  const setExpanded = (condition, type) =>
    setExpand({
      ...expanded,
      [type]: condition === null ? !expanded[type] : condition,
    });

  return (
    <Box sx={{ display: "flex" }}>
      <Drawer variant="permanent" open={open}>
        <DrawerHeader>
          <IconButton onClick={() => setOpen(!open)}>
            {!open ? <ChevronRightIcon /> : <ChevronLeftIcon />}
          </IconButton>
        </DrawerHeader>
        <Divider />
        <List>
          <ListItem key={"filter"} disablePadding sx={{ display: "block" }}>
            {expanded["filter"] ? (
              <OpenMenuFilterItem
                open={open}
                expanded={expanded}
                setExpand={() => setExpanded(null, "filter")}
              />
            ) : (
              <FilterListItem
                isDrawerOpen={open}
                setExpand={() => setExpanded(null, "filter")}
              />
            )}
          </ListItem>
          <ListItem key={"view"} disablePadding sx={{ display: "block" }}>
            <ViewListItem isDrawerOpen={open} />
          </ListItem>
          <ListItem key={"settings"} disablePadding sx={{ display: "block" }}>
            {expanded["settings"] ? (
              <OpenMenuSettingsItem
                open={open}
                expanded={expanded}
                setExpand={() => setExpanded(null, "settings")}
              />
            ) : (
              <SettingsListItem
                isDrawerOpen={open}
                setExpand={() => setExpanded(null, "settings")}
              />
            )}
          </ListItem>
        </List>
        <Divider />
        <List>
          <ListItem key={"Help"} disablePadding sx={{ display: "block" }}>
            <HelpListItem isDrawerOpen={open} />
          </ListItem>
        </List>
      </Drawer>
    </Box>
  );
};

const HelpListItem = ({ isDrawerOpen }) => (
  <ListItemWrapper Icon={() => <HelpIcon />} isDrawerOpen={isDrawerOpen}>
    <ListItemText primary={"Help"} sx={{ opacity: isDrawerOpen ? 1 : 0 }} />
  </ListItemWrapper>
);

const ViewListItem = ({ isDrawerOpen }) => (
  <ListItemWrapper Icon={() => <ViewCompactIcon />} isDrawerOpen={isDrawerOpen}>
    <ListItemText primary={"View"} sx={{ opacity: isDrawerOpen ? 1 : 0 }} />
  </ListItemWrapper>
);

const FilterListItem = ({ isDrawerOpen, setExpand }) => (
  <ListItemWrapper
    Icon={() => <FilterListIcon onClick={() => setExpand(null)} />}
    isDrawerOpen={isDrawerOpen}
  >
    <ListItemText
      onClick={() => setExpand(null)}
      primary={"Filter Data"}
      sx={{ opacity: isDrawerOpen ? 1 : 0 }}
    />
  </ListItemWrapper>
);
const SettingsListItem = ({ isDrawerOpen, setExpand }) => (
  <ListItemWrapper
    Icon={() => (
      <SettingsIcon onClick={() => setExpand(null)} key={"settings-icon"} />
    )}
    isDrawerOpen={isDrawerOpen}
  >
    <ListItemText
      key={"settings-title"}
      onClick={() => setExpand(null, "settings")}
      primary={"Settings"}
      sx={{ opacity: isDrawerOpen ? 1 : 0 }}
    />
  </ListItemWrapper>
);

export default SideDrawer;
