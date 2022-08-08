import ListItemButton from "@mui/material/ListItemButton";
import ListItemIcon from "@mui/material/ListItemIcon";
import ChevronLeftIcon from "@mui/icons-material/ChevronLeft";
import IconButton from "@mui/material/IconButton";

const ListItemWrapper = ({ Icon, isDrawerOpen, children }) => (
  <ListItemButton
    sx={{
      minHeight: 48,
      justifyContent: isDrawerOpen ? "initial" : "center",
      px: 2.5,
    }}
  >
    <ListItemIcon
      sx={{
        minWidth: 0,
        mr: isDrawerOpen ? 3 : "auto",
        justifyContent: "center",
      }}
    >
      <Icon />
    </ListItemIcon>
    {children}
    {isDrawerOpen && (
      <IconButton>
        <ChevronLeftIcon />
      </IconButton>
    )}
  </ListItemButton>
);
export default ListItemWrapper;
