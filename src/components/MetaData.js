import React from "react";

import Grid from "@mui/material/Grid";
import Typography from "@mui/material/Typography";
import Paper from "@mui/material/Paper";
import ListItemText from "@mui/material/ListItemText";
import Divider from "@mui/material/Divider";

import MenuList from "@mui/material/MenuList";
import MenuItem from "@mui/material/MenuItem";
import ListItemIcon from "@mui/material/ListItemIcon";
import ContentCut from "@mui/icons-material/ContentCut";
import ConfirmationNumberIcon from "@mui/icons-material/ConfirmationNumber";
import CircleIcon from "@mui/icons-material/Circle";
import BiotechIcon from "@mui/icons-material/Biotech";
import BubbleChartIcon from "@mui/icons-material/BubbleChart";
import LayersClearIcon from "@mui/icons-material/LayersClear";
import makeStyles from "@mui/styles/makeStyles";

const useStyles = makeStyles({
  rootWithMarginTop: {
    backgroundColor: "none",
    paddingLeft: 17,
    paddingRight: 6,
    paddingBottom: 10,
    minWidth: 275,
    marginTop: 15,
    marginLeft: 15,
    marginBottom: 0,
  },
});

const MetaData = ({ sample, hasSelection, totalCount, setHighlight }) => {
  const classes = useStyles();
  return (
    <Grid
      container
      direction="column"
      justifyContent="flex-start"
      alignItems="stretch"
    >
      <Grid item>
        <Header classes={classes} sample={sample} totalCount={totalCount} />
        <Paper
          sx={{ width: "100%", maxWidth: "90%", margin: 1, marginLeft: 8 }}
        >
          <MenuList style={{ padding: 0 }}>
            <MenuItem
              disabled={!hasSelection}
              onClick={() => {
                setHighlight();
              }}
              style={{ height: "100%" }}
            >
              <ListItemIcon>
                <LayersClearIcon fontSize="small" />
              </ListItemIcon>
              <ListItemText>Clear</ListItemText>
              <Typography variant="body2" color="text.secondary">
                ⌘X
              </Typography>
            </MenuItem>
          </MenuList>
        </Paper>
      </Grid>
    </Grid>
  );
};

const Header = ({ classes, sample, totalCount }) => (
  <div className={classes.rootWithMarginTop}>
    <Paper sx={{ width: "100%", maxWidth: "100%" }}>
      <MenuList>
        <MenuItem>
          <ListItemIcon>
            <BiotechIcon fontSize="small" />
          </ListItemIcon>
          <ListItemText>{sample}</ListItemText>
        </MenuItem>
        <MenuItem>
          <ListItemIcon>
            <BubbleChartIcon fontSize="small" />
          </ListItemIcon>
          <ListItemText>{totalCount} cells</ListItemText>
        </MenuItem>
      </MenuList>
    </Paper>
  </div>
);

export default MetaData;
