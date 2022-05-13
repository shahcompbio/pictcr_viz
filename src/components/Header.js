import React from "react";
import makeStyles from '@mui/styles/makeStyles';
import AppBar from "@mui/material/AppBar";
import Toolbar from "@mui/material/Toolbar";
import IconButton from "@mui/material/IconButton";
import Typography from "@mui/material/Typography";
import HelpOutlineIcon from "@mui/icons-material/HelpOutline";
import Grid from "@mui/material/Grid";
import BookmarkIcon from "@mui/icons-material/Bookmark";
import InfoIcon from "@mui/icons-material/Info";

const useStyles = makeStyles((theme) => ({
  root: {
    height: 100,
    //    backgroundColor: "#F7F8FB",
    //  boxShadow:
    //      "2px 2px 30px 2px rgb(0 0 0 / 20%), 0px 2px 0px 0px rgb(0 0 0 / 14%), 0px 1px 10px 0px rgb(0 0 0 / 12%)",
  },
  grow: {
    flexGrow: 1,
    fill: "#F7F8FB",
  },
  menuButton: {
    marginRight: theme.spacing(2),
  },
  title: {
    fontFamily: "Noto Sans",
    display: "none",
    [theme.breakpoints.up("sm")]: {
      display: "block",
    },
  },

  inputInput: {
    padding: theme.spacing(1, 1, 1, 0),
    // vertical padding + font size from searchIcon
    paddingLeft: `calc(1em + ${theme.spacing(4)})`,
    transition: theme.transitions.create("width"),
    width: "100%",
    [theme.breakpoints.up("md")]: {
      width: "20ch",
    },
  },
  sectionDesktop: {
    display: "none",
    [theme.breakpoints.up("md")]: {
      display: "flex",
    },
  },
  sectionMobile: {
    display: "flex",
    [theme.breakpoints.up("md")]: {
      display: "none",
    },
  },
}));

const Header = () => {
  const classes = useStyles();

  return (
    <div className={classes.grow}>
      <Grid
        container
        direction="row"
        justifyContent="flex-start"
        alignItems="stretch"
        position="static"
        className={classes.root}
      >
        <Toolbar>
          <IconButton
            edge="start"
            className={classes.menuButton}
            color="inherit"
            aria-label="open drawer"
            size="large"></IconButton>
          <Typography className={classes.title} variant="h5" noWrap>
            PICTCR
          </Typography>

          <div className={classes.grow} />
          <div className={classes.sectionDesktop}>
            <IconButton
              edge="end"
              aria-label="Help"
              aria-haspopup="true"
              color="inherit"
              size="large">
              <HelpOutlineIcon />
            </IconButton>
            <IconButton
              edge="end"
              aria-label="bookmark"
              aria-haspopup="true"
              color="inherit"
              size="large">
              <BookmarkIcon />
            </IconButton>
            <IconButton aria-label="info" color="inherit" size="large">
              <InfoIcon />
            </IconButton>
          </div>
        </Toolbar>
      </Grid>
    </div>
  );
};
export default Header;
