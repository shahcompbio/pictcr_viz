import { createTheme, adaptV4Theme } from "@mui/material/styles";
import Merriweather from "./Merriweather-Regular.ttf";
export const theme = createTheme(
  adaptV4Theme({
    props: {
      MuiSvgIcon: {
        htmlColor: "#aa0011 !important",
      },
    },
    typography: {
      fontFamily: ["Helvetica"].join(","),
    },
    palette: {
      primary: {
        main: "#95d2dc",
        dark: "#618ba0",
      },
      secondary: {
        main: "#f1c023",
      },
      error: {
        main: "#BC4746",
      },
      background: {
        default: "white",
      },
      overrides: {
        MuiFab: {
          root: {
            boxShadow: "none",
          },
        },
      },
    },

    spacing: 4,
  })
);
