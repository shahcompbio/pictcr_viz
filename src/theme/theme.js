import { createTheme, adaptV4Theme } from "@mui/material/styles";

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
        default: "#F5F5F5",
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
