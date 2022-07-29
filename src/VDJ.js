import React, { useState } from "react";

import View1 from "./View1.js";
import View2 from "./View2.js";
import { theme } from "./theme/theme";
import { ThemeProvider, StyledEngineProvider } from "@mui/material/styles";
import CssBaseline from "@mui/material/CssBaseline";
import SideDrawer from "./components/SideDrawer.js";
import Box from "@mui/material/Box";

import { useData } from "./provider/dataContext";

export const VDJ = () => {
  const [{ metadata, filters, stats }] = useData();
  const [view, setView] = useState("1");

  return (
    <StyledEngineProvider injectFirst>
      <ThemeProvider theme={theme}>
        <Box sx={{ display: "flex" }}>
          <CssBaseline />

          {filters && <SideDrawer view={view} setView={setView} />}
          {metadata &&
            (view === "1" ? (
              <View1
                metadata={metadata}
                filters={filters}
                view={view}
                setView={setView}
              />
            ) : (
              <View2
                metadata={metadata}
                filters={filters}
                view={view}
                setView={setView}
              />
            ))}
        </Box>
      </ThemeProvider>
    </StyledEngineProvider>
  );
};

export default VDJ;
