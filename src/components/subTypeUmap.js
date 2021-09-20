import React, { useEffect, useState } from "react";
import _ from "lodash";
import { INFO } from "../config.js";

import { Layout, UMAP } from "@shahlab/planetarium";
import { CircularProgress } from "@material-ui/core";

const PhenotypeUMAP = ({
  data,
  xParam,
  yParam,
  subsetParam,
  idParam,
  colorScale,
  onLasso,
  onLegendClick,
  disable,
  highlightIDs,
  options,
  Select,
}) => {
  return (
    <Layout
      title={INFO["SUBTYPEUMAP"]["title"]}
      infoText={INFO["SUBTYPEUMAP"]["text"]}
      addIcon={Select}
    >
      <UMAP
        width={700}
        height={600}
        data={data}
        xParam={xParam}
        yParam={yParam}
        subsetParam={subsetParam}
        idParam={idParam}
        colorScale={colorScale}
        onLasso={onLasso}
        fontFamily={"Noto Sans"}
        onLegendClick={onLegendClick}
        disable={disable}
        highlightIDs={highlightIDs}
      />
    </Layout>
  );
};
/*      <div
        style={{
          position: "absolute",
          width: 700,
          height: 600,
          background: "#f8f8f",
          marginTop: 250,
        }}
      >
        <CircularProgress size={150} thickness={7} />
      </div>
      <div
        style={{
          filter: "blur(3px)",
          opacity: 0.3,
        }}
      >
        {componentUmap}
      </div>*/
export default PhenotypeUMAP;
