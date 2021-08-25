import React from "react";
import _ from "lodash";
import { INFO } from "../config.js";

import { Layout, UMAP } from "@shahlab/planetarium";

import { CONSTANTS } from "../config";

const DataWrapper = ({
  data,
  chartDim,
  selectedSubtype,
  hoveredSubtype,
  setSelectedSubtype,
}) => {
  const { xParam, yParam, subtypeParam } = CONSTANTS;

  const subset = hoveredSubtype || selectedSubtype;

  const highlightIDs =
    subset === null
      ? null
      : data
          .filter((datum) => datum[subtypeParam] === subset)
          .map((datum) => datum["cell_id"]);
  return (
    <Layout
      title={INFO["SUBTYPEUMAP"]["title"]}
      infoText={INFO["SUBTYPEUMAP"]["text"]}
    >
      <UMAP
        width={chartDim["width"]}
        height={chartDim["height"]}
        data={data}
        highlightIDs={highlightIDs}
        xParam={xParam}
        yParam={yParam}
        subsetParam={subtypeParam}
        idParam={"cell_id"}
        onLegendHover={(value) => {
          setSelectedSubtype({
            hover: value,
          });
        }}
        fontFamily={{
          regular: "MyFontLight",
          bold: "MyFontBold",
          labelOffset: 3,
        }}
        onLegendClick={(value) => {
          setSelectedSubtype({ selected: value });
        }}
      />
    </Layout>
  );
};

export default DataWrapper;
