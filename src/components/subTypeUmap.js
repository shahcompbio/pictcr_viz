import React, { useState } from "react";
import _ from "lodash";
import { INFO } from "../config.js";

import { Layout, UMAP, Select } from "@shahlab/planetarium";

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
}) => {
  const [subset, setSubset] = useState(subsetParam);

  return (
    <Layout
      title={INFO["SUBTYPEUMAP"]["title"]}
      infoText={INFO["SUBTYPEUMAP"]["text"]}
    >
      <Select
        options={options}
        value={subset}
        title={"Color"}
        onSelect={setSubset}
      />
      <UMAP
        width={700}
        height={600}
        data={data}
        xParam={xParam}
        yParam={yParam}
        subsetParam={subset}
        idParam={idParam}
        colorScale={subset === subsetParam ? colorScale : null}
        onLasso={onLasso}
        onLegendClick={onLegendClick}
        disable={disable}
        highlightIDs={highlightIDs}
      />
    </Layout>
  );
};

export default PhenotypeUMAP;
