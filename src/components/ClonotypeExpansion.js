import React from "react";

import { StackedHorizontalBar, Layout } from "@shahlab/planetarium";
import { CONSTANTS, INFO } from "../config";
import _ from "lodash";

const ClonotypeExpansion = ({ data, width, height, highlightedRow }) => {
  const { clonotypeParam, subtypeParam } = CONSTANTS;

  const groupedSubtype = _.groupBy(data, subtypeParam);
  const subtypes = Object.keys(groupedSubtype).sort();

  const countedClonotypes = subtypes.reduce((countMap, subtype) => {
    const clonotypeCount = _.countBy(groupedSubtype[subtype], clonotypeParam);

    const countFreq = _.countBy(Object.values(clonotypeCount), (value) =>
      Math.min(value, 10)
    );

    return { ...countMap, [subtype]: countFreq };
  }, {});

  return (
    <Layout title={INFO["BARPLOT"]["title"]} infoText={INFO["BARPLOT"]["text"]}>
      <StackedHorizontalBar
        highlightedRow={highlightedRow}
        font={"MyFontRegular"}
        data={countedClonotypes}
        width={width}
        height={height}
        barLabels={Array.from(Array(10).keys()).map((value) => ({
          value: value + 1,
          label: value === 9 ? "≥10" : value + 1,
        }))}
      />
    </Layout>
  );
};

export default ClonotypeExpansion;
