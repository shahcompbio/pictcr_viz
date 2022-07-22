import _ from "lodash";
import * as d3 from "d3";

import { CONSTANTS, CLONOTYPE_COLORS } from "./config";

const parseClonotypeData = (data, stats) => {
  const { clonotypeParam } = CONSTANTS;
  console.log(data);
  console.log(stats);
  if (data && data.length !== 0 && stats) {
    const clonotypeCounts = _.countBy(
      data.filter((datum) => datum[clonotypeParam] !== "None"),
      clonotypeParam
    );

    console.log(stats);

    const clonotypeLabels = Object.keys(stats)
      .sort((a, b) => stats[b] - stats[a])
      .slice(0, 10)
      .map((value, index) => ({
        value,
        label: `Clone ${value}`,
        color: CLONOTYPE_COLORS[index],
      }));

    console.log(clonotypeLabels);
    const cloneColorScale = d3
      .scaleOrdinal()
      .domain(clonotypeLabels.map((label) => label["value"]))
      .range(
        CLONOTYPE_COLORS.slice(
          0,
          Math.min(clonotypeLabels.length, CLONOTYPE_COLORS.length)
        )
      )
      .unknown("[0.91,0.91,0.91,1.0]");

    return {
      clonotypeCounts: clonotypeCounts,
      clonotypeLabels: clonotypeLabels,
      cloneColorScale: cloneColorScale,
      clonotypeDataIsLoaded: true,
    };
  } else {
    return {
      clonotypeCounts: null,
      clonotypeLabels: null,
      cloneColorScale: null,
      clonotypeDataIsLoaded: false,
    };
  }
};
export default parseClonotypeData;