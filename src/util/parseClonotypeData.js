import _ from "lodash";
import * as d3 from "d3";

import { CONSTANTS, CLONOTYPE_COLORS, CLONOTYPE_COLORS_2 } from "../config";

const parseClonotypeData = (data, stats, filterNA) => {
  const { clonotypeParam } = CONSTANTS;

  if (data && data.length !== 0 && stats) {
    const clonotypeCounts = _.countBy(
      filterNA ? data.filter((datum) => datum[clonotypeParam] !== "NA") : data,
      clonotypeParam
    );

    const clonotypeLabels = Object.keys(stats)
      .filter((d) => (filterNA ? d !== "NA" : true))
      .sort((a, b) => stats[b] - stats[a])
      .slice(0, 10)
      .map((value, index) => ({
        value,
        label: `Clone ${value}`,
        color: CLONOTYPE_COLORS[index],
      }));

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

    const cloneHEXColorScale = d3
      .scaleOrdinal()
      .domain(clonotypeLabels.map((label) => label["value"]))
      .range(
        CLONOTYPE_COLORS_2.slice(
          0,
          Math.min(clonotypeLabels.length, CLONOTYPE_COLORS_2.length)
        )
      );

    return {
      clonotypeCounts: clonotypeCounts,
      clonotypeLabels: clonotypeLabels,
      cloneColorScale: cloneColorScale,
      cloneHEXColorScale: cloneHEXColorScale,
      clonotypeDataIsLoaded: true,
    };
  } else {
    return {
      clonotypeCounts: null,
      clonotypeLabels: null,
      cloneColorScale: null,
      cloneHEXColorScale: null,
      clonotypeDataIsLoaded: false,
    };
  }
};
export default parseClonotypeData;
