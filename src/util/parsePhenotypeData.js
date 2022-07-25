import _ from "lodash";
import * as d3 from "d3";

import { useData } from "../provider/dataContext";

import { PHENOTYPE_COLORS } from "../config";

const parsePhenotypeData = (data, subset) => {
  const phenotypeValues = Object.keys(_.groupBy(data, subset)).sort();
  const phenotypeColorScale = d3
    .scaleOrdinal()
    .domain(phenotypeValues)
    .range(
      PHENOTYPE_COLORS.slice(
        0,
        Math.min(phenotypeValues.length, PHENOTYPE_COLORS.length)
      )
    )
    .unknown("#e8e8e8");

  return {
    phenotypeColorScale: phenotypeColorScale,
    phenotypeValues: phenotypeValues,
  };
};

export default parsePhenotypeData;
