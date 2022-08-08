import React, {
  useState,
  useEffect,
  useRef,
  createContext,
  useReducer,
} from "react";

import { CONSTANTS } from "../config";
const { subtypeParam, clonotypeParam } = CONSTANTS;

const initialState = {
  metadataOG: null,
  metadata: null,
  filters: null,

  selectClone: null,
  selectIDs: null,
  activeGraph: null,
  subset: subtypeParam,
  selectSubset: null,
  selectFilters: null,
  selectPhenotype: null,

  settings: { filterNA: true },

  inputRefMapping: {
    1: [
      { plotName: "Sunburst", refName: "sunburstRef" },
      { plotName: "Phenotype UMAP", refName: "phenotypeRef" },
      { plotName: "Doughnut", refName: "doughnutRef" },
      //{ plotName: "Clonotype UMAP", refName: "cloneumapRef" },
      { plotName: "Heatmap", refName: "heatmapRef" },
      { plotName: "Clonotype Expansion", refName: "clonotypeExpansionRef" },
      { plotName: "Probability Histogram", refName: "probabilityHistogramRef" },
      { plotName: "Ranked Order", refName: "rankedOrderRef" },
    ],
    2: [
      { plotName: "Sankey", refName: "sankeyRef" },
      { plotName: "Phenotype UMAP", refName: "phenotypeRef" },
      { plotName: "Heatmap", refName: "heatmapRef" },
      { plotName: "Clonotype Expansion", refName: "clonotypeExpansionRef" },
      { plotName: "Probability Histogram", refName: "probabilityHistogramRef" },
      { plotName: "Ranked Order", refName: "rankedOrderRef" },
    ],
  },
};
const dataReducer = (state, action) => {
  console.log(state);
  switch (action.type) {
    case "initialSetup": {
      const initValues = action.values.reduce(
        (final, d) => ({ [d.valueType]: d.value, ...final }),
        {}
      );
      return {
        ...state,
        ...initValues,
      };
    }
    case "setFilterNA": {
      return {
        ...state,
        settings: { filterNA: !state.settings.filterNA },
        metadata: !state.settings.filterNA
          ? state.metadataOG.filter((d) => d[clonotypeParam] !== "NA")
          : state.metadataOG,
      };
    }
    case "setSettings": {
      return { ...state, settings: { ...state.settings, ...action.value } };
    }

    case "setSelectFilters": {
      return { ...state, selectFilters: action.value };
    }
    case "setSubset": {
      return { ...state, subset: action.value };
    }
    case "setSelectSubset": {
      return {
        ...state,
        selectSubset: action.value,
      };
    }
    default:
      return state;
  }
};
const actions = (action, state) => {
  switch (action.type) {
    case "set-" + action.valueType: {
      return {
        [action.valueType]: action.value,
        ...state,
      };
    }
    default:
      return state;
  }
};
export { initialState };
export default dataReducer;
