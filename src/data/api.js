import { useEffect, useState } from "react";

import * as d3 from "d3";

// import metadataSource from "./test/metadata.tsv";
// import probabilitiesSource from "./test/probabilities.tsv";
// import degsSource from "./test/deg.tsv";

// import metadataSource from "./large/metadata.tsv";
// import probabilitiesSource from "./large/probabilities.tsv";
// import degsSource from "./large/degs.tsv";

// import metadataSource from "./test2/metadata.tsv";
// import probabilitiesSource from "./test2/probabilities.tsv";
// import degsSource from "./test2/degs.tsv";

// import metadataSource from "./metadata.tsv";
// import probabilitiesSource from "./probabilities.tsv";
// import degsSource from "./degs.tsv";

// import metadataSource from "./ranked/metadata3.tsv";
// // import probabilitiesSource from "./ranked/probabilities.tsv";
// import degsSource from "./ranked/degs3.tsv";

import metadataSource from "./filters/metadata.tsv";
import degsSource from "./filters/degs.tsv";

import filters from "./filters/filters.json";

const useFetchData = () => {
  const [data, setData] = useState({});

  useEffect(() => {
    Promise.all([d3.tsv(metadataSource), d3.tsv(degsSource)]).then((data) => {
      setData({ metadata: data[0], degs: data[1], filters });
    });
  }, []);

  return data;
};

export default useFetchData;
