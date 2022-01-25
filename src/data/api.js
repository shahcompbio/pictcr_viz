import { useEffect, useState } from "react";

import * as d3 from "d3";

const url = process.env.HOST ? process.env.HOST : "http://127.0.0.1:5000";
//const url = "https://spectrum-staging.shahlab.mskcc.org";
const useFetchData = () => {
  const [data, setData] = useState({});

  useEffect(() => {
    fetch(url + "/getData/hacohen_viz.h5ad/", {
      credentials: "include",
    })
      .then((res) => res.json())
      .then((data) => {
        setData({ metadata: data["metadata"], filters: data["filters"] });
      });
  }, []);

  return data;
};

export default useFetchData;
