import React, { useEffect } from "react";

import App from "./VDJ";
import fetchFileData from "./data/api";

import { useData } from "./provider/dataContext";

const DevApp = () => {
  const data = fetchFileData();
  const [{}, dispatch] = useData();

  useEffect(() => {
    if (Object.keys(data).length !== 0) {
      dispatch({
        type: "initialSetup",
        values: [
          {
            value: data["metadata"],
            valueType: "metadata",
          },
          {
            value: data["filters"].filter((d) => d.name !== "clone_id"),
            valueType: "filters",
          },
        ],
      });
    }
  }, [data]);

  return Object.keys(data).length === 0 ? null : <App />;
};

const ProdApp = () => <App data={window.isablData} />;

//export default Test;
export default process.env.NODE_ENV === "development" ? DevApp : ProdApp;
