import React from "react";

import App from "./VDJ";
import fetchFileData from "./data/api";

const DevApp = () => {
  const data = fetchFileData();

  return Object.keys(data).length === 0 ? null : <App data={data} />;
};

const ProdApp = () => <App data={window.isablData} />;

export default process.env.NODE_ENV === "development" ? DevApp : ProdApp;
