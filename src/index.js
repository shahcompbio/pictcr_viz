import React from "react";
import ReactDOM from "react-dom";
import "./index.css";
import App from "./App";
import * as serviceWorker from "./serviceWorker";
import dataReducer, { initialState } from "./provider/dataReducer";
import { DataProvider } from "./provider/dataContext";
//containerId, history
window.renderpicture = (containerId) => {
  ReactDOM.render(
    <DataProvider initialState={initialState} reducer={dataReducer}>
      <App />
    </DataProvider>,
    document.getElementById(containerId)
  );
  serviceWorker.unregister();
};

window.unmountpicture = (containerId) => {
  ReactDOM.unmountComponentAtNode(document.getElementById("picture-container"));
};

if (!document.getElementById("picture-container")) {
  ReactDOM.render(
    <DataProvider initialState={initialState} reducer={dataReducer}>
      <App />
    </DataProvider>,
    document.getElementById("root")
  );
  serviceWorker.unregister();
}

// If you want your app to work offline and load faster, you can change
// unregister() to register() below. Note this comes with some pitfalls.
// Learn more about service workers: https://bit.ly/CRA-PWA
//serviceWorker.unregister();
