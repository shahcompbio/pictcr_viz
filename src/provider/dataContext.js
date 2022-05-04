import React, { useReducer, useContext, createContext } from "react";
const DataContext = createContext();

function DataProvider({ reducer, children, initialState }) {
  return (
    <DataContext.Provider value={useReducer(reducer, initialState)}>
      {children}
    </DataContext.Provider>
  );
}
function useData() {
  return useContext(DataContext);
}
export { DataProvider, useData };
