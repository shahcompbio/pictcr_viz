import React, { useState, useEffect } from "react";
import * as d3 from "d3";
import _ from "lodash";
import { CONSTANTS, INFO } from "../config";
import { Layout, SearchIcon, DownloadIcon } from "@shahlab/planetarium";

import Button from "@mui/material/Button";
import ClearIcon from "@mui/icons-material/Clear";
import TextField from "@mui/material/TextField";
import Grid from "@mui/material/Grid";
import CircularProgress from "@mui/material/CircularProgress";
import Box from "@mui/material/Box";

import * as d3Dsv from "d3-dsv";
import DataTable from "react-data-table-component";
const formatCols = ["adj_pval", "log_fc"];
const formatDecimal = [(num) => num.toExponential(2), d3.format(",.4f")];
const url = "http://127.0.0.1:5000";
const dir = "Users/bojilovv/Downloads/hacohen_viz.h5ad/";

//log_fc: "log fold change",
//subtype: "Phenotype",
const formatHeader = {
  cell_id: "Cell Name",
  cell_idx: "Cell ID",
  gene_idx: "Gene ID",
  gene_id: "Gene",
  expression: `Adjusted p-value`,
};
const customStyles = {
  rows: {
    style: {
      fontFamily: "Noto Sans",
    },
  },
  headCells: {
    style: {
      fontFamily: "Noto Sans",
      fontWeight: "bold",
    },
  },
  cells: {
    style: {
      fontFamily: "Noto Sans",
    },
  },
};

const DEGTable = ({ chartDim, selectedSubtype, selection }) => {
  const [filterText, setFilterText] = useState("");
  const { subtypeParam } = CONSTANTS;
  const [error, setError] = useState(null);
  const [isLoaded, setIsLoaded] = useState(false);
  const [data, setData] = useState([]);
  //  var fd = new FormData();
  //  fd.append("cells", selection.join(","));

  useEffect(() => {
    fetch("http://127.0.0.1:5000/testing/", {
      credentials: "include",
    })
      .then((res) => res.json())
      .then(
        (result) => {
          if (result.length > 0) {
            setIsLoaded(true);
            setData(result);
            return {};
          } else {
            return fetch(`${url}/load/${dir}`, {
              credentials: "include",
            });
          }
        },
        (error) => {
          setIsLoaded(true);
          setError(error);
        }
      )
      .then((result) => (_.isEmpty(result) ? {} : result.json()))
      .then((response) => {
        if (response.length) {
          setIsLoaded(true);
          setData(response);
        }
      });
  }, []);

  const columns = data && data.length > 0 ? Object.keys(data[0]) : [];
  /*  const dataSource = selectedSubtype
      ? data.filter((row) => row[subtypeParam] === selectedSubtype)
      : data
    ;

  const filteredItems = dataSource.filter(
    (item) =>
      item["gene"] &&
      item["gene"].toLowerCase().includes(filterText.toLowerCase())
  );*/
  const handleClear = () => {
    if (filterText) {
      setFilterText("");
    }
  };

  const download = () => {
    const dataSource = new Blob([d3Dsv.tsvFormat(data)], {
      type: "text/tsv",
    });
    const tsvURL = window.URL.createObjectURL(dataSource);
    const tempLink = document.createElement("a");
    tempLink.href = tsvURL;
    tempLink.setAttribute("download", "filename.tsv");
    tempLink.click();
  };

  return (
    <Layout
      title={INFO["TABLE"]["title"]}
      infoText={INFO["TABLE"]["text"]}
      addIcon={[
        <SearchIcon>
          <FilterComponent
            onFilter={(e) => setFilterText(e.target.value)}
            onClear={handleClear}
            filterText={filterText}
            data={data}
          />
        </SearchIcon>,
        <DownloadIcon download={download} />,
      ]}
    >
      <Grid
        item
        style={{
          overflowY: "hidden",
          width: chartDim["width"],
          height: chartDim["height"],
          padding: 15,
        }}
      >
        {isLoaded ? (
          <DataTable
            fixedHeader
            dense
            noHeader
            defaultSortAsc
            overflowY
            customStyles={customStyles}
            compact
            columns={columns.map((col) => {
              const formatIndex = formatCols.indexOf(col);

              return formatIndex !== -1
                ? {
                    name: <b>{formatHeader[col]}</b>,
                    selector: col,
                    sortable: true,
                    right: true,
                    cell: (row) => (
                      <span>
                        {formatDecimal[formatIndex](parseFloat(row[col]))}
                      </span>
                    ),
                  }
                : {
                    name: <b>{formatHeader[col]}</b>,
                    selector: col,
                    sortable: true,
                    right: true,
                  };
            })}
            data={data}
          />
        ) : (
          <LoadingComponent />
        )}
      </Grid>
    </Layout>
  );
};

const FilterComponent = ({ filterText, onFilter, onClear, data }) => (
  <div style={{ display: "flex", width: "100%" }}>
    <TextField
      color="primary"
      type="text"
      id="searchGenes"
      size="small"
      autoFocus
      margin="dense"
      placeholder="Filter By Gene"
      aria-label="Search Input"
      value={filterText}
      onChange={onFilter}
      InputProps={{
        endAdornment: (
          <Button
            label="Clear"
            color="primary"
            size="small"
            onClick={onClear}
            style={{
              marginLeft: 15,
              marginBottom: 5,
            }}
          >
            <ClearIcon />
          </Button>
        ),
      }}
    />
  </div>
);
const LoadingComponent = () => (
  <Box sx={{ display: "flex" }}>
    <CircularProgress />
  </Box>
);
export default DEGTable;
