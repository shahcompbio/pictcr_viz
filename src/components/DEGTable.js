import React, { useState } from "react";
import * as d3 from "d3";
import { CONSTANTS, INFO } from "../config";
import { Layout, SearchIcon, DownloadIcon } from "@shahlab/planetarium";

import Button from "@mui/material/Button";
import ClearIcon from "@mui/icons-material/Clear";
import TextField from "@mui/material/TextField";
import Grid from "@mui/material/Grid";
import Switch from "@mui/material/Switch";
import QueryStatsIcon from "@mui/icons-material/QueryStats";
import Avatar from "@mui/material/Avatar";
import Tooltip from "@mui/material/Tooltip";

import * as d3Dsv from "d3-dsv";
import DataTable from "react-data-table-component";
const formatCols = ["adj_pval", "log_fc"];
const formatDecimal = [(num) => num.toExponential(2), d3.format(",.4f")];

const formatHeader = {
  gene: "Gene",
  adj_pval: `Adjusted p-value`,
  log_fc: "log fold change",
  subtype: "Phenotype",
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

const DEGTable = ({ data, chartDim, selectedSubtype }) => {
  const [filterText, setFilterText] = useState("");
  const { subtypeParam } = CONSTANTS;

  const columns = Object.keys(data[0]);
  const dataSource = selectedSubtype
    ? data.filter((row) => row[subtypeParam] === selectedSubtype)
    : data;
  const filteredItems = dataSource.filter(
    (item) =>
      item["gene"] &&
      item["gene"].toLowerCase().includes(filterText.toLowerCase())
  );
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
        <Avatar
          style={{ width: 30, height: 30, marginRight: 10 }}
          variant="rounded"
        >
          <Tooltip title={"T-Test"} arrow>
            <QueryStatsIcon sx={{ fontSize: 30 }} />
          </Tooltip>
        </Avatar>,
      ]}
    >
      <Grid
        item
        style={{
          overflowY: "hidden",
          width: chartDim["width"],
          height: chartDim["height"],
          marginTop: 10,
        }}
      >
        <DataTable
          subHeader
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
          data={filteredItems}
        />
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
export default DEGTable;
