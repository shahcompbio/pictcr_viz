import React, { useState, useMemo } from "react";
import * as d3 from "d3";
import { CONSTANTS, INFO } from "../config";
import { Layout } from "@shahlab/planetarium";

import Button from "@material-ui/core/Button";
import ClearIcon from "@material-ui/icons/Clear";
import TextField from "@material-ui/core/TextField";
import Grid from "@material-ui/core/Grid";
import GetAppIcon from "@material-ui/icons/GetApp";
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

  const subHeaderComponentMemo = useMemo(() => {
    const handleClear = () => {
      if (filterText) {
        setFilterText("");
      }
    };

    return (
      <FilterComponent
        onFilter={(e) => setFilterText(e.target.value)}
        onClear={handleClear}
        filterText={filterText}
        data={data}
      />
    );
  }, [filterText, data]);
  return (
    <Layout title={INFO["TABLE"]["title"]} infoText={INFO["TABLE"]["text"]}>
      <Grid
        item
        style={{
          overflowY: "hidden",
          width: chartDim["width"],
          height: chartDim["height"],
        }}
      >
        <DataTable
          subHeader
          fixedHeader
          dense
          noHeader
          defaultSortAsc
          overflowY
          subHeaderComponent={subHeaderComponentMemo}
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
      placeholder="Filter By Gene"
      aria-label="Search Input"
      value={filterText}
      onChange={onFilter}
      InputProps={{
        endAdornment: (
          <Button
            label="Clear"
            color="primary"
            variant="outlined"
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
    <Button
      variant="outlined"
      label="Download"
      id="tsv-download"
      onClick={() => {
        const dataSource = new Blob([d3Dsv.tsvFormat(data)], {
          type: "text/tsv",
        });
        const tsvURL = window.URL.createObjectURL(dataSource);
        const tempLink = document.createElement("a");
        tempLink.href = tsvURL;
        tempLink.setAttribute("download", "filename.tsv");
        tempLink.click();
      }}
      color="secondary"
      style={{ marginLeft: 15, right: 0, position: "absolute" }}
    >
      <GetAppIcon />
    </Button>
  </div>
);
export default DEGTable;
