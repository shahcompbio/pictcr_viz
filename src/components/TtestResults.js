import React, { useState, useRef, Fragment, useEffect } from "react";

import * as d3 from "d3";
import { useD3, DownloadIcon } from "@shahlab/planetarium";
import _ from "lodash";

import makeStyles from "@mui/styles/makeStyles";

import { Layout } from "@shahlab/planetarium";

import {
  Box,
  CircularProgress,
  Typography,
  Table,
  TableBody,
  TableHead,
  TableCell,
  TableRow,
  TableContainer,
} from "@mui/material";
import { ClearBox } from "./SideDrawer/Filters";

import { INFO } from "../config";
const url = process.env.HOST ? process.env.HOST : "http://localhost:5000";

const TtestResults = ({ count, cells, idParam, resetLasso, width, height }) => {
  const [data, settTestData] = useState([]);

  useEffect(() => {
    if (cells && cells.length > 0) {
      const param = cells.map((d) => d[idParam]).join(",");
      fetch(url + "/ttest/", {
        method: "POST",
        credentials: "include",
        body: JSON.stringify({ data: param }),
      })
        .then((res) => res.json())
        .then((result) => {
          if (result.data) {
            settTestData(result.data);
          }
        });
    }
  }, [cells]);

  return count ? (
    <div
      style={{
        height: height,
        width: width,
      }}
    >
      <div style={{ marginBottom: "10px" }}>
        <ClearBox
          disabled={false}
          onClick={resetLasso}
          clearText={"Clear Lasso"}
        />
      </div>
      <Typography color="textSecondary" gutterBottom>
        <b>{count}</b> cells selected
      </Typography>
      <Typography color="textSecondary" gutterBottom>
        Top genes
      </Typography>
      {count > 0 && data.length === 0 ? (
        <Box
          sx={{
            display: "flex",
            height: 250,
            width: 150,
            paddingTop: 25,
            justifyContent: "space-evenly",
          }}
        >
          <CircularProgress />
        </Box>
      ) : (
        <TableContainer style={{ maxHeight: height - 100 }}>
          <Table
            stickyHeader
            size="small"
            style={{ overflowY: "scroll", overflowX: "none" }}
          >
            <TableHead>
              <TableRow>
                <TableCell style={{ fontSize: "0.7rem" }}>Gene</TableCell>
                <TableCell style={{ fontSize: "0.7rem" }} align="right">
                  Adjusted P value
                </TableCell>
              </TableRow>
            </TableHead>
            <TableBody>
              {data.map((gene) => (
                <TableRow key={gene["gene"]} style={{ paddingBottom: "5px" }}>
                  <TableCell style={{ borderBottom: "none" }}>
                    <b>{gene["gene"]}</b>
                  </TableCell>
                  <TableCell style={{ borderBottom: "none" }}>
                    {d3.format("0.3")(gene["p"])}
                  </TableCell>
                </TableRow>
              ))}
            </TableBody>
          </Table>
        </TableContainer>
      )}
    </div>
  ) : (
    <div />
  );
};

export default TtestResults;
