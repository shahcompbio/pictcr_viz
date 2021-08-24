import React from "react";

import * as d3 from "d3";
import { useD3, sortAlphanumeric } from "@shahlab/planetarium";
import _ from "lodash";

import { Layout } from "@shahlab/planetarium";

import { CONSTANTS, INFO } from "../config";

const Doughnut = ({ data, colors, width, height, subsetParam, type }) => {
  const totalCount = data.length;
  const subsetValues = _.uniq(data.map((datum) => datum[subsetParam])).sort(
    sortAlphanumeric
  );
  //console.log(data);
  const allSubsets = _.groupBy(data, (datum) => datum[subsetParam]);

  const subsetCounts = Object.keys(allSubsets)
    .sort((a, b) => allSubsets[b].length - allSubsets[a].length)
    .slice(0, 10);

  const colorScale = d3
    .scaleOrdinal()
    .domain(subsetCounts)
    .range(colors.slice(0, Math.min(subsetCounts.length, colors.length)));

  const subsets = subsetCounts.map((d) => ({
    key: d,
    value: allSubsets[d],
  }));

  const radius = Math.min(width, height) / 2;
  const arc = d3
    .arc()
    .innerRadius(radius * 0.3)
    .outerRadius(radius * 0.65);

  const arcLabel = d3
    .arc()
    .innerRadius(radius * 0.9)
    .outerRadius(radius * 0.9);

  const drawArea = (svg) => {
    const pie = d3
      .pie()
      .value((d) => {
        console.log(d);
        return d["value"].length;
      })
      .sort(null);
    const arcs = pie(subsets);

    const path = svg
      .append("g")
      .attr("transform", "translate(" + width / 2 + "," + height / 2 + ")")
      .attr("stroke", "white")
      .selectAll("path")
      .data(arcs);

    path
      .join("path")
      .attr("fill", (d) => colorScale(d["data"].key))
      .attr("d", arc)
      .on("mouseover", function (d) {
        d3.select(this)
          .transition()
          .duration(500)
          .attr("transform", function (d) {
            var x;
            var y;
            if (d.data._translate) {
              x = d.data._translate.x;
              y = d.data._translate.y;
            } else if (!d.data._expanded) {
              d.data._expanded = true;
              var a =
                d.startAngle + (d.endAngle - d.startAngle) / 2 - Math.PI / 2;
              x = Math.cos(a) * 10;
              y = Math.sin(a) * 10;
              d.data._translate = { x: x, y: y };
            }

            return "translate(" + x + "," + y + ")";
          });
        d3.select("#label-" + d.data.key).style("font-weight", "bold");
      })
      .on("mouseout", function (d) {
        d3.select(this)
          .transition()
          .delay(200)
          .duration(500)
          .attr("transform", function (d) {
            d.data._expanded = false;
          });
        d3.select("#label-" + d.data.key).style("font-weight", "normal");
      });

    svg
      .append("g")
      .attr("transform", "translate(" + width / 2 + "," + height / 2 + ")")
      .attr("font-family", "MyFontLight")
      .attr("font-size", 12)
      .attr("text-anchor", "middle")
      .selectAll("text")
      .data(arcs)
      .join("text")
      .attr("id", (d) => "label-" + d.data.key)
      .attr("transform", (d, i) => {
        return `translate(${arcLabel.centroid(d, i)})`;
      })
      .call((text) =>
        text
          .append("tspan")
          .attr("y", "-0.4em")
          .attr("font-weight", "bold")
          .text((d) => d.data.name)
      )
      .call((text) =>
        text
          .filter((d) => d.endAngle - d.startAngle > 0.25)
          .append("tspan")

          .attr("fill-opacity", 0.7)
          .text(
            (d) =>
              d.data.key.toLocaleString() +
              " (" +
              d3.format(".0%")(d.data.value.length / totalCount) +
              ")"
          )
      );
  };

  const ref = useD3(
    (svg) => {
      drawArea(svg);
    },
    width,
    height,
    [data]
  );

  return (
    <Layout title={INFO[type]["title"]} infoText={INFO[type]["text"]}>
      <svg ref={ref} />
    </Layout>
  );
};

export default Doughnut;
