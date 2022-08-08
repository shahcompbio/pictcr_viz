import { useEffect, useState } from "react";

import { CONSTANTS } from "../config";
const url = process.env.HOST ? process.env.HOST : "http://127.0.0.1:5000";
const { clonotypeParam } = CONSTANTS;
//const url = "https://spectrum-staging.shahlab.mskcc.org";
const useFetchData = () => {
  const [data, setData] = useState({});

  useEffect(() => {
    fetch(url + "/getStats/", {
      credentials: "include",
    })
      .then((res) => res.json())
      .then((data) => {
        setData({
          ...data,
          stats: JSON.parse(data["stats"][clonotypeParam + "-stats"]),
        });
      });

    fetch(url + "/getFilterData/", {
      credentials: "include",
    })
      .then((res) => res.json())
      .then((data) => {
        setData({ ...data, filters: data["filters"] });
      });
    fetch(url + "/getStreamData/", {
      credentials: "include",
    })
      .then((response) => response.body)
      .then((rb) => {
        const reader = rb.getReader();

        return new ReadableStream({
          start(controller) {
            // each data chunk
            function push() {
              // "done" is a Boolean and value a "Uint8Array"
              reader.read().then(({ done, value }) => {
                if (done) {
                  //  console.log("done", done);
                  controller.close();
                  return;
                }
                controller.enqueue(value);
                //  console.log(done, value);
                push();
              });
            }

            push();
          },
        });
      })
      .then((stream) => {
        return new Response(stream, {
          headers: { "Content-Type": "text/html" },
        }).text();
      })
      .then((result) => {
        const metadata = result.split("&/").reduce((final, d) => {
          if (d.length !== 0) {
            return [...final, ...JSON.parse(d)];
          } else {
            return [...final];
          }
        }, []);

        setData({
          ...data,
          metadata: [...metadata.filter((d) => d[clonotypeParam] !== "NA")],
          metadataOG: [...metadata],
        });
      });
  }, []);

  return data;
};

export default useFetchData;
