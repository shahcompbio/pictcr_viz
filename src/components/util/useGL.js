import { useRef, useCallback } from "react";
import * as d3 from "d3";

export function useGL(width, height, ids = ["umapCanvas"]) {
  const ref = useRef(null);
  const setRef = useCallback(
    (node) => {
      if (ref.current) {
      }

      if (node) {
        ids.map((id) => {
          const canvas = node.appendChild(document.createElement("canvas"));
          canvas.id = id;
          const context = canvas.getContext("webgl");
          let scale = window.devicePixelRatio;
          canvas.style.width = width + "px";
          canvas.style.height = height + "px";
          canvas.width = width * scale;
          canvas.height = height * scale;
        });
      }

      ref.current = node;
    },
    [width, height]
  );

  return [setRef];
}
