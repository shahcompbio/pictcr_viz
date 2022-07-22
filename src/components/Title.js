import React from "react";
import styled from "styled-components";
import Typography from "@mui/material/Typography";

export const StyledTitle = styled(Typography)({
  fontWeight: "bold",
  fontSize: 28,
  color: "rgba(106,141,149, 0.8)",
  textShadow: "0px 4px 4px #fff, 0 0 0 #000, 0px 4px 4px #fff",
});
export const StyledInfo = styled(Typography)({
  fontWeight: "bold",
  fontSize: 20,
  marginLeft: 10,
  marginTop: 5,
  color: "rgba(106,141,149, 0.8)",
  textShadow: "0px 4px 4px #fff, 0 0 0 #000, 0px 4px 4px #fff",
});
const Title = ({ title }) => {
  return (
    <div
      id="scrollbar-wrapper"
      style={{
        //  position: "fixed",
        width: "300px",
        height: "100px",
        display: "flex",
        //marginTop: 50,
        //  background: "white",
      }}
    >
      <StyledTitle>{title}</StyledTitle>
      <StyledInfo>ⓘ</StyledInfo>
    </div>
  );
};

export default Title;
