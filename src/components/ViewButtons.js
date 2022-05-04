import Button from "@mui/material/Button";

import { styled } from "@mui/material/styles";

const ViewButtons = ({ view, setView }) => {
  const MyButton = styled(Button)(({ variant }) => ({
    ...(variant === "contained" && {
      width: "130px",
      backgroundColor: "#c2ccc2",
      marginRight: "10px",
      fontWeight: "bold",
      boxShadow: 0,
      "&:hover": {
        backgroundColor: "#fff",
        color: "#666c66",
        border: "#666c66 1px solid",
      },
    }),
    ...(variant === "outlined" && {
      width: "100px",
      backgroundColor: "#fff",
      marginRight: "10px",
      fontWeight: "bold",
      color: "#666c66",
      border: "#666c66 1px solid",
      boxShadow: 0,
      "&:hover": {
        backgroundColor: "#edf0ed",
        color: "#666c66",
        border: "#666c66 1px solid",
      },
    }),
  }));

  return (
    <div style={{ marginBottom: "20px", marginLeft: "10px" }}>
      <MyButton
        variant={view === "1" ? "contained" : "outlined"}
        disableElevation
        onClick={() => setView("1")}
      >
        View 1
      </MyButton>
      <MyButton
        variant={view === "2" ? "contained" : "outlined"}
        disableElevation
        onClick={() => setView("2")}
      >
        View 2
      </MyButton>
    </div>
  );
};
export default ViewButtons;
