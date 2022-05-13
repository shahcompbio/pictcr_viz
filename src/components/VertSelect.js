import * as React from "react";
import Box from "@mui/material/Box";
import InputLabel from "@mui/material/InputLabel";
import MenuItem from "@mui/material/MenuItem";
import FormControl from "@mui/material/FormControl";
import Select from "@mui/material/Select";

const ColorSelect = ({
  options,
  value,
  onSelect = (value) => {},
  title = "",
  width = 300,
  id = "select",
}) => {
  console.log(options);
  return (
    <FormControl sx={{ width: width }}>
      <InputLabel id={id}>{title}</InputLabel>
      <Select
        labelId="demo-simple-select-label"
        id="demo-simple-select"
        value={value}
        label="Age"
        onChange={(e) => {
          onSelect(e.target.value);
        }}
      >
        {options.map((option) => (
          <MenuItem value={option}>{option}</MenuItem>
        ))}
      </Select>
    </FormControl>
  );
};

export default ColorSelect;
