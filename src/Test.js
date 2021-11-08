import React, { useEffect, useState } from "react";

const App = () => {
  const [error, setError] = useState(null);
  const [isLoaded, setIsLoaded] = useState(false);
  const [data, setData] = useState({});

  useEffect(() => {
    fetch(
      "http://127.0.0.1:5000/load/Users/leungs1/Desktop/pictcr_viz/src/data/hacohen_viz.h5ad"
    )
      .then((res) => res.json())
      .then(
        (result) => {
          setIsLoaded(true);
          setData(result);
        },
        (error) => {
          setIsLoaded(true);
          setError(error);
        }
      );
  }, []);

  if (error) {
    return <div>Error: {error.message}</div>;
  } else if (!isLoaded) {
    return <div>Loading...</div>;
  }

  console.log(data);
  return null;
};

export default App;
