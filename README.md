# PICTCR

This project is the visualization application for Pictcr.

## Available scripts


### Developing locally

To start the development server:

```
yarn start
```


## Compiling locally

To compile your changes into a local build file

```
yarn build
```

To then populate it with data (given as an h5ad file)

```
python render.py </path/to/h5ad>
```

## Compile and distribute

To share your build file, you can compile it into a docker image:

```
docker build . -t pictcr
```

Once this has been made, you can run it to generate the html file

```
docker run -v </path/to/h5ad>:/usr/src/app/ndv.h5ad pictcr
```

Then copy it into your local directory

```
docker cp <container name>:/usr/src/app/build/pictcr.html .
```

