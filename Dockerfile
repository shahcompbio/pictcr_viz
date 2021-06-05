FROM node as builder

WORKDIR /usr/src/app
COPY . .
RUN yarn install
RUN yarn build



FROM python:3.6
WORKDIR /usr/src/app

COPY render.py .
COPY requirements.txt .
RUN pip install -r requirements.txt

COPY --from=builder /usr/src/app/build ./build

CMD ["python", "render.py", "ndv.h5ad"]
