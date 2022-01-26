FROM node:12 as builder

WORKDIR /usr/src/app

COPY . .
RUN yarn install

CMD ["yarn", "run", "start"]
