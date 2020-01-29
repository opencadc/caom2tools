FROM python:3.8-alpine

RUN apk --no-cache add make gcc musl-dev libffi-dev libxslt-dev openssl-dev

WORKDIR /usr/src/app

COPY . /usr/src/app/

RUN pip install .

RUN apk --no-cache del make gcc
