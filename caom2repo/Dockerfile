ARG OPENCADC_PYTHON_VERSION=3.12
FROM opencadc/matplotlib:${OPENCADC_PYTHON_VERSION}-slim as builder

RUN apt-get update --no-install-recommends && \
    apt-get dist-upgrade -y && \
    apt-get install -y build-essential git libcfitsio-bin&& \
    rm -rf /var/lib/apt/lists/ /tmp/* /var/tmp/*

WORKDIR /usr/src/app

ARG OPENCADC_BRANCH=main
ARG OPENCADC_REPO=opencadc

RUN git clone https://github.com/${OPENCADC_REPO}/caom2tools.git && \
    cd caom2tools && \
    git checkout ${OPENCADC_BRANCH} && \
    pip install ./caom2repo && \
    cd ..

#COPY ./docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh
#ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]
