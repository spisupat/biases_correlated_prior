# hash:sha256:6a967d4dcfc3b8a65c3ac877137d7a89081d9a9803e90051b4a5dbd3377c18b6
FROM registry.codeocean.com/codeocean/matlab:2020b-ubuntu20.04

ARG DEBIAN_FRONTEND=noninteractive
ARG MLM_LICENSE_FILE

RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys \
        0xAD2323F17326AE31401037733E05EBFF05441C52 \
    && apt-get update \
    && apt-get install -y --no-install-recommends software-properties-common \
    && add-apt-repository -y 'deb http://deb.codeocean.com/rstudio-server-bionic/ ubuntu main' \
    && apt-get purge -y --autoremove software-properties-common \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential=12.8ubuntu1.1 \
        libgit2-dev=0.28.4+dfsg.1-2 \
        libssl-dev=1.1.1f-1ubuntu2.13 \
        pandoc=2.5-3build2 \
        pkg-config=0.29.1-0ubuntu4 \
        r-base=3.6.3-2 \
        rstudio-server=1.2.5033 \
    && rm -rf /var/lib/apt/lists/*
