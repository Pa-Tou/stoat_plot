# Use a stable R base image
FROM ubuntu:24.04

# System dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libz-dev \
    liblzma-dev \
    r-base \
    r-base-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy package source
COPY . /opt/StoatPlot
WORKDIR /opt/StoatPlot

RUN R -f install.R
