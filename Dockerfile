# Use a stable R base image
FROM rocker/r-ver:4.3.2

# System dependencies (minimal for ggplot2 stack)
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /usr/local/src

# Copy package source
COPY . /usr/local/src/StoatPlot

# Install R package dependencies
RUN R -e "install.packages(c( \
    'ggplot2', \
    'dplyr', \
    'tidyr', \
    'tidyselect', \
    'magrittr', \
    'rlang' \
), repos='https://cloud.r-project.org')"

# Install StoatPlot
RUN R -e "install.packages('StoatPlot', repos=NULL, type='source')"

# Default command
CMD ["R"]
