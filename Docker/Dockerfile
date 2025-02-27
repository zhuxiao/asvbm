# Use the official Ubuntu 22.04 image as the base image
FROM ubuntu:22.04

# Set up a non-interactive front end
ENV DEBIAN_FRONTEND=noninteractive

# Maintainer information
MAINTAINER xzhu <xzhu@ytu.edu.cn>

# Set working directory
WORKDIR /usr

# Create a directory to store your code and add source code and test data
RUN mkdir /usr/local/asvbm
ADD asvbm_v1.3.0.tar.gz /usr/local/asvbm/
ADD test_data.tar.gz /usr/local/asvbm/

# Install build tools, C++ compilers, dependencies, and R environments
RUN apt-get update && apt-get install -y \
    g++ \
    make \
    gnuplot \
    libhts-dev \
    r-base \
    r-base-dev \
    && apt-get clean -y

# Install the UpSetR package
RUN R -e "install.packages('UpSetR', repos='http://cran.r-project.org')"

# Set the working directory to the unzipped code directory
WORKDIR /usr/local/asvbm/asvbm_v1.3.0

# Build C++ applications
RUN make
