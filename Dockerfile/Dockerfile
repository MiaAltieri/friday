# USE UBUNTU
FROM ubuntu:18.04
MAINTAINER Kishwar Shafin

# Install some basic utilities
RUN apt-get update && apt-get install -y \
    curl \
    ca-certificates \
    git \
    bzip2 \
    libx11-6 \
    wget \
    autoconf \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    software-properties-common \
 && rm -rf /var/lib/apt/lists/*

# INSTALL CMAKE
RUN wget -qO- "https://cmake.org/files/v3.12/cmake-3.12.1-Linux-x86_64.tar.gz" | \
  tar --strip-components=1 -xz -C /usr/local

# Create a Python 3.6 environment
RUN add-apt-repository ppa:jonathonf/python-3.6
RUN apt-get update -y && apt-get install -y build-essential python3.6 python3.6-dev python3-pip python3.6-venv


# CLONE FRIDAY
RUN mkdir /software
WORKDIR /software

RUN git clone https://github.com/kishwarshafin/friday.git
RUN cd friday \
    && git checkout friday_v2 \
    && ./build.sh

# install pytorch
RUN pip3 install http://download.pytorch.org/whl/cpu/torch-0.4.1-cp36-cp36m-linux_x86_64.whl \
    && pip3 install torchvision

# install other python dependencines
RUN pip3 install pysam \
    && pip3 install tqdm \
    && pip3 install pandas \
    && pip3 install h5py \
    && pip3 install torchnet \
    && pip3 install graphviz