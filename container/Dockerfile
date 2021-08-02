FROM debian:stable

MAINTAINER Houcemeddine Othman "houcemeddine.othman@wits.ac.za"

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update

# Download and install miniconda
RUN apt-get install -y  wget  \
	build-essential \ 
	git 

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash Miniconda3-latest-Linux-x86_64.sh -bfp /usr/local && rm -rf Miniconda3-latest-Linux-x86_64.sh

# Insatall dependencies
COPY dependencies.yml  /

RUN conda env update -n base  --file dependencies.yml && conda clean -a -y

RUN wget https://freesasa.github.io/freesasa-2.0.3.tar.gz \
	 && tar xvf freesasa-2.0.3.tar.gz \
	 && cd freesasa-2.0.3 && ./configure --disable-xml --disable-json && make && make install && cd / && rm -rf freesasa-2.0.3 freesasa-2.0.3.tar.gz

	 