FROM debian:stable

MAINTAINER Houcemeddine Othman "houcemeddine.othman@wits.ac.za"

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update

COPY dependencies.yml  /

# Download dependencies
RUN apt-get install -y  wget  \
	build-essential \ 
	git \
	autoconf \
	libgsl-dev \
	libboost-all-dev

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh 
RUN bash Miniconda3-latest-Linux-x86_64.sh -bfp /usr/local 

RUN conda env update -n base  --file dependencies.yml 
RUN conda clean -a -y

# Install freesasa
RUN wget https://freesasa.github.io/freesasa-2.0.3.tar.gz \
	 && tar xvf freesasa-2.0.3.tar.gz \
	 && cd freesasa-2.0.3  \
	 && ./configure --disable-xml --disable-json \
	 && make && make install && cd /  \
	 && rm -rf freesasa-2.0.3 freesasa-2.0.3.tar.gz 

# Install Encom
RUN git clone https://github.com/NRGlab/ENCoM.git \
	&& cd ENCoM/build/ \
	&& make -f Makefile.Unix \
	&& cd /
ENV PATH /ENCoM/bin:$PATH

# Install stride
RUN mkdir stride \
	&& cd stride \
	&& wget http://webclu.bio.wzw.tum.de/stride/stride.tar.gz \
	&& tar xvf stride.tar.gz \
	&& make \
	&& cp stride  /usr/local/bin/ \
	&& cd / && rm -rf stride 

# Install dssp 
RUN git clone https://github.com/cmbi/dssp.git \
	&& cd dssp \
	&& ./autogen.sh \ 
	&& ./configure \
	&& make && make install \
	&& cd /

# Install foldx 
RUN mkdir foldx 

COPY foldxLinux64.tar.gz /foldx 
RUN cd foldx \ 
	&& tar xvf foldxLinux64.tar.gz 

ENV PATH /foldx:$PATH

