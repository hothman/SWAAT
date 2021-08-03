FROM debian:stable

MAINTAINER Houcemeddine Othman "houcemeddine.othman@wits.ac.za"
LABEL workflow.name="SWAAT"
LABEL github.repo="https://github.com/hothman/SWAAT"
LABEL container.size="3.73GB"

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update

# Download dependencies
RUN apt-get install -y  wget  \
	build-essential \ 
	git \
	autoconf \
	libgsl-dev \
	libboost-all-dev

# Install conda 
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh 
RUN bash Miniconda3-latest-Linux-x86_64.sh -bfp /usr/local 

# Install conda dependencies
COPY dependencies.yml  /
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
# source tar file must be renamed to foldxLinux64.tar.gz
RUN mkdir foldx 
COPY foldxLinux64.tar.gz /foldx 
RUN cd foldx \ 
	&& tar xvf foldxLinux64.tar.gz 

# add foldx and encom to PATH 
ENV PATH /foldx:/ENCoM/bin:$PATH

# pull SWAAT
RUN git clone https://github.com/hothman/SWAAT.git

# clean the container
RUN apt-get -qq -y remove wget git build-essential \
	&& apt-get -qq -y autoremove \
	&& apt-get autoclean \ 
	&& rm -rf /var/lib/apt/lists/* /var/log/dpkg.log \ 
	&& rm /foldx/foldxLinux64.tar.gz \
	&& rm /Miniconda3-latest-Linux-x86_64.sh \
	&& conda clean --all --yes