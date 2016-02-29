# vim: filetype=dockerfile:

FROM debian:testing
MAINTAINER Jozsef Bakosi <jbakosi@lanl.gov>

# From behind LANL firewall
ENV http_proxy http://proxyout.lanl.gov:8080/
ENV https_proxy https://proxyout.lanl.gov:8080/

# Install system-wide prerequisites, including gcc
RUN echo "deb http://ftp.us.debian.org/debian squeeze main non-free contrib \n deb http://ftp.us.debian.org/debian testing main non-free contrib \n deb http://security.debian.org testing/updates main contrib non-free" > /etc/apt/sources.list
RUN apt-get -y update && apt-get install -y m4 zlib1g-dev procps cpio git cmake gcc g++ gfortran environment-modules libpugixml-dev libpstreams-dev libboost-all-dev liblapack-dev liblapacke-dev ninja-build gmsh

# Switch to bash
RUN rm /bin/sh && ln -s /bin/bash /bin/sh

# Populate and setup environment modules
RUN mkdir -p /usr/share/modules/modulefiles/gnu /usr/share/modules/modulefiles/openmpi/1.10.2/gnu
COPY modules/gnu/system /usr/share/modules/modulefiles/gnu/
COPY modules/openmpi/1.10.2/gnu/system /usr/share/modules/modulefiles/openmpi/1.10.2/gnu/

# Install OpenMPI with gnu
ADD openmpi-1.10.2.tar.bz2 /install
RUN source /usr/share/modules/init/bash && module load gnu/system && cd /install/openmpi-1.10.2 && ./configure --prefix=/opt/openmpi/1.10.2/gnu/system && make -sj36 install
RUN rm -rf /install/openmpi-1.10.2

# Create a non-root user 'quinoa'
RUN mkdir /home/quinoa && groupadd -r quinoa -g 433 && useradd -u 431 -r -g quinoa -d /home/quinoa -s /sbin/nologin -c "Quinoa user" quinoa && chown -R quinoa:quinoa /home/quinoa
# Run init script to initialize shell
COPY init-gnu.sh /
ENTRYPOINT ["/init-gnu.sh"]
# Switch default user to 'quinoa' and change default work directory
USER quinoa
WORKDIR /home/quinoa
# Set bash as default interactive shell
CMD ["/bin/bash"]

# Clone quinoa
RUN git clone https://github.com/jbakosi/quinoa.git

# Build TPLs
RUN cd quinoa && mkdir -p tpl/build/gnu && cd tpl/build/gnu && cmake -D CMAKE_BUILD_TYPE=Release ../.. && make -sj36
# Build quinoa
RUN cd quinoa && mkdir -p build/gnu && cd build/gnu && cmake -G Ninja -D CMAKE_BUILD_TYPE=Release ../../src && ninja
