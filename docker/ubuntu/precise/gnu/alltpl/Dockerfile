FROM ubuntu:12.04
MAINTAINER Jozsef Bakosi <jbakosi@lanl.gov>

# Install prereqs
RUN apt-get -u update && apt-get -y install python-software-properties
RUN add-apt-repository -y ppa:ubuntu-toolchain-r/test
RUN add-apt-repository -y ppa:robotics/ppa
RUN add-apt-repository -y ppa:kalakris/cmake
RUN apt-get -y update && apt-get install -y git cmake make gcc-4.9 g++-4.9 gfortran openmpi-bin libopenmpi-dev libhdf5-serial-dev libnetcdf-dev libexodusii-dev libpstreams-dev libboost-all-dev liblapack-dev liblapacke-dev
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.9 20
RUN update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-4.9 20

CMD ["/bin/bash"]
