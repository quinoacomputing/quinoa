FROM debian:testing
MAINTAINER Jozsef Bakosi <jbakosi@lanl.gov>

# Install prereqs
RUN apt-get -y update && apt-get install -y wget
RUN wget -O - http://llvm.org/apt/llvm-snapshot.gpg.key | apt-key add -
RUN echo "deb http://llvm.org/apt/wheezy/ llvm-toolchain-wheezy main" >> /etc/apt/sources.list
RUN apt-get -y update && apt-get install -y git cmake clang-3.5 libc++-dev libc++abi-dev gfortran libhdf5-serial-dev libnetcdf-dev libexodusii-dev libpstreams-dev libboost-all-dev liblapack-dev liblapacke-dev

# Download and build OpenMPI
ENV OPENMPI=$HOME/openmpi
RUN wget http://www.open-mpi.org/software/ompi/v1.8/downloads/openmpi-1.8.6.tar.gz
RUN tar xzf openmpi-1.8.6.tar.gz
RUN cd openmpi-1.8.6 && ./configure --prefix=$OPENMPI CC=clang-3.5 CXX=clang++-3.5 FC=gfortran && make -sj4 install

# Export OpenMPI environment
ENV PATH=$OPENMPI/bin:$PATH
ENV LD_LIBRARY_PATH=$OPENMPI/lib:$LD_LIBRARY_PATH
ENV CC=$OPENMPI/bin/mpicc
ENV CXX=$OPENMPI/bin/mpic++
ENV FC=$OPENMPI/bin/mpif90

CMD ["/bin/bash"]
