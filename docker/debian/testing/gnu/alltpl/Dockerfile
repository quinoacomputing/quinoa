FROM debian:testing
MAINTAINER Jozsef Bakosi <jbakosi@lanl.gov>

# Install prereqs
RUN apt-get -y update && apt-get install -y git cmake gcc g++ gfortran openmpi-bin libopenmpi-dev libhdf5-serial-dev libnetcdf-dev libexodusii-dev libpstreams-dev libboost-all-dev liblapack-dev liblapacke-dev

CMD ["/bin/bash"]
