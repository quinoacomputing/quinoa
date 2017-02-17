#!/bin/bash

#-----------------------------------------------------------------------------
# Simple build script with options
#-----------------------------------------------------------------------------

INC_PATH="-I. -I../../src"
CXX="g++"
CXXFLAGS="-Wall"

CXX_SOURCES="explicit_main.cpp explicit_test_host.cpp"
CXX_SOURCES="${CXX_SOURCES} ../../src/impl/*.cpp"
CXX_SOURCES="${CXX_SOURCES} ../../src/Host/KokkosArray_Host_Impl.cpp"
CXX_SOURCES="${CXX_SOURCES} ../../src/Host/KokkosArray_Host_MemoryManager.cpp"

#-----------------------------------------------------------------------------

while [ -n "${1}" ] ; do

ARG="${1}"
shift 1

case ${ARG} in
#-------------------------------
#----------- OPTIONS -----------
CUDA | Cuda | cuda ) HAVE_CUDA=1 ;;
HWLOC | hwloc ) HAVE_HWLOC=${1} ; shift 1 ;;
OPT | opt | O3 | -O3 ) OPTFLAGS="-O3" ;;
DBG | dbg | g | -g )   OPTFLAGS="-g" ;;
#-------------------------------
#---------- COMPILERS ----------
GNU | gnu | g++ )
  CXX="g++"
  CXXFLAGS="-Wall"
  ;;
INTEL | intel | icc )
  CXX="icc"
  # -xW = use SSE and SSE2 instructions
  CXXFLAGS="-Wall -xW"
  LIB="${LIB} -lstdc++"
  ;;
#-------------------------------
*) echo 'unknown option: ' ${ARG} ; exit -1 ;;
esac
done

#-----------------------------------------------------------------------------

if [ -n "${HAVE_CUDA}" ] ;
then
  TEST_MACRO="${TEST_MACRO} -DTEST_KOKKOS_CUDA"
  NVCC_SOURCES="../../src/Cuda/*.cu explicit_test_cuda.cu"
  LIB="${LIB} -L/usr/local/cuda/lib64 libCuda.a -lcudart -lcuda -lcusparse"
  nvcc -arch=sm_20 -lib -o libCuda.a ${OPTFLAGS} ${INC_PATH} ${NVCC_SOURCES} ;
fi

#-----------------------------------------------------------------------------

if [ -n "${HAVE_HWLOC}" ] ;
then
  CXX_SOURCES="${CXX_SOURCES} ../../src/Host/KokkosArray_Host_hwloc.cpp"
  LIB="${LIB} -L${HAVE_HWLOC}/lib -lhwloc"
  INC_PATH="${INC_PATH} -I${HAVE_HWLOC}/include"
else
  CXX_SOURCES="${CXX_SOURCES} ../../src/Host/KokkosArray_Host_hwloc_unavailable.cpp"
fi

#-----------------------------------------------------------------------------
# Option for PTHREAD or WINTHREAD eventually

HAVE_PTHREAD=1

if [ -n "${HAVE_PTHREAD}" ] ;
then
  CXX_SOURCES="${CXX_SOURCES} ../../src/Host/KokkosArray_Host_pthread.cpp"
  LIB="${LIB} -lpthread"
else
  CXX_SOURCES="${CXX_SOURCES} ../../src/Host/KokkosArray_Host_nothread.cpp"
fi

#-----------------------------------------------------------------------------

echo "Building regular files as: " ${CXX} ${CXXFLAGS} ${OPTFLAGS}

${CXX} ${CXXFLAGS} ${OPTFLAGS} ${INC_PATH} ${TEST_MACRO} -o explicit_dynamics.x ${CXX_SOURCES} ${LIB}

rm -f *.o *.a ThreadPool_config.h

#-----------------------------------------------------------------------------

