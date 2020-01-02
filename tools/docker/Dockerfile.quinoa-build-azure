################################################################################
# vim: filetype=dockerfile:
#
# \file      tools/docker/Dockerfile.quinoa-build-azure
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \brief     Dockerfile for continuous integration with Microsoft Azure
# \details   Dockerfile for continuous integration with Microsoft Azure
#
################################################################################

FROM quinoacomputing/buildenv:debian

ARG COMPILER
ARG BUILD
ARG SHARED_LIBS
ARG DOC
ARG SMP

ENV PATH=/opt/openmpi/${COMPILER}/bin:$PATH
ENV RUNNER_ARGS="--bind-to none -oversubscribe"
ENV OMPI_MCA_btl_base_warn_component_unused="0"

USER quinoa
COPY . /home/quinoa/
RUN mkdir build
WORKDIR build
RUN cmake \
    -GNinja \
    -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_CXX_COMPILER=mpicxx \
    -DTPL_DIR="/home/quinoa/external/${COMPILER}${SMP:+-smp}${SHARED_LIBS:+-static}" \
    -DCMAKE_CXX_FLAGS="-Werror -ftemplate-depth=1024" \
    -DCMAKE_BUILD_TYPE="${BUILD}" \
    -DCMAKE_CXX_FLAGS_DEBUG=-g0 \
    ${SHARED_LIBS:+-DBUILD_SHARED_LIBS=$SHARED_LIBS} \
    -DRUNNER=mpirun \
    -DRUNNER_NCPUS_ARG=-n \
    -DRUNNER_ARGS="${RUNNER_ARGS}" \
    ${SMP:+-DPOSTFIX_RUNNER_ARGS="+setcpuaffinity"} \
    ../src

RUN if [ "${DOC}" ]; then \
      ninja doc; \
    else \
      ninja && if [ "${SMP}" ]; then mpirun -n 1 ${RUNNER_ARGS} Main/unittest -v -q +ppn 2; else mpirun -n 2 ${RUNNER_ARGS} Main/unittest -v -q; fi && ctest -j2 --output-on-failure -LE extreme; \
    fi

USER root
RUN [ ${DOC} ] || ninja install
USER quinoa
RUN [ ${DOC} ] || (export PATH=/usr/local/bin:$PATH && unittest -h && inciter && rngtest && meshconv && walker)
WORKDIR /home/quinoa
