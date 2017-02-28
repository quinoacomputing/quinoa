FROM quinoacomputing/buildenv:gentoo

ARG COVERAGE
ARG CC
ARG CXX
ENV OMPI_CXX=${CXX}
ARG CXXFLAGS

#for coverage
ENV CI true
ENV TRAVIS true
ARG TRAVIS_BRANCH
ARG TRAVIS_JOB_NUMBER
ARG TRAVIS_PULL_REQUEST
ARG TRAVIS_JOB_ID
ARG TRAVIS_TAG
ARG TRAVIS_REPO_SLUG
ARG TRAVIS_COMMIT

RUN rm -rf /home/quinoa/.ccache
COPY quinoa /home/quinoa/quinoa
COPY ccache/ /home/quinoa/.ccache
USER root
RUN chown -R quinoa:quinoa /home/quinoa/quinoa /home/quinoa/.ccache
USER quinoa

WORKDIR /home/quinoa/quinoa
RUN mkdir build

WORKDIR build
RUN ccache -z
RUN cmake \
    ${COVERAGE:+-DLCOV=/bin/true} \
    -DGENHTML=/bin/true \
    -DRUNNER_ARGS="" \
    ../src
RUN make VERBOSE=1 -j2
RUN ccache -s
RUN ctest -V -R walker:Dirichlet_pe1
RUN make install DESTDIR=${PWD}
RUN cd .. && if [ ${COVERAGE} ]; then \
  if [ ${CC} = clang ]; then \
    $HOME/.local/bin/codecov --gcov-exec "llvm-cov gcov"; \
  else \
    $HOME/.local/bin/codecov; \
  fi; \
fi && cd -
USER root
RUN make install
USER quinoa
WORKDIR /home/quinoa
