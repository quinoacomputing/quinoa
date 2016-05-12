#!/bin/bash -eux
# vim: filetype=sh:
################################################################################
# 
# \file      script/update_doc.sh
# \author    J. Bakosi
# \date      Wed 11 May 2016 11:59:39 AM MDT
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Regenerate doc and test coverage and upload to github pages
# \details   This script clones the github repository, builds the third-party
#   libraries and the code in DEBUG, generates the documentation, and uploads to
#   github.
# 
################################################################################

# trap and cleanup on signals INT, TERM, and EXIT
trap cleanup INT TERM EXIT

# Set work directory
WORKDIR=/tmp/q
rm -rf ${WORKDIR} && mkdir ${WORKDIR}

# Remove work directory
cleanup()
{
  rm -rf ${WORKDIR}
  exit
}

# Clone from github
cd ${WORKDIR}
git clone git@github.com:quinoacomputing/quinoa.git
cd ${WORKDIR}/quinoa

# Get git commit sha1 for code and doc branches
CODE_SHA=$(git rev-parse --verify master)
DOC_SHA=$(git rev-parse --verify origin/gh-pages)

# See if doc commit sha1 equals that of code (if so, no need to regenerate doc)
if [ $CODE_SHA != $DOC_SHA ]; then

  # Query number of CPUs
  case "$OSTYPE" in
    darwin*)  CPUS=`sysctl -n hw.ncpu`;;
    linux*)   CPUS=`cat /proc/cpuinfo | grep MHz | wc -l`;;
  esac

  # Build third-party libraries
  mkdir tpl/build && cd tpl/build && cmake .. && make -sj$CPUS

  # Set build directory
  BUILDDIR=${WORKDIR}/quinoa/gh-pages

  # Change to work directory
  mkdir ${BUILDDIR}
  cd ${BUILDDIR}

  # Generate unit test coverage report, move it to ${WORKDIR}, and clean
  cmake -G Ninja ../src
  ninja unittest_coverage
  mv doc/html/unittest_coverage ${WORKDIR}
  rm * .ninja_* -rf

  # Generate regression test coverage report, move it to ${WORKDIR}, and clean
  cmake -G Ninja ../src
  ninja regression_coverage
  mv doc/html/regression_coverage ${WORKDIR}
  rm * .ninja_* -rf

  # Generate full test coverage report, move it to ${WORKDIR}, and clean
  cmake -G Ninja ../src
  ninja test_coverage
  mv doc/html/test_coverage ${WORKDIR}
  rm * .ninja_* -rf

  # Start out in empty build-dir with a fresh clone of the gh-pages branch
  git clone git@github.com:quinoacomputing/quinoa.git --branch gh-pages --single-branch doc/html
  cd doc/html
  git rm -rf .
  cd -

  # Generate documentation, move code coverage reports in place, and push
  cmake -G Ninja ../src
  ninja doc
  cd doc/html
  mv ${WORKDIR}/unittest_coverage .
  mv ${WORKDIR}/regression_coverage .
  mv ${WORKDIR}/test_coverage .
  touch .nojekyll
  cp ../../../doc/images/quinoa.svg .
  git add .
  git commit --no-gpg-sign -m "Documentation for changeset ${CODE_SHA}"
  git push origin gh-pages

fi

cleanup
