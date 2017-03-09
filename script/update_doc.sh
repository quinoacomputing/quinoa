#!/bin/bash -eux
# vim: filetype=sh:
################################################################################
# 
# \file      script/update_doc.sh
# \author    J. Bakosi
# \date      Thu 09 Mar 2017 07:05:49 AM MST
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Regenerate doc and test coverage and upload to github pages
# \details   This script assumes a new clone in /tmp/q, builds the third-party
#   libraries and the code in DEBUG, generates the documentation, and uploads to
#   github.
# 
################################################################################

# trap and cleanup on signals INT, TERM, and EXIT
trap cleanup INT TERM EXIT

# Set work directory
WORKDIR=/tmp/q

# Remove work directory
cleanup()
{
  rm -rf ${WORKDIR}
  exit
}

# Clone from github
cd ${WORKDIR}/quinoa

# Get git commit sha1 for latest code commit and for the commit which doc was
# most recently generated
CODE_SHA=$(git rev-parse --verify master)
DOC_SHA=$(git log -1 origin/gh-pages --pretty=%B | awk '{print $4}')

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
  cmake -DCOVERAGE=on ../src
  make -sj$CPUS unittest_coverage
  mv doc/html/unittest_coverage ${WORKDIR}
  rm * -rf

  # Generate regression test coverage report, move it to ${WORKDIR}, and clean
  cmake -DCOVERAGE=on ../src
  make -sj$CPUS regression_coverage
  mv doc/html/regression_coverage ${WORKDIR}
  rm * -rf

  # Generate full test coverage report, move it to ${WORKDIR}, and clean
  cmake -DCOVERAGE=on ../src
  make -sj$CPUS test_coverage
  mv doc/html/test_coverage ${WORKDIR}
  rm * -rf

  # Start out in empty build-dir with a fresh clone of the gh-pages branch
  git clone git@github.com:quinoacomputing/quinoa.git --branch gh-pages --single-branch doc/html
  cd doc/html
  git rm -rf .
  cd -

  # Generate documentation, move code coverage reports in place, and push
  cmake ../src
  make -sj$CPUS doc
  cd doc/html
  mv ${WORKDIR}/unittest_coverage .
  mv ${WORKDIR}/regression_coverage .
  mv ${WORKDIR}/test_coverage .
  touch .nojekyll
  cp ../../../doc/images/* .
  git add .
  git commit --no-gpg-sign -m "Documentation for changeset ${CODE_SHA}"
  git push origin gh-pages

fi

cleanup
