#!/bin/bash -eux
# vim: filetype=sh:
################################################################################
# 
# \file      script/update_doc.sh
# \author    J. Bakosi
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Regenerate doc and test coverage and upload to github pages
# \details   This script assumes a new clone of
#   (1) github.com:quinoacomputing/quinoa.git and
#   (2) github.com:quinoacomputing/quinoacomputing.github.io.git in WORKDIR
#   (defined below), builds the third-party libraries and the code in DEBUG,
#   generates API documentation, and various coverage reports, and uploads the
#   results to (2).
# 
################################################################################

# trap and cleanup on signals INT, TERM, and EXIT
trap cleanup INT TERM EXIT

# Remove work directory
cleanup()
{
  rm -rf ${WORKDIR}
  exit
}

# Set and enter work directory
WORKDIR=/tmp/q
cd ${WORKDIR}

# Get git commit sha1 for latest code commit and for the commit which doc was
# most recently generated
cd quinoa && CODE_SHA=$(git rev-parse --verify master) && cd -
cd quinoacomputing.github.io && DOC_SHA=$(git log -1 master --pretty=%B | awk '{print $4}') && cd -

# See if doc commit sha1 equals that of code (if so, no need to regenerate doc)
if [ $CODE_SHA != $DOC_SHA ]; then

  # Query number of CPUs
  case "$OSTYPE" in
    darwin*)  CPUS=`sysctl -n hw.ncpu`;;
    linux*)   CPUS=`grep -c processor /proc/cpuinfo`;;
  esac

  # Build third-party libraries
  cd ${WORKDIR}/quinoa && mkdir tpl/build && cd tpl/build && cmake .. && make -sj$CPUS && cd -

  # Set build directory
  BUILDDIR=${WORKDIR}/quinoa/build

  # Change to work directory
  mkdir -p ${BUILDDIR} && cd ${BUILDDIR}

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

  # Generate cppcheck static analysis report and move it to ${WORKDIR}
  cmake ../src
  make cppcheck
  mv doc/cppcheck ${WORKDIR}
  rm * -rf

  # Generate documentation, combine doc with code coverage reports
  cd ${BUILDDIR}
  cmake ../src
  make -sj$CPUS doc
  cd doc/html
  mv ${WORKDIR}/unittest_coverage .
  mv ${WORKDIR}/regression_coverage .
  mv ${WORKDIR}/test_coverage .
  mv ${WORKDIR}/cppcheck .
  touch .nojekyll
  cp ../../../doc/images/* .
  cp ${WORKDIR}/quinoacomputing.github.io/README.md .

  # Push as new commit to web repo
  cd ${WORKDIR}/quinoacomputing.github.io
  git rm -rf .
  mv ${BUILDDIR}/doc/html/* ${BUILDDIR}/doc/html/.nojekyll .
  git add .
  git commit --no-gpg-sign -m "Documentation for changeset ${CODE_SHA}"
  git push

fi

cleanup
