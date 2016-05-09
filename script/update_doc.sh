################################################################################
# 
# \file      script/update_doc.sh
# \author    J. Bakosi
# \date      Sun 02 Aug 2015 10:06:34 PM MDT
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Update Sun 08 May 2016 06:26:26 PM MDT
# \details   This script assumes that
#   - a clone already exists and the TPLs are already built
#   - the environment is setup for running cmake
# 
################################################################################

# Set quinoa directory
QUINOA=/home/jbakosi/code/quinoa
# Set work directory
WORKDIR=${QUINOA}/gh-pages

# Get commit we are generating documentation for
cd ${QUINOA}
DOC4COMMIT=$(git rev-parse --verify HEAD)

# Change to work directory
mkdir ${WORKDIR}
cd ${WORKDIR}

# Generate unit test coverage report, move it to /tmp, and clean
cmake -G Ninja ../src
ninja unittest_coverage
mv doc/html/unittest_coverage /tmp
rm * .ninja_* -rf

# Generate regression test coverage report, move it to /tmp, and clean
cmake -G Ninja ../src
ninja regression_coverage
mv doc/html/regression_coverage /tmp
rm * .ninja_* -rf

# Generate full test coverage report, move it to /tmp, and clean
cmake -G Ninja ../src
ninja test_coverage
mv doc/html/test_coverage /tmp
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
mv /tmp/unittest_coverage .
mv /tmp/regression_coverage .
mv /tmp/test_coverage .
touch .nojekyll
cp ../../../doc/images/quinoa.svg .
git add .
git commit -m "Documentation for changeset ${DOC4COMMIT}"
git push origin gh-pages

# Remove work directory
cd ${QUINOA}
rm -rf ${WORKDIR}
