################################################################################
# 
# \file      script/update_doc.sh
# \author    J. Bakosi
# \date      Sun 02 Aug 2015 10:06:34 PM MDT
# \copyright 2012-2015, Jozsef Bakosi.
# \brief     Update documentation and upload to github pages
# 
################################################################################

cd /home/jbakosi/code/quinoa
DOC4COMMIT=$(git rev-parse --verify HEAD)
cd build/gnu
rm -rf doc/html
git clone git@github.com:jbakosi/quinoa.git --branch gh-pages --single-branch doc/html
cd doc/html
git rm -rf .
cd -
ninja unittest_coverage
ninja regression_coverage
ninja doc
cd doc/html
touch .nojekyll
git add .
git commit -m "Automated documentation build for changeset ${DOC4COMMIT}"
git push origin gh-pages
