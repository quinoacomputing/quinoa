#!/bin/bash -e
################################################################################
#
# \file      tools/recent_mods.sh
# \brief     Generate a list of recently modified doc pages
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \details   This script runs git and queries recently modified documentation
# pages and replaces this list in mainpage.dox, which then doxygen incorporates
# into the documentation. The final product is then renderd under section
# "Recently modified pages on the main documentation pages.
#
# A single command line argument ise required: the root of the quinoa git
# repository.
################################################################################

echo "Generating recently modified documentation pages..."

cd $1
export mods=$(paste -d "\n" <(git ls-tree -r --name-only HEAD doc/pages | grep -v directories | while read filename; do echo "$(git log -1 --format="%ai" -- $filename) $filename"; done | sort -r | grep dox$ | head -n 5 | awk '{print $4}' | xargs grep -i -- '@page\|\\page\|@mainpage\|\\mainpage' | awk '{ if (!length($3)) print "  @ref index"; else print "  @ref " $3 }') <(git ls-tree -r --name-only HEAD doc/pages | grep -v directories | while read filename; do echo "$(git log -1 --format="%ai" -- $filename) $filename"; done | sort -r | grep dox$ | head -n 5 | awk '{print $1}' | xargs -i date -d "{}" +"%b %d, %Y" | awk '{print "  @m_span{m-text m-dim} <em> " $0 " </em> @m_endspan <br>" }'))
cd -

perl -i -p0e 's/(\@par Recently modified pages).*(\@section mainpage_tools Tools)/$1\n$ENV{mods}\n\n$2/s' pages/mainpage.dox
