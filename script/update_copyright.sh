################################################################################
#
# \file      script/update_copyright.sh
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Switch copyright year in all files
#
################################################################################

find <dirs> -type f -exec grep -Iq . {} \; -and -print | cut -f1 -d: | xargs sed -i 's/opyright 2012-2015, Jozsef Bakosi, 2016/opyright 2012-2015, J. Bakosi, 2016-2017/g'
find <dirs> -type f -exec grep -Iq . {} \; -and -print | cut -f1 -d: | xargs sed -i 's/Copyright (c) 2016/Copyright (c) 2017/g'
find <dirs> -type f -exec grep -Iq . {} \; -and -print | cut -f1 -d: | xargs sed -i 's/Copyright 2016/Copyright 2017/g'
