################################################################################
# 
# \file      script/update_copyright.sh
# \author    J. Bakosi
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Switch copyright year in all files in src and doc
# \date      Thu 05 May 2016 12:09:02 PM MDT
# 
################################################################################

find src doc script -type f -exec grep -Iq . {} \; -and -print | cut -f1 -d: | xargs sed -i 's/opyright 2012-2015, Jozsef Bakosi, 2016/opyright 2012-2015, Jozsef Bakosi, 2016-2017/g'
