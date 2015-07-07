################################################################################
# 
# \file      script/update_copyright.sh
# \author    J. Bakosi
# \date      Thu 29 Jan 2015 09:56:00 PM MST
# \copyright 2012-2015, Jozsef Bakosi.
# \brief     Update copyright year in all files in src and doc
# 
################################################################################

find src doc -type f -print | xargs file | grep text | cut -f1 -d: | xargs sed -i 's/opyright 2012-2014/opyright 2012-2015/g'
