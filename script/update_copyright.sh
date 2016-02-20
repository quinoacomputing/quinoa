################################################################################
# 
# \file      script/update_copyright.sh
# \author    J. Bakosi
# \date      Sat 20 Feb 2016 12:13:27 PM MST
# \copyright 2012-2015, Jozsef Bakosi.
# \brief     Update copyright year in all files in src and doc
# 
################################################################################

find src doc script -type f -print | xargs file | grep text | cut -f1 -d: | xargs sed -i 's/opyright 2012-2015/opyright 2012-2016/g'
