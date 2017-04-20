################################################################################
#
# \file      script/update_copyright.sh
# \author    J. Bakosi
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Switch copyright year in all files in relevant directories
#
################################################################################

find src doc script cmake docker -type f -exec grep -Iq . {} \; -and -print | cut -f1 -d: | xargs sed -i 's/opyright 2012-2015, Jozsef Bakosi, 2016/opyright 2012-2015, J. Bakosi, 2016-2017/g'
