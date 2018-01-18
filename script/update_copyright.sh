################################################################################
#
# \file      script/update_copyright.sh
# \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
# \brief     Switch copyright year in all files
#
################################################################################

declare -a dirs=( "doc" "docker" "regression" "script" "src" )

for i in "${dirs[@]}"
do
  find $i -type f -exec grep -Iq . {} \; -and -print | cut -f1 -d: | grep -v update_copyright.sh | xargs sed -i 's/opyright 2012-2015, J. Bakosi, 2016-2017/opyright 2012-2015, J. Bakosi, 2016-2018/g'
  find $i -type f -exec grep -Iq . {} \; -and -print | cut -f1 -d: | grep -v update_copyright.sh | xargs sed -i 's/Copyright (c) 2016-2017/Copyright (c) 2016-2018/g'
  find $i -type f -exec grep -Iq . {} \; -and -print | cut -f1 -d: | grep -v update_copyright.sh | xargs sed -i 's/Copyright 2016-2017/Copyright 2016-2018/g'
done

find LICENSE -type f -exec grep -Iq . {} \; -and -print | cut -f1 -d: | xargs sed -i 's/Copyright (c) 2016-2017/Copyright (c) 2016-2018/g'
find LICENSE -type f -exec grep -Iq . {} \; -and -print | cut -f1 -d: | xargs sed -i 's/Copyright 2016-2017/Copyright 2016-2018/g'
