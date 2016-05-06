################################################################################
#
# \file      cmake/DetectOS.cmake
# \author    J. Bakosi
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Detect operating system
# \date      Fri 06 May 2016 06:41:22 AM MDT
#
################################################################################

function(detect_os)
  if(APPLE)
    set(os "mac")
  elseif(WIN32)
    set(os "windows")
  else()
    execute_process(COMMAND cat /etc/os-release
                    COMMAND grep ^ID
                    COMMAND awk -F= "{ print $2 }"
                    COMMAND tr "\n" " "
                    COMMAND sed "s/ //"
                    OUTPUT_VARIABLE os)
  endif()
  set(HOST_OS ${os} CACHE STRING "Name of the operating system (or Linux distribution)." FORCE)
  mark_as_advanced(HOST_OS)
endfunction()
