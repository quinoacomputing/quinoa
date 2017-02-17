# @HEADER
# ************************************************************************
#
#            TriBITS: Tribal Build, Integrate, and Test System
#                    Copyright 2013 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ************************************************************************
# @HEADER

# This script is a cmake -P script that should be called with the following
# variables defined:
#
#    -D binary_dir=${CMAKE_CURRENT_BINARY_DIR}
#    -D source_dir=${CMAKE_CURRENT_SOURCE_DIR}
#    -D ctest_type=${ctest_type}
#    -D scriptname=${scriptname}
#    -D TD_BASE_DIR=${TD_BASE_DIR}
#    -D testname=${testname}
#
# It looks recursively under ${TD_BASE_DIR}/tools/cmake-${ctest_type}
# for a ctest executable and uses it to drive a ctest -S script to run
# a dashboard for the project. It has to be run this way indirectly
# through a cmake -P script because the desired ctest executable
# location is unknown at CMake configure time of the driver project.

if(WIN32)
  set(ctest_filename "ctest.exe")
else()
  set(ctest_filename "ctest")
endif()

set(TRIBITS_TDD_USE_SYSTEM_CTEST $ENV{TRIBITS_TDD_USE_SYSTEM_CTEST})
message("TRIBITS_TDD_USE_SYSTEM_CTEST='${TRIBITS_TDD_USE_SYSTEM_CTEST}'")

if (TRIBITS_TDD_USE_SYSTEM_CTEST STREQUAL "1")

  find_program(CTEST_EXE ${ctest_filename})

elseif(NOT CTEST_EXE)

  message("Selecting '${ctest_filename}' of type '${ctest_type}'...")

  set(CTEST_EXE
    "${TD_BASE_DIR}/tools/cmake-${ctest_type}/bin/${ctest_filename}")

  if(NOT CTEST_EXE)
    message(FATAL_ERROR "error: '${ctest_type}' ctest could not be found...")
  endif()

endif()

if(NOT EXISTS "${CTEST_EXE}")
  message(FATAL_ERROR "error: CTEST_EXE='${CTEST_EXE}' does not exist...")
endif()

message("CTEST_EXE='${CTEST_EXE}'")
execute_process(COMMAND ${CTEST_EXE} --version)

message("=========== variables ===========")
message("binary_dir='${binary_dir}'")
message("source_dir='${source_dir}'")
message("ctest_type='${ctest_type}'")
message("scriptname='${scriptname}'")
message("TD_BASE_DIR='${TD_BASE_DIR}'")
message("testname='${testname}'")
message("=================================")

message("========== environment ==========")
execute_process(COMMAND ${CMAKE_COMMAND} -E environment)
message("=================================")

message("============ script =============")
message("executing ctest -S '${scriptname}' for test '${testname}'...")

execute_process(COMMAND ${CTEST_EXE}
  -S
  "${source_dir}/${scriptname}"
  -V
  --output-log
  "${binary_dir}/${testname}.log"
  RESULT_VARIABLE rv
  )

if(NOT "${rv}" STREQUAL "0")
  message("warning: calling ctest -S script failed with '${rv}'")
endif()

message("=================================")
