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

INCLUDE(ConcatStrings)

#
# @FUNCTION: APPEND_STRING_VAR_WITH_SEP()
#
# Append strings to a given string variable, joining them using a separator
# string.
#
# Usage::
#
#   APPEND_STRING_VAR_WITH_SEP(<stringVar> "<sepStr>" "<str0>" "<str1>" ...)
#
# Each of the strings ``<stri>`` are appended to ``<stringVar>`` using the
# separation string ``<sepStr>``.
#
FUNCTION(APPEND_STRING_VAR_WITH_SEP  STRING_VAR  SEP_STR)
  #MESSAGE("APPEND_STRING_VAR: '${STRING_VAR}' '${SEP_STR}' ${ARGN}")
  #PRINT_VAR(STRING_VAR)
  #PRINT_VAR(${STRING_VAR})
  IF (${STRING_VAR})
    CONCAT_STRINGS( TMP_STRING "${${STRING_VAR}}${SEP_STR}" ${ARGN} )
  ELSE()
    CONCAT_STRINGS( TMP_STRING ${ARGN} )
  ENDIF()
  #PRINT_VAR( TMP_STRING )
  SET(${STRING_VAR} "${TMP_STRING}" PARENT_SCOPE)
  #SET(${STRING_VAR} "${${STRING_VAR}}${LINE}" PARENT_SCOPE)
  #PRINT_VAR(STRING_VAR)
ENDFUNCTION()
