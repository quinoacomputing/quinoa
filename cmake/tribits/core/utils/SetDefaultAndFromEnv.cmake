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



INCLUDE(SetDefault)
INCLUDE(PrintVar)

#
# @MACRO: SET_DEFAULT_AND_FROM_ENV()
#
# Set a default value for a local variable and override from an environment
# variable of the same name if it is set.
#
# Usage::
#
#   SET_DEFAULT_AND_FROM_ENV(<varName> <defaultVal>)
#
# First calls ``SET_DEFAULT(<varName> <defaultVal>)`` and then looks for an
# environment variable named ``<varName>``, and if non-empty then overrides
# the value of the local variable ``<varName>``.
#
# This macro is primarily used in CTest code to provide a way to pass in the
# value of CMake variables.  Older versions of ``ctest`` did not support the
# option ``-D <var>:<type>=<value>`` to allow variables to be set through the
# command-line like ``cmake`` always allowed.
#
MACRO(SET_DEFAULT_AND_FROM_ENV  VAR  DEFAULT_VAL)

  SET_DEFAULT(${VAR} "${DEFAULT_VAL}")

  SET(ENV_${VAR} $ENV{${VAR}})
  IF (NOT "${ENV_${VAR}}" STREQUAL "")
    PRINT_VAR(ENV_${VAR})
    SET(${VAR} ${ENV_${VAR}})
  ENDIF()

  PRINT_VAR(${VAR})

ENDMACRO()
