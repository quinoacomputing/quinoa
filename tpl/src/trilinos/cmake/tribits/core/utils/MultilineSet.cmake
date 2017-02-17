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


#
# @FUNCTION: MULTILINE_SET()
#
# Function to set a single string by concatenating a list of separate strings
#
# Usage::
#
#   MULTILINE_SET(<outputStrVar>
#     "<string0>"
#     "<string1>"
#     ...
#     )
#
# On output, the local variables ``<outputStrVar>`` is set to::
#
#   "<string0><string1>..."
#
# The purpose of this is function to make it easier to set longer strings over
# multiple lines.
#
# This function is exactly the same as `CONCAT_STRINGS()`_ and should not even
# exist :-(
#
FUNCTION(MULTILINE_SET VARAIBLE_NAME)

  SET(MULTILINE_SET_LOCAL_STR "")

  FOREACH(LINE_STR ${ARGN})
    SET(MULTILINE_SET_LOCAL_STR "${MULTILINE_SET_LOCAL_STR}${LINE_STR}")
  ENDFOREACH()

  SET(${VARAIBLE_NAME} "${MULTILINE_SET_LOCAL_STR}" PARENT_SCOPE)

ENDFUNCTION()
