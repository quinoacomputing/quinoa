# @HEADER
# ************************************************************************
#
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001) Sandia Corporation
#
#
# Copyright (2001) Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of this
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
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
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in this data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
#
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
#
# ************************************************************************
# @HEADER

INCLUDE(ParseVariableArguments)
INCLUDE(GlobalSet)


#
# Perform a single equality check and update overall test statistics
#

FUNCTION(UNITTEST_COMPARE_CONST VAR_NAME CONST_VAL)

  MATH( EXPR NUMRUN ${UNITTEST_OVERALL_NUMRUN}+1 )
  GLOBAL_SET(UNITTEST_OVERALL_NUMRUN ${NUMRUN})

  MESSAGE(
    "\nCheck:\n"
    "    ${VAR_NAME} =\n"
    "    [${${VAR_NAME}}]\n"
    "  EQUALS:\n"
    "    [${CONST_VAL}]"
    )

  IF ("${${VAR_NAME}}" STREQUAL "${CONST_VAL}")
    MESSAGE("  [PASSED]\n")
    MATH( EXPR NUMPASSED ${UNITTEST_OVERALL_NUMPASSED}+1 )
    GLOBAL_SET(UNITTEST_OVERALL_NUMPASSED ${NUMPASSED})
  ELSE()
    MESSAGE("  [FAILED]\n")
    GLOBAL_SET(UNITTEST_OVERALL_PASS FALSE)
  ENDIF()

ENDFUNCTION()


#
# Perform a series regex of given strings and update overall test statistics
#

FUNCTION(UNITTEST_STRING_REGEX INPUT_STRING)

  PARSE_ARGUMENTS(
     #prefix
     PARSE
     #lists
     "REGEX_STRINGS"
     #options
     ""
     ${ARGN}
     )

  FOREACH(REGEX ${PARSE_REGEX_STRINGS})

    MATH( EXPR NUMRUN ${UNITTEST_OVERALL_NUMRUN}+1 )
    GLOBAL_SET(UNITTEST_OVERALL_NUMRUN ${NUMRUN})

    STRING(REGEX MATCH "${REGEX}" REGEX_MATCH_RESULT "${INPUT_STRING}")

    IF (REGEX_MATCH_RESULT)
      MESSAGE("  Searching for REGEX {${REGEX}}:  [PASSED]\n")
      MATH( EXPR NUMPASSED ${UNITTEST_OVERALL_NUMPASSED}+1 )
      GLOBAL_SET(UNITTEST_OVERALL_NUMPASSED ${NUMPASSED})
    ELSE()
      MESSAGE("  Searching for REGEX {${REGEX}}:  [FAILED]\n")
      GLOBAL_SET(UNITTEST_OVERALL_PASS FALSE)
    ENDIF()

  ENDFOREACH()

ENDFUNCTION()


#
# Perform a series regex of given strings and update overall test statistics
#

FUNCTION(UNITTEST_FILE_REGEX  INPUT_FILE)
  MESSAGE("\nRegexing for strings in the file '${INPUT_FILE}':\n")
  FILE(READ "${INPUT_FILE}" INPUT_FILE_STRING)
  UNITTEST_STRING_REGEX("${INPUT_FILE_STRING}" ${ARGN})
ENDFUNCTION()


#
# Print final statstics from all tests and assert final pass/fail
#

FUNCTION(UNITTEST_FINAL_RESULT  EXPECTED_NUMPASSED)
   MESSAGE("\nFinal UnitTests Result: num_run = ${UNITTEST_OVERALL_NUMRUN}\n")
  IF (UNITTEST_OVERALL_PASS)
    IF (UNITTEST_OVERALL_NUMPASSED EQUAL ${EXPECTED_NUMPASSED})  
      MESSAGE("Final UnitTests Result: PASSED (num_passed = ${UNITTEST_OVERALL_NUMPASSED})")
    ELSE()
      MESSAGE("\nError: num_passed = ${UNITTEST_OVERALL_NUMPASSED} != num_expected = ${EXPECTED_NUMPASSED}")
      MESSAGE("\nFinal UnitTests Result: FAILED\n")
      MESSAGE(SEND_ERROR "FAIL")
    ENDIF()
  ELSE()
    MESSAGE("\nFinal UnitTests Result: FAILED\n")
  ENDIF()
ENDFUNCTION()

