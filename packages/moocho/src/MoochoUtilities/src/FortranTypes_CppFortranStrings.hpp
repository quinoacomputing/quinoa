// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef CPP_FORTRAN_STRINGS_H
#define CPP_FORTRAN_STRINGS_H

#include "Moocho_ConfigDefs.hpp"
#include "Teuchos_F77_wrappers.h"

namespace FortranTypes {

/** @name Pass strings to and from Fortran and C++.
  *
  * This is a portable method to pass strings from C++ to Fortran.
  * The basic idea is to convert a C++ (or C) null terminated #char[]#
  * string to a Fortran compatable integer array #f_int[]# with an
  * integer length #f_int#.  Then in Fortran this integer array is
  * converted into a Fortran CHARACTER* string.  For this to work 
  * the ASCII numeric representation in the target C++ and Fortran
  * compilers should be the same.  If they are not then these
  * encapsulation routines may be specialized to make the mapping work
  * , either  on the C++ side or the Fortran side.
  *
  * The Mapping routines are:
  *
  * C++\\
  * ----------\\
  \verbatim
  int convert_to_f_int_string( const char string[], f_int** i_string
    , f_int* string_len );
  int convert_from_f_int_string( const f_int i_string[], f_int string_len
    , char string[] );
  \endverbatim
  *
  * Fortran\\
  * ----------\\
  \verbatim
        INTEGER FUNCTION CONVERT_FROM_C_INT_STR( I_STRING, STRING_LEN, STRING )
        INTEGER I_STRING(*), STRING_LEN
      CHARACTER*(*) STRING(*)
      ...
      END

        INTEGER FUNCTION CONVERT_TO_C_INT_STR( STRING, I_STRING, STRING_LEN )
      CHARACTER*(*) STRING(*)
        INTEGER I_STRING(*), STRING_LEN
      ...
      END
  \endverbatim
  *
  *	When using this care must be taken that enough memory is allocated to
  * hold the target strings and integer arrays before calling these converstion
  * precedures in each language.
  *
  * For example, here is how a you can pass a string from C++ to Fortran in
  * order to open up a Fortran file:
  *
  * C++ \\
  * ----------\\
  \verbatim
  char file_name = "myfile.opt";
  f_int i_file_name[20];
  f_int file_name_len;
  convert_to_f_int_string( file_name, &i_file_name, &file_name_len );
  f_int f_file_num;
  F_OPEN_FILE( i_file_name, file_name_len, &f_file_num );
  // Now we can call Fortran routines to read and write to the
  // Fortran file with the unit number f_file_num  
  \endverbatim
  *
  * Fortran \\
  * ----------\\
  \verbatim
        SUBROUTINE F_OPEN_FILE( I_FILE_NAME, FILE_NAME_LEN, I_FILE_NUM )
      *** Formal parameters
      INTEGER I_FILE_NAME(*), FILE_NAME_LEN, I_FILE_NUM
      *** Local variables
      CHARACTER*20 FILE_NAME
      CALL CONVERT_FROM_C_INT_STR( I_FILE_NAME, FILE_NAME_LEN, FILE_NAME )
      OPEN( I_FILE_NUM, FILE = FILE_NAME )
      END 
  \endverbatim
  */
//@{

/** \brief Convert from a C++ (or C) null terminated string to an integer array
  * suitable for passing to Fortran and later conversion to a CHARACTER
  * string.
  *
  * @param	string	[I] Null terminated C++ (C) string.
  *	@param	i_string	[O]	Integer array containing the string
  *	@param	string_len [O] Length of the string.
  *
  * @return Returns 0 if successful.  Otherwise it returns the indice (1,...)
  *			of the ith character that is not ASCII convertable.
  */
int convert_to_f_int_string( const char string[], f_int i_string[]
  , f_int* string_len );

/** \brief Convert from a C++ (or C) null terminated string to an integer array
  * suitable for passing to Fortran and later conversion to a CHARACTER
  * string.
  *
  *	@param	i_string	[I]	Integer array containing the string
  *	@param	string_len [I] Length of the string.
  * @param	string	[O] Null terminated C++ (C) string.
  *
  * @return Returns 0 if successful.  Otherwise it returns the indice (1,...)
  *			of the ith character that is not ASCII convertable.
  */
int convert_from_f_int_string( const f_int i_string[], f_int string_len
  , char string[] );

//@}


}	// end namespace FortranTypes

#endif // CPP_FORTRAN_STRINGS_H
