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


#include "Moocho_Config.h"


#ifdef HAVE_MOOCHO_FORTRAN


#include "FortranTypes_f_open_file.hpp"
#include "FortranTypes_CppFortranStrings.hpp"
#include "Teuchos_Assert.hpp"


typedef FortranTypes::f_int	f_int;


extern "C" {


FORTRAN_FUNC_DECL_UL_(f_int,FOR_OPEN_FILE,for_open_file) ( const f_int& iunit
  , const f_int i_file[], const f_int& i_file_len, const f_int& istatus
  , const f_int& iform, const f_int& iblank, const f_int& iaccess
  , const f_int& irecl );


FORTRAN_FUNC_DECL_UL_(void,FOR_CLOSE_FILE,for_close_file) (
  const f_int& iunit, const f_int& keep );


}	// end extern "C"


void FortranTypes::f_open_file( const f_int iunit, const char file[]
  , EOpenStatus status, EOpenForm form, EOpenBlank blank
  , EOpenAccess access, f_int recl )
{
  int result;
  // Convert the file name to an integer to pass to Fortran.
  FortranTypes::f_int int_file[100], int_file_len;
  result = convert_to_f_int_string( file, int_file, &int_file_len );
  TEUCHOS_TEST_FOR_EXCEPTION(
    result, InvalidFileNameException
    ,"f_open_file(...) : Error, the "
    << -result << " Character of the file name \""
    << file << "\" is not a valid ASCII character." );

  if(
    //result = FORTRAN_FUNC_CALL_UL_(F_OPEN_FILE,f_open_file)(
    result = FORTRAN_FUNC_CALL_UL_(FOR_OPEN_FILE,for_open_file)(
      iunit, int_file, int_file_len, status, form, blank, access, recl
      )
    )
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      result < 0, InvalidFileNameException
      ,"f_open_file(...) : Error, the "
      << -result << " Character of the file name \""
      << file << "\" is not a valid ASCII character." );
    TEUCHOS_TEST_FOR_EXCEPTION(
      result > 0, OpenException
      ,"f_open_file(...) : Error, the file named \""
      << file << "\" could not be opened and OPEN(...) "
      << "returned and IOSTAT = " << result );
  }
}


void FortranTypes::f_close_file( const f_int iunit, bool keep )
{
  FORTRAN_FUNC_CALL_UL_(FOR_CLOSE_FILE,for_close_file)( iunit, keep ? 1 : 0 ); 
}


#endif // HAVE_MOOCHO_FORTRAN
