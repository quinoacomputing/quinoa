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

#ifndef F_OPEN_FILE_H
#define F_OPEN_FILE_H

#include "Moocho_ConfigDefs.hpp"
#include "Teuchos_F77_wrappers.h"

namespace FortranTypes {

/** @name Open a Fortran file.
  */
//@{

/** \brief . */
enum EOpenStatus { OPEN_OLD = 0, OPEN_NEW = 1, OPEN_SCRATCH = 2
  , OPEN_UNKNOWN = 3 };
/** \brief . */
enum EOpenForm { OPEN_FORMATTED = 0, OPEN_UNFORMATTED = 1 };
/** \brief . */
enum EOpenBlank { OPEN_NULL = 0, OPEN_ZERO = 1 };
/** \brief . */
enum EOpenAccess { OPEN_SEQUENTIAL = 0, OPEN_DIRECT = 1 };

/** Open a Fortran file given its name and unit number.
  *
  * If successful #iunit# is returned for the opened file.
  *
  * The standard options to Fortran OPEN(...) are included.  The
  * only mandatory option to set is for the file name in FILE.
  *
  * The length of the file name must be <= 100 characters long.
  *
  * Note that all of the options are optional but if
  * you set access == OPEN_DIRECT you must set recl to some
  * value greater than zero.  See Metcalf, 1990, p. 127.
  *
  * If #file# is not a valid ASCII string then the exception
  * InvalidFileNameException will be thrown
  *
  * If the file could not be opened for some reason then
  * the exception OpenException will be thrown.
  */
void f_open_file( const f_int iunit, const char file[]
  , EOpenStatus status = OPEN_UNKNOWN, EOpenForm form = OPEN_FORMATTED
  , EOpenBlank blank = OPEN_NULL, EOpenAccess access = OPEN_SEQUENTIAL
  , f_int recl = -1 );

/** Close a Fortran file given its unit number.
  *
  * If successful #iunit# is returned for the opened file.
  *
  * The standard options to Fortran OPEN(...) are included.  The
  * only mandatory option to set is for the file name in FILE.
  *
  * The length of the file name must be <= 100 characters long.
  *
  * Note that all of the options are optional but if
  * you set access == OPEN_DIRECT you must set recl to some
  * value greater than zero.  See Metcalf, 1990, p. 127.
  *
  * If #file# is not a valid ASCII string then the exception
  * InvalidFileNameException will be thrown
  *
  * If the file could not be opened for some reason then
  * the exception OpenException will be thrown.
  */
void f_close_file( const f_int iunit, bool keep = true );

//@}

/// Thrown if the file name is not a valid ASCII string
class InvalidFileNameException : public std::logic_error
{public: InvalidFileNameException(const std::string& what_arg) : std::logic_error(what_arg) {}};

/// Thrown if the open operation fails
class OpenException : public std::logic_error
{public: OpenException(const std::string& what_arg) : std::logic_error(what_arg) {}};

}	// end namespace FortranTypes 

#endif // F_OPEN_FILE_H
