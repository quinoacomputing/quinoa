// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

/*! \file AnasaziConfigDefs.hpp
  \brief Anasazi header file which uses auto-configuration information to include
  necessary C++ headers
*/

#ifndef ANASAZI_CONFIGDEFS_HPP
#define ANASAZI_CONFIGDEFS_HPP

#include "Teuchos_ConfigDefs.hpp"

#ifndef __cplusplus
#  define __cplusplus
#endif

#ifndef TRILINOS_NO_CONFIG_H

/*
 * The macros PACKAGE, PACKAGE_NAME, etc, get defined for each package and need to
 * be undef'd here to avoid warnings when this file is included from another package.
 * KL 11/25/02
 */
#  ifdef PACKAGE
#    undef PACKAGE
#  endif

#  ifdef PACKAGE_NAME
#    undef PACKAGE_NAME
#  endif

#  ifdef PACKAGE_BUGREPORT
#    undef PACKAGE_BUGREPORT
#  endif

#  ifdef PACKAGE_STRING
#    undef PACKAGE_STRING
#  endif

#  ifdef PACKAGE_TARNAME
#    undef PACKAGE_TARNAME
#  endif

#  ifdef PACKAGE_VERSION
#    undef PACKAGE_VERSION
#  endif

#  ifdef VERSION
#    undef VERSION
#  endif

#  include <Anasazi_config.h>

#  ifdef HAVE_MPI
#    ifndef EPETRA_MPI
#      define EPETRA_MPI
#    endif
#  endif

#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#include <cctype>
#include <numeric>
#include <complex>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <cmath>
#include <functional>

#else /*TRILINOS_NO_CONFIG_H is defined*/

#  include <iterator>
#  include <iostream>
#  include <string>

#  if defined(SGI) || defined(SGI64) || defined(SGI32) || defined(CPLANT) || defined (TFLOP)
#    include <stdlib.h>
#    include <stdio.h>
#    include <math.h>
#  else
#    include <cstdlib>
#    include <cstdio>
#    include <cmath>
#  endif

#  include <vector>
#  include <map>
#  include <deque>
#  include <algorithm>
#  include <numeric>
#  include <functional>

#endif /*ndef TRILINOS_NO_CONFIG_H*/

/* Define some macros */
#define ANASAZI_MAX(x,y) (( (x) > (y) ) ? (x)  : (y) )     /* max function  */
#define ANASAZI_MIN(x,y) (( (x) < (y) ) ? (x)  : (y) )     /* min function  */
#define ANASAZI_SGN(x)   (( (x) < 0.0 ) ? -1.0 : 1.0 )     /* sign function */

#ifdef HAVE_TEUCHOS_COMPLEX
#  if defined(HAVE_COMPLEX)
#    define ANSZI_CPLX_CLASS std::complex
#  elif  defined(HAVE_COMPLEX_H)
#    define ANSZI_CPLX_CLASS ::complex
#  endif
#endif

#include "Anasazi_DLLExportMacro.h"

/*
 * Anasazi_Version() method 
 */
namespace Anasazi {
  ANASAZI_LIB_DLL_EXPORT std::string Anasazi_Version();
}

#endif /*ANASAZI_CONFIGDEFS_HPP*/
