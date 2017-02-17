// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//              PyTrilinos: Python Interface to Trilinos
//                 Copyright (2010) Sandia Corporation
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

// Teuchos_Array.i is a SWIG interface file that provides SWIG
// directives to handle Teuchos Array types.  These classes are not
// wrapped, but instead typemaps are defined so that the python user
// can use NumPy arrays instead.  Currently, the following classes are
// handled:
//
//     Teuchos::ArrayView< T >
%{
#include "Teuchos_ArrayView.hpp"
using Teuchos::ArrayView;
%}
#define REFCOUNTPTR_INLINE
%import  "Teuchos_ArrayView.hpp"

////////////////////////////////////////////////////////////////////////
// The philosophy is that wherever Teuchos Array classes are used in
// C++, NumPy arrays will be used in python.  Thus we need the NumPy
// SWIG directives.
%include "numpy.i"

////////////////////////////////////////////////////////////////////////
// Define a macro that takes a C++ data type (TYPE) and a
// corresponding NumPy typecode (TYPECODE) and define all of the
// typemaps needed to handle that TYPE array.
%define %teuchos_array_typemaps(TYPE, TYPECODE)

// If an ArrayView argument has a const TYPE, then we know that the
// argument is input only.  Therefore we allow any type of sequence to
// be converted to a PyArrayObject and then extract the resulting data
// pointer to construct the ArrayView.  If the conversion creates a
// new PyArrayObject, then we have to be sure to decrement its
// reference count once the ArrayView has been used.
%typemap(in) Teuchos::ArrayView< const TYPE > (int is_new = 0,
					       PyArrayObject * npArray = NULL)
{
  npArray = obj_to_array_contiguous_allow_conversion($input, TYPECODE, &is_new);
  if (!npArray) SWIG_fail;
  $1 = Teuchos::arrayView( (TYPE*) array_data(npArray), array_size(npArray, 0));
}
%typemap(freearg) Teuchos::ArrayView< const TYPE >
{
  if (is_new$argnum) Py_DECREF(npArray$argnum);
}

// If an ArrayView argument has a non-const TYPE, then the default
// behavior is to assume that the array is input/output.  Therefore
// the input python argument must be a NumPy array.
%typemap(in) Teuchos::ArrayView< TYPE >
{
  PyArrayObject * npArray = obj_to_array_no_conversion($input, TYPECODE);
  if (!npArray) SWIG_fail;
  $1 = Teuchos::arrayView( (TYPE*) array_data(npArray), array_size(npArray, 0));
}

// If an ArrayView is output, with either a const or non-const TYPE,
// convert the underlying data to a NumPy array of correct type.
%typemap(out) Teuchos::ArrayView< TYPE >
{
  npy_intp dims[1] = { $1.size() };
  $result = PyArray_SimpleNewFromData(1, dims, TYPECODE, (void*) $1.getRawPtr());
  if (!$result) SWIG_fail;
}
%typemap(out) Teuchos::ArrayView< const TYPE >
{
  npy_intp dims[1] = { $1.size() };
  $result = PyArray_SimpleNewFromData(1, dims, TYPECODE, (void*) $1.getRawPtr());
  if (!$result) SWIG_fail;
}

%enddef

////////////////////////////////////////////////////////////////////////
// Call the %teuchos_array_typemaps() macro for specific data types
// that are supported by NumPy
%teuchos_array_typemaps(signed char       , NPY_BYTE     )
%teuchos_array_typemaps(unsigned char     , NPY_UBYTE    )
%teuchos_array_typemaps(short             , NPY_SHORT    )
%teuchos_array_typemaps(unsigned short    , NPY_USHORT   )
%teuchos_array_typemaps(int               , NPY_INT      )
%teuchos_array_typemaps(unsigned int      , NPY_UINT     )
%teuchos_array_typemaps(long              , NPY_LONG     )
%teuchos_array_typemaps(unsigned long     , NPY_ULONG    )
%teuchos_array_typemaps(long long         , NPY_LONGLONG )
%teuchos_array_typemaps(unsigned long long, NPY_ULONGLONG)
%teuchos_array_typemaps(float             , NPY_FLOAT    )
%teuchos_array_typemaps(double            , NPY_DOUBLE   )
