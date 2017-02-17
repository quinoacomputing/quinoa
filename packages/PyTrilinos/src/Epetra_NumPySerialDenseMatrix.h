// @HEADER
// ***********************************************************************
//
//              PyTrilinos: Python Interface to Trilinos
//                 Copyright (2005) Sandia Corporation
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

#ifndef EPETRA_NUMPYSERIALDENSEMATRIX_H
#define EPETRA_NUMPYSERIALDENSEMATRIX_H

#define NO_IMPORT_ARRAY
#include "numpy_include.h"

#include "PyTrilinos_PythonException.h"
#ifdef HAVE_INTTYPES_H
#undef HAVE_INTTYPES_H
#endif
#ifdef HAVE_STDINT_H
#undef HAVE_STDINT_H
#endif
#include "Epetra_SerialDenseMatrix.h"

namespace PyTrilinos
{

class Epetra_NumPySerialDenseMatrix : public Epetra_SerialDenseMatrix
{
public:

  // Constructors
  Epetra_NumPySerialDenseMatrix(int set_object_label=1);
  Epetra_NumPySerialDenseMatrix(int numRows, int numCols, int set_object_label=1);
  Epetra_NumPySerialDenseMatrix(PyObject * pyObject, int set_object_label=1);
  Epetra_NumPySerialDenseMatrix(const Epetra_SerialDenseMatrix & src);

  // Destructor
  ~Epetra_NumPySerialDenseMatrix();

  // Overridden Epetra_SerialDenseMatrix methods.  These are overriden
  // for one of two reasons: (1) to provide a more python-like
  // signature, or (2) to maintain synchronization between the
  // Epetra_SerialDenseMatrix and the numpy array.
  double     operator() (int rowIndex, int colIndex);
  int        Shape(  int numRows, int numCols);
  int        Reshape(int numRows, int numCols);
  PyObject * A();

  // Static cleanup function, to be called when python exceptions are
  // encountered
  static void cleanup();

private:

  // These private static pointers are for use with constructors only.
  // Outside of a constructor, they should be NULL.  they should only
  // be set by the static helper functions below, and when a
  // constructor sets these pointers, it should reset them to NULL
  // when done.
  static PyArrayObject * tmp_array;
  static bool            tmp_bool;

  // Static helper functions.  These are intended to be called from
  // the constructors, specifically to compute arguments in the
  // Epetra_SerialDenseMatrix constructors called in the constructor
  // initialization lists.  They all assume that if tmp_array is
  // already set, it has been set by the same PyObject.
  static double * getArray(  PyObject *, int);
  static int      getNumRows(PyObject *, int);
  static int      getNumCols(PyObject *, int);
  static bool     getBool(   PyObject *, int);

  // Private method.  This method is typically called after an
  // Epetra_SerialDenseMatrix constructor has been called to
  // synchronize the internal PyArrayObject.
  void setArray(bool copy=false);

  // Private data
  PyArrayObject * array;

};

}  // Namespace PyTrilinos

#endif
