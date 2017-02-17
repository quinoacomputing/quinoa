/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef SUPPORT_DH
#define SUPPORT_DH

#include "euclid_common.h"

#ifdef __cplusplus
extern "C"
{
#endif

  extern void shellSort_int (const int n, int *x);
  extern void shellSort_float (int n, double *v);

/*
extern void shellSort_int_int(const int n, int *x, int *y);
extern void shellSort_int_float(int n, int *x, double *v);
extern void shellSort_int_int_float(int n, int *x, int *y, double *v);
*/
#ifdef __cplusplus
}
#endif
#endif
