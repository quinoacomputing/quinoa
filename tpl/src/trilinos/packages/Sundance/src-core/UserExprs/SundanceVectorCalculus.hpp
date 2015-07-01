/* @HEADER@ */
// ************************************************************************
// 
//                             Sundance
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Kevin Long (kevin.long@ttu.edu)
// 

/* @HEADER@ */

#ifndef SUNDANCE_VECTORCALCULUS_H
#define SUNDANCE_VECTORCALCULUS_H

#include "SundanceExpr.hpp"

namespace Sundance
{


/** \relates Expr */
Expr gradient(int dim);
  
/** \relates Expr */
Expr div(const Expr& f);
  
/** \relates Expr */
Expr cross(const Expr& a, const Expr& b);
  
/** \relates Expr */
Expr curl(const Expr& f);

/** \relates Expr 
 * Compute the colon (Frobenius) product of two matrices. The result
 * is a scalar. 
 **/
Expr colonProduct(const Expr& A, const Expr& B);

/** \relates Expr 
 * Compute the outer (Kronecker) product of two vectors. The result
 * is a dim(a) by dim(b) rectangular matrix.
 **/
Expr outerProduct(const Expr& a, const Expr& b);

/** \relates Expr
 * Indicate whether the given expression is a square matrix. If so,
 * return by reference argument the size of the matrix.
 *
 * A scalar is a 1 by 1 matrix.
 */
bool isSquareMatrix(const Expr& x, int& N);

/** \relates Expr
 * Indicate whether the given expression is a vector. If so,
 * return by reference argument the size of the vector.
 *
 * A scalar is a vector of dimension 1.
 * */
bool isVector(const Expr& x, int& N);

}

#endif
