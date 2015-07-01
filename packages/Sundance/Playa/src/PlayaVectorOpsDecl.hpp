/* @HEADER@ */
// ************************************************************************
// 
//                 Playa: Programmable Linear Algebra
//                 Copyright 2012 Sandia Corporation
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

#ifndef PLAYA_VECTOROPSDECL_HPP
#define PLAYA_VECTOROPSDECL_HPP

#include "PlayaDefs.hpp"

namespace Playa
{

template <class Scalar> class Vector;

/** \relates Vector \brief Return minimum element and its location */
template <class Scalar> 
Scalar minloc(const Vector<Scalar>& x, int& gni);

/** \relates Vector \brief Return maximum element and its location*/
template <class Scalar> 
Scalar maxloc(const Vector<Scalar>& x, int& gni);

/** \relates Vector \brief Return minimum element greater than a specified
bound, and its location*/
template <class Scalar> 
Scalar minlocWithBound(const Scalar& lowerBound, 
  const Vector<Scalar>& x, int& gni);

/** \relates Vector \brief Return maximum element less than a specified
bound, and its location*/
template <class Scalar> 
Scalar maxlocWithBound(const Scalar& upperBound, 
  const Vector<Scalar>& x, int& gni);

/** \relates Vector \brief Compute the Euclidean norm of a vector */
template <class Scalar>
Scalar norm2(const Vector<Scalar>& x);

/** \relates Vector \brief Compute the one-norm of a vector */
template <class Scalar>
Scalar norm1(const Vector<Scalar>& x);

/** \relates Vector \brief Compute the infinity norm of a vector */
template <class Scalar>
Scalar normInf(const Vector<Scalar>& x);

/** \relates Vector 
 * \brief Compute the Euclidean distance between two vectors */
template <class Scalar>
Scalar norm2Dist(const Vector<Scalar>& x, const Vector<Scalar>& y);

/** \relates Vector \brief Compute the one-norm distance between two vectors */
template <class Scalar>
Scalar norm1Dist(const Vector<Scalar>& x, const Vector<Scalar>& y);

/** \relates Vector
 *  \brief Compute the infinity-norm distance between two vectors */
template <class Scalar>
Scalar normInfDist(const Vector<Scalar>& x, const Vector<Scalar>& y);

}

 

#endif
