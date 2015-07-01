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


#include "PlayaDefs.hpp"

#ifdef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION

#include "PlayaVectorImpl.hpp"
#include "PlayaLinearCombinationImpl.hpp"

namespace Playa
{

template class Vector<double>;

template LoadableVector<double>* loadable(Vector<double> vec);

template 
double* dataPtr(Vector<double> vec) ;

template 
const double* dataPtr(const Vector<double>& vec) ;

template class LCN<double, 1>;
template class LCN<double, 2>;
template class LCN<double, 3>;
template class LCN<double, 4>;

template Vector<double> operator*(
  const LinearOperator<double>& A,
  const Vector<double>& x);

template Vector<double> operator*(
  const LinearOperator<double>& A,
  const LCN<double, 1>& x);

template Vector<double> operator*(
  const LinearOperator<double>& A,
  const LCN<double, 2>& x);

template Vector<double> operator*(
  const LinearOperator<double>& A,
  const LCN<double, 3>& x);

template Vector<double>& Vector<double>::operator+=(const LCN<double, 3>& x);
template Vector<double>& Vector<double>::operator-=(const LCN<double, 3>& x);

template double norm1(const LCN<double, 1>& x);
template double norm1(const LCN<double, 2>& x);
template double norm1(const LCN<double, 3>& x);

template double norm2(const LCN<double, 1>& x);
template double norm2(const LCN<double, 2>& x);
template double norm2(const LCN<double, 3>& x);

template double normInf(const LCN<double, 1>& x);
template double normInf(const LCN<double, 2>& x);
template double normInf(const LCN<double, 3>& x);

template double min(const LCN<double, 1>& x);
template double min(const LCN<double, 2>& x);
template double min(const LCN<double, 3>& x);

template double max(const LCN<double, 1>& x);
template double max(const LCN<double, 2>& x);
template double max(const LCN<double, 3>& x);

template Vector<double> abs(const LCN<double, 1>& x);
template Vector<double> abs(const LCN<double, 2>& x);
template Vector<double> abs(const LCN<double, 3>& x);

template Vector<double> reciprocal(const LCN<double, 1>& x);
template Vector<double> reciprocal(const LCN<double, 2>& x);
template Vector<double> reciprocal(const LCN<double, 3>& x);

template LCN<double, 1> operator*(const double& a, const Vector<double>& x);
template LCN<double, 1> operator*(const Vector<double>& x, const double& a);
template LCN<double, 1> operator/(const Vector<double>& x, const double& a);

template LCN<double, 1> operator*(const double& a, const LCN<double, 1>& x);
template LCN<double, 1> operator*(const LCN<double, 1>& x, const double& a);
template LCN<double, 1> operator/(const LCN<double, 1>& x, const double& a);

template LCN<double, 2> operator*(const double& a, const LCN<double, 2>& x);
template LCN<double, 2> operator*(const LCN<double, 2>& x, const double& a);
template LCN<double, 2> operator/(const LCN<double, 2>& x, const double& a);

template LCN<double, 3> operator*(const double& a, const LCN<double, 3>& x);
template LCN<double, 3> operator*(const LCN<double, 3>& x, const double& a);
template LCN<double, 3> operator/(const LCN<double, 3>& x, const double& a);

template LCN<double, 2> 
operator+(const Vector<double>& y, const LCN<double, 1>& x);
template LCN<double, 2> 
operator+(const LCN<double, 1>& x, const Vector<double>& y);
template LCN<double, 2> 
operator+(const LCN<double, 1>& x, const LCN<double, 1>& y);

template LCN<double, 2> 
operator-(const Vector<double>& y, const LCN<double, 1>& x);
template LCN<double, 2> 
operator-(const LCN<double, 1>& x, const Vector<double>& y);
template LCN<double, 2> 
operator-(const LCN<double, 1>& x, const LCN<double, 1>& y);


template LCN<double, 3> 
operator+(const Vector<double>& y, const LCN<double, 2>& x);
template LCN<double, 3> 
operator+(const LCN<double, 2>& x, const Vector<double>& y);
template LCN<double, 3> 
operator+(const LCN<double, 2>& x, const LCN<double, 1>& y);
template LCN<double, 3> 
operator+(const LCN<double, 1>& x, const LCN<double, 2>& y);

template LCN<double, 3> 
operator-(const Vector<double>& y, const LCN<double, 2>& x);
template LCN<double, 3> 
operator-(const LCN<double, 2>& x, const Vector<double>& y);
template LCN<double, 3> 
operator-(const LCN<double, 2>& x, const LCN<double, 1>& y);
template LCN<double, 3> 
operator-(const LCN<double, 1>& x, const LCN<double, 2>& y);



}

#endif
