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


#ifndef PLAYA_LINEARCOMBINATIONTESTER_HPP
#define PLAYA_LINEARCOMBINATIONTESTER_HPP

#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaSimpleComposedOpDecl.hpp"
#include "PlayaSimpleScaledOpDecl.hpp"
#include "PlayaSimpleAddedOpDecl.hpp"
#include "PlayaSimpleIdentityOpDecl.hpp"
#include "PlayaOut.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaTesterBase.hpp"
#include "PlayaRandomSparseMatrixBuilderDecl.hpp"
#include "PlayaTestSpecifier.hpp"
#include "Teuchos_ScalarTraits.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaRandomSparseMatrixBuilderImpl.hpp"
#include "PlayaSimpleComposedOpImpl.hpp"
#include "PlayaSimpleScaledOpImpl.hpp"
#include "PlayaSimpleAddedOpImpl.hpp"
#include "PlayaSimpleTransposedOpImpl.hpp"
#include "PlayaSimpleIdentityOpDecl.hpp"
#endif

using namespace Playa;
using namespace Teuchos;





#define TESTER(form1, form2)\
  {\
    Out::os() << "testing " #form1 << std::endl;\
    Vector<Scalar> _val1 = form1;\
    Out::os() << "testing " #form2 << std::endl;\
    Vector<Scalar> _val2 = form2;\
    Out::os() << "done testing... checking error" << std::endl;\
    ScalarMag err = (_val1-_val2).norm2();\
    if (!this->checkTest(spec_, err, "[" #form1 "] == [" #form2 "]")) pass = false;\
  }
    






namespace Playa
{

/** */
template <class Scalar>
class LinearCombinationTester : public TesterBase<Scalar>
{
public:
  /** \brief Local typedef for promoted scalar magnitude */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** */
  LinearCombinationTester(int nLocalRows, double onProcDensity, 
    double offProcDensity,
    const VectorType<Scalar>& vecType,
    const TestSpecifier<Scalar>& spec);

  /** */
  bool runAllTests() const ;


private:

  /** */
  bool nonModifyingOpTests() const ;

  /** */
  bool selfModifyingOpTests() const ;

  TestSpecifier<Scalar> spec_;

  int nLocalRows_;
  double onProcDensity_;
  double offProcDensity_;
  VectorType<Scalar> vecType_;

};

template <class Scalar> 
inline LinearCombinationTester<Scalar>
::LinearCombinationTester(int nLocalRows, double onProcDensity, 
  double offProcDensity,
  const VectorType<Scalar>& vecType,
  const TestSpecifier<Scalar>& spec)
  : TesterBase<Scalar>(),
    spec_(spec),
    nLocalRows_(nLocalRows),
    onProcDensity_(onProcDensity),
    offProcDensity_(offProcDensity),
    vecType_(vecType)
{;}

template <class Scalar> 
inline bool LinearCombinationTester<Scalar>
::runAllTests() const
{
  bool pass = true;

  pass = nonModifyingOpTests() && pass;
  pass = selfModifyingOpTests() && pass;

  return pass;
}

template <class Scalar> 
inline bool LinearCombinationTester<Scalar>
::nonModifyingOpTests() const
{
  bool pass = true;

  RandomSparseMatrixBuilder<double> ABuilder(nLocalRows_, nLocalRows_, 
    onProcDensity_, offProcDensity_, 
    vecType_);
  RandomSparseMatrixBuilder<double> BBuilder(nLocalRows_, nLocalRows_, 
    onProcDensity_, offProcDensity_, 
    vecType_);
  RandomSparseMatrixBuilder<double> CBuilder(nLocalRows_, nLocalRows_, 
    onProcDensity_, offProcDensity_, 
    vecType_);


  LinearOperator<double> A = ABuilder.getOp();
  LinearOperator<Scalar> B = BBuilder.getOp();
  LinearOperator<Scalar> C = CBuilder.getOp();

  Vector<Scalar> x = A.domain().createMember();
  Vector<Scalar> y = A.domain().createMember();
  Vector<Scalar> z = A.domain().createMember();

  this->randomizeVec(x);
  this->randomizeVec(y);
  this->randomizeVec(z);


  TESTER(x*2.0, 2.0*x);

  TESTER(2.0*(x + y), 2.0*x + 2.0*y);

  TESTER(2.0*(x - y), 2.0*x - 2.0*y);

  TESTER((2.0*x + y) - y, 2.0*x);

  TESTER(-1.0*y + (2.0*x + y), 2.0*x);

  TESTER((x + y)*2.0, 2.0*x + 2.0*y);

  TESTER((x - y)*2.0, 2.0*x - 2.0*y);

  TESTER(2.0*(x - y), -2.0*(y - x));

  TESTER(0.25*(2.0*(x + y) - 2.0*(x - y)), y);

  TESTER((2.0*A)*x, 2.0*(A*x));

  TESTER(2.0*(A*x), (A*x)*2.0);

  TESTER(A*(B*x), (A*B)*x);

  TESTER(2.0*A*(B*x), A*(B*(2.0*x)));

  TESTER(3.0*(2.0*A)*x, 6.0*(A*x));

  TESTER(y + A*x, A*x + y);


  TESTER(z + (A*x + B*y), (B*y + A*x) + z);


  TESTER(z - (A*x + B*y), -1.0*((B*y + A*x) - z));

  TESTER(C*z + (A*x + B*y), (B*y + A*x) + C*z);

  TESTER(C*z - (A*x + B*y), -1.0*((B*y + A*x) - C*z));

  TESTER(2.0*z + (A*x + B*y), (B*y + A*x) + 2.0*z);

  TESTER(2.0*z - (A*x + B*y), -1.0*((B*y + A*x) - 2.0*z));

  TESTER(A*x - y, -1.0*(y - A*x));

  TESTER(A*x + B*y, B*y + A*x);

  TESTER(A*x - B*y - A*x + B*y +z, z);

  TESTER(2.0*(A*x + y), 2.0*A*x + 2.0*y);

  TESTER(2.0*(A*x + B*y), A*x + B*y + A*x + B*y);

  TESTER(2.0*(y + A*x), 2.0*y + 2.0*(A*x));

  TESTER(x + 2.0*A*y, x + 2.0*(A*y));

  TESTER(2.0*A*y + x, 2.0*(A*y) + x);

  TESTER(2.0*A*(3.0*B)*y, 6.0*(A*B*y));

  TESTER(2.0*A*(3.0*B)*y, 6.0*(A*(B*y)));

  TESTER(2.0*A*(3.0*B + 2.0*A)*x, 6.0*A*B*x + 4.0*A*A*x );

  TESTER(2.0*(A*x + B*y), 2.0*A*x + 2.0*B*y);

  TESTER(2.0*(A*x - B*y), 2.0*A*x - 2.0*B*y);

  TESTER(2.0*(A*x + B*y + z), 2.0*A*x + 2.0*B*y + 2.0*z);

  TESTER(2.0*(A*x + 3.0*B*y), 2.0*A*x + 6.0*B*y);

  TESTER(2.0*(A*x + 3.0*(z + B*y)), 2.0*A*x + 6.0*B*y + 6.0*z);

  TESTER(2.0*(z + A*x + B*y + z), 2.0*A*x + 2.0*B*y + 4.0*z);

  TESTER(2.0*(3.0*(z + A*x) + B*y), 6.0*z + 6.0*A*x + 2.0*B*y);

  TESTER(2.0*(3.0*(z + A*x) + 4.0*(B*y + z)), 6.0*z + 6.0*A*x + 8.0*B*y + 8.0*z);
    
  TESTER((A*x + B*y) + (A*y + B*x), (A + B)*x + (A+B)*y);
  TESTER((A*x + B*y) - (A*y + B*x), A*x - A*y + B*y - B*x);

  TESTER((A*x + B*y) + 2.0*(A*y + B*x), A*(x + 2.0*y) + B*(2.0*x + y));
  TESTER((A*x + B*y) - 2.0*(A*y + B*x), A*(x - 2.0*y) + B*(y - 2.0*x));


  return pass;
}


template <class Scalar> 
inline bool LinearCombinationTester<Scalar>
::selfModifyingOpTests() const
{
  bool pass = true;

  RandomSparseMatrixBuilder<double> ABuilder(nLocalRows_, nLocalRows_, 
    onProcDensity_, offProcDensity_, 
    vecType_);
  RandomSparseMatrixBuilder<double> BBuilder(nLocalRows_, nLocalRows_, 
    onProcDensity_, offProcDensity_, 
    vecType_);
  RandomSparseMatrixBuilder<double> CBuilder(nLocalRows_, nLocalRows_, 
    onProcDensity_, offProcDensity_, 
    vecType_);


  LinearOperator<double> A = ABuilder.getOp();
  LinearOperator<double> B = BBuilder.getOp();
  LinearOperator<double> C = CBuilder.getOp();

  VectorSpace<Scalar> vs = A.domain();
  A = identityOperator(vs);
  B = identityOperator(vs);
  C = identityOperator(vs);

  Vector<Scalar> x = A.domain().createMember();
  Vector<Scalar> y = A.domain().createMember();
  Vector<Scalar> z = A.domain().createMember();

  this->randomizeVec(x);
  this->randomizeVec(y);
  this->randomizeVec(z);

  Vector<Scalar> a = x.copy();
  Vector<Scalar> b = y.copy();
  Vector<Scalar> c = z.copy();
    

  Out::os() << "starting linear combination tests" << std::endl;

  x = 2.0*A*x;
  Scalar err = (x - 2.0*A*a).norm2();
  if (!this->checkTest(spec_, err, "x=2.0*A*x")) pass = false;

  a = x.copy();
  x = x + y;
  err = (x - (a + y)).norm2();
  if (!this->checkTest(spec_, err, "x=x+y")) pass = false;

  a = x.copy();
  x = y + x;
  err = (x - (y + a)).norm2();
  if (!this->checkTest(spec_, err, "x=y+x")) pass = false;

  a = x.copy();
  x = A*x + B*y;
  err = (x - (A*a + B*y)).norm2();
  if (!this->checkTest(spec_, err, "x=A*x+B*y")) pass = false;

  a = x.copy();
  x = B*y + A*x;
  err = (x - (B*y + A*a)).norm2();
  if (!this->checkTest(spec_, err, "x=B*y+A*x")) pass = false;

  a = x.copy();
  x = A*x + (B*y + C*x);
  err = (x - (A*a + (B*y + C*a))).norm2();
  if (!this->checkTest(spec_, err, "x=A*x + (B*y + C*x)")) pass = false;

  a = x.copy();
  x = (A*x + B*y) + C*x;
  err = (x - ((A*a + B*y) + C*a)).norm2();
  if (!this->checkTest(spec_, err, "x=(A*x + B*y) + C*x")) pass = false;

  /* test assignment of OpTimesLC into empty and non-empty vectors */
  Vector<Scalar> u;
  u = 2.0*A*B*x;
  err = (u - 2.0*A*B*x).norm2();
  if (!this->checkTest(spec_, err, "(empty) u=2*A*B*x")) pass = false;

  u = 2.0*A*B*x;
  err = (u - 2.0*A*B*x).norm2();
  if (!this->checkTest(spec_, err, "(non-empty) u=2*A*B*x")) pass = false;

  /* test assignment of LC2 into empty and non-empty vectors */
  Vector<Scalar> v;
  v = 2.0*x + 3.0*y;
  err = (v - (2.0*x + 3.0*y)).norm2();
  if (!this->checkTest(spec_, err, "(empty) v=2*x + 3*y")) pass = false;

  v = 2.0*x + 3.0*y;
  err = (v - (2.0*x + 3.0*y)).norm2();
  if (!this->checkTest(spec_, err, "(non-empty) v=2*x + 3*y")) pass = false;

  /* test assignment of LC3 into empty and non-empty vectors */
  Vector<Scalar> w;
  w = 2.0*x + 3.0*y + 5.0*z;
  err = (w - (2.0*x + 3.0*y + 5.0*z)).norm2();
  if (!this->checkTest(spec_, err, "(empty) w=2*x + 3*y + 5*z")) pass = false;

  w = 2.0*x + 3.0*y + 5.0*z;
  err = (w - (2.0*x + 3.0*y + 5.0*z)).norm2();
  if (!this->checkTest(spec_, err, "(non-empty) w=2*x + 3*y + 5*z")) pass = false;

  /* test assignment of LC4 into empty and non-empty vectors */
  Vector<Scalar> w2;
  w2 = 2.0*x + 3.0*y + 5.0*z + 7.0*u;
  err = (w2 - (2.0*x + 3.0*y + 5.0*z + 7.0*u)).norm2();
  if (!this->checkTest(spec_, err, 
      "(empty) w2=2*x + 3*y + 5*z + 7*u")) pass = false;

  w2 = 2.0*x + 3.0*y + 5.0*z + 7.0*u;
  err = (w2 - (2.0*x + 3.0*y + 5.0*z + 7.0*u)).norm2();
  if (!this->checkTest(spec_, err, 
      "(non-empty) w2=2*x + 3*y + 5*z + 7*u")) pass = false;

  /* test assignment of LC3 into one of the operands */
  x = 2.0*x + 3.0*y + 5.0*z;
  err = (w - x).norm2(); // Note: w contains 2x + 3y + 5z 
  if (!this->checkTest(spec_, err, "x=2*x + 3*y + 5*z")) pass = false;


  return pass;
}

  
}
#endif
