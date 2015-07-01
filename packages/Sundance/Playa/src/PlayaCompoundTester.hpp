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


#ifndef PLAYA_COMPOUNDTESTER_HPP
#define PLAYA_COMPOUNDTESTER_HPP

#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaSimpleComposedOpDecl.hpp"
#include "PlayaSimpleScaledOpDecl.hpp"
#include "PlayaSimpleAddedOpDecl.hpp"
#include "PlayaSimpleDiagonalOpDecl.hpp"
#include "PlayaTesterBase.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "PlayaLinearCombinationImpl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaSimpleComposedOpImpl.hpp"
#include "PlayaSimpleTransposedOpImpl.hpp"
#include "PlayaSimpleScaledOpImpl.hpp"
#include "PlayaSimpleAddedOpImpl.hpp"
#include "PlayaSimpleDiagonalOpImpl.hpp"
#include "PlayaRandomSparseMatrixBuilderImpl.hpp"
#endif

using namespace Playa;
using namespace Teuchos;



namespace Playa
{

/** */
template <class Scalar>
class CompoundTester : public TesterBase<Scalar>
{
public:
  /** \brief Local typedef for promoted scalar magnitude */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** */
  CompoundTester(const LinearOperator<Scalar>& A,
    const LinearOperator<Scalar>& B,
    const TestSpecifier<Scalar>& sumSpec,
    const TestSpecifier<Scalar>& composedSpec,
    const TestSpecifier<Scalar>& scaledSpec,
    const TestSpecifier<Scalar>& diagSpec);

  /** */
  bool runAllTests() const ;

  /** */
  bool sumTest() const ;

  /** */
  bool composedTest() const ;

  /** */
  bool scaledTest() const ;

  /** */
  bool diagTest() const ;


private:

  LinearOperator<Scalar> A_;

  LinearOperator<Scalar> B_;

  TestSpecifier<Scalar> sumSpec_;

  TestSpecifier<Scalar> composedSpec_;

  TestSpecifier<Scalar> scaledSpec_;

  TestSpecifier<Scalar> diagSpec_;

};

template <class Scalar> 
inline CompoundTester<Scalar>
::CompoundTester(const LinearOperator<Scalar>& A,
  const LinearOperator<Scalar>& B,
  const TestSpecifier<Scalar>& sumSpec,
  const TestSpecifier<Scalar>& composedSpec,
  const TestSpecifier<Scalar>& scaledSpec,
  const TestSpecifier<Scalar>& diagSpec)
  : TesterBase<Scalar>(), 
    A_(A),
    B_(B),
    sumSpec_(sumSpec),
    composedSpec_(composedSpec),
    scaledSpec_(scaledSpec),
    diagSpec_(diagSpec)
{;}

template <class Scalar> 
inline bool CompoundTester<Scalar>
::runAllTests() const
{
  bool pass = true;

  pass = sumTest() && pass;
  pass = composedTest() && pass;
  pass = scaledTest() && pass;
  pass = diagTest() && pass;

  return pass;
}

template <class Scalar> 
inline bool CompoundTester<Scalar>
::sumTest() const 
{
  if (sumSpec_.doTest())
  {
    Out::root() << "running operator addition test..." << std::endl;
    LinearOperator<Scalar> sum = A_ + B_;

    Vector<Scalar> x = A_.domain().createMember();
    this->randomizeVec(x);
    Out::root() << "computing y1 = (A+B)*x..." << std::endl;
    Vector<Scalar> y1 = sum*x;
    Out::root() << "computing y2 = A*x + B*x..." << std::endl;
    Vector<Scalar> y2 = A_*x + B_*x;
    
    ScalarMag err = (y1 - y2).norm2();

    Out::root() << "|y1-y2| = " << err << std::endl;
        
    return this->checkTest(sumSpec_, err, "operator addition");
  }
  Out::root() << "skipping operator addition test..." << std::endl;
  return true;
}


template <class Scalar> 
inline bool CompoundTester<Scalar>
::composedTest() const 
{
  if (composedSpec_.doTest())
  {
    Out::root() << "running operator composition test..." << std::endl;
    LinearOperator<Scalar> composed = A_ * B_;

    Vector<Scalar> x = B_.domain().createMember();
    this->randomizeVec(x);
    Out::root() << "computing y1 = (A*B)*x..." << std::endl;
    Vector<Scalar> y1 = composed*x;
    Out::root() << "computing y2 = B*x..." << std::endl;
    Vector<Scalar> y2 = B_*x;
    Out::root() << "computing y3 = A*y2..." << std::endl;
    Vector<Scalar> y3 = A_*y2;

    ScalarMag err = (y1 - y3).norm2();

    Out::root() << "|y1-y3| = " << err << std::endl;
    return this->checkTest(composedSpec_, err, "operator composition");
  }
  Out::root() << "skipping operator composition test..." << std::endl;
  return true;
}


template <class Scalar> 
inline bool CompoundTester<Scalar>
::scaledTest() const 
{
  if (scaledSpec_.doTest())
  {
    Out::root() << "running operator scaling test..." << std::endl;
    Scalar alpha = sqrt(2.0);
    LinearOperator<Scalar> scaled = alpha*A_;

    Vector<Scalar> x = A_.domain().createMember();
    this->randomizeVec(x);
    Out::root() << "computing y1 = (alpha*A)*x..." << std::endl;
    Vector<Scalar> y1 = scaled*x;
    Out::root() << "computing y2 = A*x..." << std::endl;
    Vector<Scalar> y2 = A_*x;
    Out::root() << "computing y3 = alpha*y2..." << std::endl;
    Vector<Scalar> y3 = alpha*y2;

    ScalarMag err = (y1 - y3).norm2();

    Out::root() << "|y1-y3| = " << err << std::endl;
    return this->checkTest(composedSpec_, err, "operator scaling");
  }
  Out::root() << "skipping operator scaling test..." << std::endl;
  return true;
}

  

template <class Scalar> 
inline bool CompoundTester<Scalar>
::diagTest() const 
{
  if (diagSpec_.doTest())
  {
    Out::root() << "running diagonal operator test..." << std::endl;

    Vector<Scalar> x = A_.domain().createMember();
    this->randomizeVec(x);

    Vector<Scalar> d = A_.domain().createMember();
    this->randomizeVec(d);
        
    LinearOperator<Scalar> D = diagonalOperator(d);

    Out::root() << "computing y1 = D*x..." << std::endl;
    Vector<Scalar> y1 = D*x;
    Out::root() << "computing y2 = d .* x..." << std::endl;
    Vector<Scalar> y2 = x.dotStar(d);

    ScalarMag err = (y1 - y2).norm2();

    Out::root() << "|y1-y2| = " << err << std::endl;
    return this->checkTest(diagSpec_, err, "diagonal operator");
  }
  Out::root() << "skipping diagonal operator test..." << std::endl;
  return true;
}

  
  
}
#endif
