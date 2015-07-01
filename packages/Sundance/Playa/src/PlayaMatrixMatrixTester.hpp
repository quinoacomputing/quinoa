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


#ifndef PLAYA_MATRIXMATRIXTESTER_HPP
#define PLAYA_MATRIXMATRIXTESTER_HPP

#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaEpetraMatrixMatrixProduct.hpp"
#include "PlayaEpetraMatrixMatrixSum.hpp"
#include "PlayaEpetraMatrixOps.hpp"
#include "PlayaEpetraMatrix.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "PlayaLinearCombinationImpl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaSimpleComposedOpImpl.hpp"
#include "PlayaSimpleScaledOpImpl.hpp"
#include "PlayaSimpleAddedOpImpl.hpp"
#include "PlayaSimpleDiagonalOpImpl.hpp"
#include "PlayaSimpleTransposedOpImpl.hpp"
#include "PlayaRandomSparseMatrixBuilderImpl.hpp"
#endif

using namespace Playa;
using namespace Teuchos;



namespace Playa
{

/** */
template <class Scalar>
class MatrixMatrixTester : public TesterBase<Scalar>
{
public:
  /** \brief Local typedef for promoted scalar magnitude */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** */
  MatrixMatrixTester(const LinearOperator<Scalar>& A,
    const LinearOperator<Scalar>& B,
    const TestSpecifier<Scalar>& sumSpec,
    const TestSpecifier<Scalar>& prodSpec,
    const TestSpecifier<Scalar>& diagSpec,
    const TestSpecifier<Scalar>& diagLeftProdSpec,
    const TestSpecifier<Scalar>& diagRightProdSpec);

  /** */
  bool runAllTests() const ;

  /** */
  bool sumTest() const ;

  /** */
  bool prodTest() const ;

  /** */
  bool diagTest() const ;

  /** */
  bool diagLeftProdTest() const ;

  /** */
  bool diagRightProdTest() const ;


private:

  LinearOperator<Scalar> A_;

  LinearOperator<Scalar> B_;

  TestSpecifier<Scalar> sumSpec_;

  TestSpecifier<Scalar> prodSpec_;

  TestSpecifier<Scalar> diagSpec_;

  TestSpecifier<Scalar> diagLeftProdSpec_;

  TestSpecifier<Scalar> diagRightProdSpec_;

};

template <class Scalar> 
inline MatrixMatrixTester<Scalar>
::MatrixMatrixTester(const LinearOperator<Scalar>& A,
  const LinearOperator<Scalar>& B,
  const TestSpecifier<Scalar>& sumSpec,
  const TestSpecifier<Scalar>& prodSpec,
  const TestSpecifier<Scalar>& diagSpec,
  const TestSpecifier<Scalar>& diagRightProdSpec,
  const TestSpecifier<Scalar>& diagLeftProdSpec)
  : TesterBase<Scalar>(), 
    A_(A),
    B_(B),
    sumSpec_(sumSpec),
    prodSpec_(prodSpec),
    diagSpec_(diagSpec),
    diagLeftProdSpec_(diagLeftProdSpec),
    diagRightProdSpec_(diagRightProdSpec)
{;}

template <class Scalar> 
inline bool MatrixMatrixTester<Scalar>
::runAllTests() const
{
  bool pass = true;

  pass = this->sumTest() && pass;
  pass = this->prodTest() && pass;
  pass = this->diagTest() && pass;
  pass = this->diagLeftProdTest() && pass;
  pass = this->diagRightProdTest() && pass;

  return pass;
}

template <class Scalar> 
inline bool MatrixMatrixTester<Scalar>
::diagTest() const 
{
  if (diagSpec_.doTest())
  {
    Vector<Scalar> x = B_.domain().createMember();
    Vector<Scalar> d = B_.domain().createMember();
    this->randomizeVec(x);
    this->randomizeVec(d);
    LinearOperator<Scalar> D0 = diagonalOperator(d);
    LinearOperator<Scalar> D = makeEpetraDiagonalMatrix(d);
    Vector<Scalar> d1 = getEpetraDiagonal(D);

    Out::root() << "computing implicit product y1 = D*x..." << std::endl;
    Vector<Scalar> y1 = D0*x;
    Out::root() << "computing explicit product y2 = D*x..." << std::endl;
    Vector<Scalar> y2 = D*x;

    ScalarMag err = (y1 - y2).norm2();

    Out::root() << "|y1-y2| = " << err << std::endl;
    
    Out::root() << "comparing recovered and original diagonals" << std::endl;
    ScalarMag err2 = (d - d1).norm2();
    Out::root() << "|d1-d2| = " << err2 << std::endl;
    
    return this->checkTest(prodSpec_, err+err2, "matrix-matrix multiply");
    
  }
  Out::root() << "skipping matrix-matrix multiply test..." << std::endl;
  return true;
}



template <class Scalar> 
inline bool MatrixMatrixTester<Scalar>
::sumTest() const 
{
  if (sumSpec_.doTest())
  {
    /* skip incompatible matrices. This will occur when we're testing
     * multiplication of rectangular matrices */
    if (A_.range() != B_.range() || A_.domain() != B_.domain())
    {
      Out::root() << "skipping sum on incompatible matrices" << std::endl;
      return true;
    }
    /* If here, the sum should work */
    Out::root() << "running matrix-matrix multiply test..." << std::endl;
    LinearOperator<Scalar> implicitAdd = A_ + B_;
    LinearOperator<Scalar> explicitAdd = epetraMatrixMatrixSum(A_, B_);

    Vector<Scalar> x = B_.domain().createMember();
    this->randomizeVec(x);
    Out::root() << "computing implicit sum y1 = (A+B)*x..." << std::endl;
    Vector<Scalar> y1 = implicitAdd*x;
    Out::root() << "computing explicit sum y2 = (A+B)*x..." << std::endl;
    Vector<Scalar> y2 = explicitAdd*x;

    ScalarMag err = (y1 - y2).norm2();

    Out::root() << "|y1-y2| = " << err << std::endl;
    return this->checkTest(prodSpec_, err, "matrix-matrix multiply");
    
  }
  Out::root() << "skipping matrix-matrix multiply test..." << std::endl;
  return true;
}


template <class Scalar> 
inline bool MatrixMatrixTester<Scalar>
::prodTest() const 
{
  if (prodSpec_.doTest())
  {
    Out::root() << "running matrix-matrix multiply test..." << std::endl;
    LinearOperator<Scalar> composed = A_ * B_;
    LinearOperator<Scalar> multiplied = epetraMatrixMatrixProduct(A_, B_);

    Vector<Scalar> x = B_.domain().createMember();
    this->randomizeVec(x);
    Out::root() << "computing implicit product y1 = (A*B)*x..." << std::endl;
    Vector<Scalar> y1 = composed*x;
    Out::root() << "computing explicit product y2 = (A*B)*x..." << std::endl;
    Vector<Scalar> y2 = multiplied*x;

    ScalarMag err = (y1 - y2).norm2();

    Out::root() << "|y1-y2| = " << err << std::endl;
    return this->checkTest(prodSpec_, err, "matrix-matrix multiply");
  }
  Out::root() << "skipping matrix-matrix multiply test..." << std::endl;
  return true;
}


template <class Scalar> 
inline bool MatrixMatrixTester<Scalar>
::diagLeftProdTest() const 
{
  if (diagLeftProdSpec_.doTest())
  {
    Out::root() << "running diagonal*matrix multiplication test..." << std::endl;

    Vector<Scalar> x = A_.domain().createMember();
    this->randomizeVec(x);

    Vector<Scalar> d = A_.range().createMember();
    this->randomizeVec(d);
        
    LinearOperator<Scalar> D = diagonalOperator(d);
    LinearOperator<Scalar> DA = epetraLeftScale(d, A_);

    Out::root() << "computing implicit y1 = D*A*x..." << std::endl;
    Vector<Scalar> y1 = D*A_*x;
    Out::root() << "computing explicit y2 = D*A*x..." << std::endl;
    Vector<Scalar> y2 = DA*x;

    ScalarMag err = (y1 - y2).norm2();

    Out::root() << "|y1-y2| = " << err << std::endl;

    return this->checkTest(diagLeftProdSpec_, err, "diagonal*matrix multiplication");
  }
  Out::root() << "skipping diagonal matrix-matrix test..." << std::endl;
  return true;
}

  
template <class Scalar> 
inline bool MatrixMatrixTester<Scalar>
::diagRightProdTest() const 
{
  if (diagRightProdSpec_.doTest())
  {
    Out::root() << "running diagonal*matrix multiplication test..." << std::endl;

    Vector<Scalar> x = A_.domain().createMember();
    this->randomizeVec(x);

    Vector<Scalar> d = A_.domain().createMember();
    this->randomizeVec(d);
        
    LinearOperator<Scalar> D = diagonalOperator(d);
    LinearOperator<Scalar> AD = epetraRightScale(A_, d);

    Out::root() << "computing implicit y1 = A*D*x..." << std::endl;
    Vector<Scalar> y1 = A_*D*x;
    Out::root() << "computing explicit y2 = A*D*x..." << std::endl;
    Vector<Scalar> y2 = AD*x;

    ScalarMag err = (y1 - y2).norm2();

    Out::root() << "|y1-y2| = " << err << std::endl;

    return this->checkTest(diagLeftProdSpec_, err, "matrix*diagonal multiplication");
  }
  Out::root() << "skipping diagonal matrix-matrix test..." << std::endl;
  return true;
}

  
  
}
#endif
