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


#ifndef PLAYA_TESTERBASE_HPP
#define PLAYA_TESTERBASE_HPP

#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaTestSpecifier.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "PlayaOut.hpp"


using namespace Teuchos;



namespace Playa
{

/** */
template <class Scalar>
class TesterBase
{
public:
  /** \brief Local typedef for promoted scalar magnitude */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** */
  TesterBase(){;}

  /** */
  virtual ~TesterBase(){;}

  /** */
  virtual bool runAllTests() const = 0 ;


  /** */
  bool checkTest(const TestSpecifier<Scalar>& spec,
    const ScalarMag& err, 
    const std::string& testName) const ;

  /** */
  void randomizeVec(Vector<Scalar>& x) const ;

};

template <class Scalar> 
inline void TesterBase<Scalar>
::randomizeVec(Vector<Scalar>& x) const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  x.randomize();
    
}

template <class Scalar> 
inline bool TesterBase<Scalar>
::checkTest(const TestSpecifier<Scalar>& spec,
  const ScalarMag& err, 
  const std::string& testName) const 
{
  bool rtn = true;
  if (err > spec.errorTol())
  {
    Out::root() << testName << " test FAILED: err=" << err << ", tol = " 
                << spec.errorTol() << std::endl;
    rtn = false;
  }
  else if (err > spec.warningTol())
  {
    Out::root() << "WARNING: " << testName << " test err="
                << err << " could not beat tol = " 
                << spec.warningTol() << std::endl;
  }
  else
  {
    Out::root() << "test " << testName << " PASSED with tol=" << spec.errorTol() << std::endl;
  }
  return rtn;
}
  
}
#endif
