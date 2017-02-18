// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_SCALARMINIMIZATIONTEST_H
#define ROL_SCALARMINIMIZATIONTEST_H

/** \class ROL::ScalarMinimizationTest
    \brief Tests the minimization of scalar functions.
*/

#include "ROL_BrentsScalarMinimization.hpp"
#include "ROL_BisectionScalarMinimization.hpp"
#include "ROL_GoldenSectionScalarMinimization.hpp"
#include "ROL_ScalarFunction.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include <iostream>

namespace ROL { 

template<class Real>
class ScalarMinimizationTest {
private:
  Teuchos::RCP<ScalarMinimization<Real> > algo_;

public:
  virtual ~ScalarMinimizationTest(void) {}

  ScalarMinimizationTest(Teuchos::ParameterList &parlist) {
    std::string type = parlist.sublist("Scalar Minimization").get("Type","Brent's");
    if ( type == "Brent's" ) {
      algo_ = Teuchos::rcp(new BrentsScalarMinimization<Real>(parlist));
    }
    else if ( type == "Bisection" ) {
      algo_ = Teuchos::rcp(new BisectionScalarMinimization<Real>(parlist));
    }
    else if ( type == "Golden Section" ) {
      algo_ = Teuchos::rcp(new GoldenSectionScalarMinimization<Real>(parlist));
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> (ROL::ScalarMinimizationTest): Undefined ScalarMinimization type!");
    }
  }

  virtual bool test(std::ostream &stream = std::cout) = 0;

protected:
  void run(Real &fx, Real &x, int &nfval, int &ngrad,
            ScalarFunction<Real> &f,
           const Real &A, const Real &B) {
    algo_->run(fx, x, nfval, ngrad, f, A, B);
  }
};

}

#endif
