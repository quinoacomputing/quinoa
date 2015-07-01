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


// $Id$ 
// $Source$ 


//   


#ifndef NOX_STATUSTEST_SAFECOMBO_H
#define NOX_STATUSTEST_SAFECOMBO_H

#include "NOX_StatusTest_Generic.H" // base class
#include "NOX_Common.H"		// class data element (vector)
#include "Teuchos_RefCountPtr.hpp"

namespace NOX {

  namespace StatusTest {

    /*!
      \brief Arbitrary combination of status tests, implemented with safe memory management.

      In the \c AND (see NOX::StatusTest::Combo::ComboType) combination, the
      result is \c Unconverged (see NOX::StatusTest::StatusType) if \e any of
      the tests is \c Unconverged. Otherwise, the result is equal to the
      result of the \e first test in the list that is either \c Converged
      or \c Failed. It is not recommended to mix \c Converged and \c
      Failed tests in an \c AND combination.

      In the \c OR combination, the result is \c Unconverged if \e all of
      the tests are \c Unconverged. Otherwise, it is the result of the \e
      first test in the list that is either \c Converged or \c
      Failed. Therefore, it will generally make sense to put the \c Failed
      -type tests at the end of the \c OR list.

      \note We call checkStatus on \e every convergence test, though some
      may be called with the NOX::StatusTest::None option.

      \note This is identical in operation to StatusTest::Combo. but stores
      the constituent subtests using a reference-counted pointer. This ensures
      safe behavior if the constituent subtests go out of scope before the combo.
      - KL 24 August 2004.

      \author Tammy Kolda (SNL 8950)
      \author Kevin Long (SNL 8962)
    */
    class SafeCombo : public Generic {

    public:

      /*! 
        \brief The test can be either the AND of all the component tests,
        or the OR of all the component tests.
      */
      enum ComboType {
        //! Logically "AND" together the results of the tests in this combination
        AND, 
        //! Logically "OR" together the results of the tests in this combination
        OR
      };

      //! Constructor
      SafeCombo(ComboType t);

      //! Constructor with a single test.
      SafeCombo(ComboType t, const Teuchos::RCP<Generic>& a);

      //! Constructor with two tests.
      SafeCombo(ComboType t, const Teuchos::RCP<Generic>& a,
            const Teuchos::RCP<Generic>& b);

      //! Add another test to this combination. 
      /*!
        Calls isSafe() to determine if it is safe to add \c a to the combination.
      */
      virtual SafeCombo& addStatusTest(const Teuchos::RCP<Generic>& a);

      //! Destructor
      virtual ~SafeCombo();

      /*!
        \brief Calls checkStatus(problem, NOX::StatusTest::Minimal)
      */
      virtual StatusType checkStatus(const NOX::Solver::Generic& problem);

      /*!
        \brief Tests stopping criterion.
    
        See addOp() and orOp() for details.
      */
#ifdef TRILINOS_6
      virtual StatusType checkStatusEfficiently(const NOX::Solver::Generic& problem, NOX::StatusTest::CheckType checkType);
#else
      virtual StatusType checkStatus(const NOX::Solver::Generic& problem, NOX::StatusTest::CheckType checkType);
#endif

      virtual StatusType getStatus() const;

      virtual std::ostream& print(std::ostream& stream, int indent = 0) const;
  
    protected:

      //! Use this for checkStatus when this is an OR type combo. Updates NOX::StatusTest::Combo::status.
      /*!  
        If there is a combination of NOX::StatusTest::Failed and
        NOX::StatusTest::Converged in the tests that are OR'd together,
        the value of status for this test is set to the status of the
        first test it encounters which is not NOX::Status::Unconvered. The
        tests are evaluated in the order that they were added to the
        combination.

        \note We compute the status of all tests in the combination for
        the sake of completeness, even if we could determine the status of
        this combination test without that check.

      */
      virtual void orOp(const Solver::Generic& problem, NOX::StatusTest::CheckType checkType);

      //! Use this for checkStatus when this is an AND type combo. Updates NOX::StatusTest::Combo::status.
      /*!  

      If any tests are NOX::StatusTest::Unconverged, then the status of
      this test is NOX::StatusTest::Unconverged.  If there is a
      combination of NOX::StatusTest::Failed and
      NOX::StatusTest::Converged in the tests that are AND'd together,
      the value of status for this test is set to the status of the
      first test it encounters.  The tests are evaluated in the
      order that they were added to the combination.

      \note We compute the status of all tests in the combination for
      the sake of completeness, even if we could determine the status of
      this combination test without that check.
      */
      virtual void andOp(const Solver::Generic& problem, NOX::StatusTest::CheckType checkType);

      /*! \brief Check whether or not it is safe to add \c a to this list
        of tests.

        This is necessary to avoid any infinite recursions 
        (i.e., a test cannot own a copy of itself).
      */
      bool isSafe(const Teuchos::RCP<Generic>& a);

    private:

      //! Type of test
      const ComboType type;

      //! Vector of generic status tests
      std::vector<Teuchos::RCP<Generic> > tests;

      //! %Status
      StatusType status;

    }; // class Combo

  } // namespace Status
} // namespace NOX

#endif
