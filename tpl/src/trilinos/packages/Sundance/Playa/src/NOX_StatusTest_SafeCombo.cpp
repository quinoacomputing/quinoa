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


#include "PlayaDefs.hpp"
#include "NOX_StatusTest_SafeCombo.hpp"
#include "NOX_Utils.H"

using namespace Teuchos;
using std::runtime_error;
using std::cout;
using std::ostream;
using std::vector;



NOX::StatusTest::SafeCombo::SafeCombo(ComboType t) :
  type(t)
{
  status = Unevaluated;
}

NOX::StatusTest::SafeCombo::SafeCombo(ComboType t, 
                                      const Teuchos::RCP<Generic>& a) :
  type(t)
{
  tests.push_back(a);
  status = Unevaluated;
}

NOX::StatusTest::SafeCombo::SafeCombo(ComboType t, 
                                      const Teuchos::RCP<Generic>& a, 
                                      const Teuchos::RCP<Generic>& b) :
  type(t)
{
  tests.push_back(a);
  addStatusTest(b);
  status = Unevaluated;
}

NOX::StatusTest::SafeCombo& NOX::StatusTest::SafeCombo::addStatusTest(const Teuchos::RCP<Generic>& a)
{
  if (isSafe(a))
    tests.push_back(a);
  else 
  {
    const int indent = 2;
    cout << "\n*** WARNING! ***\n";
    cout << "This combo test currently consists of the following:\n";
    this->print(cout, indent);
    cout << "Unable to add the following test:\n";
    a->print(cout, indent);
    cout << "\n";
  }
  return *this;
}

bool NOX::StatusTest::SafeCombo::isSafe(const Teuchos::RCP<Generic>& a)
{
  // Are we trying to add "this" to "this"? This would result in an infinite recursion.
  if (a.get() == this)
    return false;
  
  // Recursively test that we're not adding something that's already
  // in the list because that can also lead to infinite recursions.
  for (vector<Teuchos::RCP<Generic> >::iterator i = tests.begin(); 
       i != tests.end(); ++i) 
  {
    
    SafeCombo* ptr = dynamic_cast<SafeCombo*>((*i).get());
    if (ptr != NULL)
      if (!ptr->isSafe(a))
	return false;
  }

  // Otherwise, it's safe to add a to the list.
  return true;
}

NOX::StatusTest::SafeCombo::~SafeCombo()
{
}

NOX::StatusTest::StatusType NOX::StatusTest::SafeCombo::checkStatus(const Solver::Generic& problem)
{
#ifdef TRILINOS_6
  return checkStatusEfficiently(problem, NOX::StatusTest::Minimal);
#else
  return checkStatus(problem, NOX::StatusTest::Minimal);
#endif
}

#ifdef TRILINOS_6
NOX::StatusTest::StatusType NOX::StatusTest::SafeCombo
::checkStatusEfficiently(const Solver::Generic& problem, 
                         NOX::StatusTest::CheckType checkType)
#else
NOX::StatusTest::StatusType NOX::StatusTest::SafeCombo
::checkStatus(const Solver::Generic& problem, 
              NOX::StatusTest::CheckType checkType)
#endif
{
  if (type == OR)
    orOp(problem, checkType);
  else
    andOp(problem, checkType);

  return status;
}

NOX::StatusTest::StatusType NOX::StatusTest::SafeCombo::getStatus() const
{
  return status;
}

void NOX::StatusTest::SafeCombo::orOp(const Solver::Generic& problem, NOX::StatusTest::CheckType checkType)
{
  if (checkType == NOX::StatusTest::None)
    status = Unevaluated;
  else
    status = Unconverged;

  // Checks the status of each test. The first test it encounters, if
  // any, that is unconverged is the status that it sets itself too.
  for (vector<Teuchos::RCP<Generic> >::const_iterator i = tests.begin(); i != tests.end(); ++i) 
  {
#ifdef TRILINOS_6
    NOX::StatusTest::StatusType s = (*i)->checkStatusEfficiently(problem, checkType);
#else
    NOX::StatusTest::StatusType s = (*i)->checkStatus(problem, checkType);
#endif
    if ((status == Unconverged) && (s != Unconverged)) 
    {
      status = s;

      // Turn off checking for the remaining tests
      if (checkType == NOX::StatusTest::Minimal)
        checkType = NOX::StatusTest::None;
    }

  }

  return;
}

void NOX::StatusTest::SafeCombo::andOp(const Solver::Generic& problem, NOX::StatusTest::CheckType checkType)
{
  if (checkType == NOX::StatusTest::None)
    status = Unevaluated;
  else
    status = Unconverged;

  bool isUnconverged = false;

  for (vector<Teuchos::RCP<Generic> >::const_iterator i = tests.begin(); i != tests.end(); ++i) {


#ifdef TRILINOS_6
    NOX::StatusTest::StatusType s = (*i)->checkStatusEfficiently(problem, checkType);
#else
    NOX::StatusTest::StatusType s = (*i)->checkStatus(problem, checkType);
#endif

    // If any of the tests are unconverged, then the AND test is
    // unconverged.
    if (s == Unconverged) 
    {
      isUnconverged = true;
      status = Unconverged;

      // Turn off checking for the remaining tests
      if (checkType == NOX::StatusTest::Minimal)
	checkType = NOX::StatusTest::None;
    }

    // If this is the first test and it's converged/failed, copy its
    // status to the combo status.
    if ((!isUnconverged) && (status == Unconverged)) 
    {
      status = s;
    }

  }

  return;
}


ostream& NOX::StatusTest::SafeCombo::print(std::ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
//   stream << setiosflags(ios::left) << setw(13) << setfill('.');
//   if (status == Unconverged) 
//     stream << "**";
//   else if (status == Failed)
//     stream << "Failed";
//   else
//     stream << "Converged";
  stream << status;
  stream << ((type == OR) ? "OR" : "AND");
  stream << " Combination";
  stream << " -> " << std::endl;

  for (vector<Teuchos::RCP<Generic> >::const_iterator i = tests.begin(); i != tests.end(); ++i) 
    (*i)->print(stream, indent+2);
    
  return stream;
}
