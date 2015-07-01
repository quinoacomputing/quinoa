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


#include "PySundanceLinearSolver.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaLinearSolverImpl.hpp"
#endif


using namespace Playa;


namespace Playa
{
  SolverState<double> 
  PySundanceLinearSolver_solve(const PySundanceLinearSolver* solver,
                               const LinearOperator<double>& op,
                               const Vector<double>& rhs,
                               Vector<double>& soln);
}

PySundanceLinearSolver::PySundanceLinearSolver(PyObject* functor) 
  : LinearSolverBase<double>(ParameterList()), 
    py_functor_(functor), py_solve_()

{
  // Increment the reference count
  Py_XINCREF(py_functor_);

  // If the python object has a "solve" attribute, set it
  // to be the PySundanceLinearSolver solve() callback function
  if (PyObject_HasAttrString (py_functor_,
			      "solve")) {
    setSolve(PyObject_GetAttrString(py_functor_,
                                    "solve"));
  }
}

PySundanceLinearSolver::~PySundanceLinearSolver() {
  // Decrement the reference count
  Py_XDECREF(py_functor_);
}


SolverState<double> PySundanceLinearSolver::solve(const LinearOperator<double>& op,
                                                  const Vector<double>& rhs,
                                                  Vector<double>& soln) const
{
  return PySundanceLinearSolver_solve(this, op, rhs, soln);
}


PyObject * PySundanceLinearSolver::setSolve(PyObject * p_pyObject)
{
  return py_solve_.setFunction(p_pyObject);
}

PyObject* PySundanceLinearSolver::pySolve(PyObject* opObj, PyObject* rhsObj, 
                                          PyObject* solnObj) const 
{
  PyObject* arglist = Py_BuildValue("(OOO)", opObj, rhsObj, solnObj);

  PyObject* result = PyEval_CallObject(py_solve_.getFunction(), arglist);
  
  Py_DECREF(arglist);  // All done with argument list

  return result;
}
