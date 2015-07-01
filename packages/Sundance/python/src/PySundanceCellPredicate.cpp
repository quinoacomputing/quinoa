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


#include "PySundanceCellPredicate.hpp"

using namespace Sundance;
using namespace Sundance;

PySundanceCellPredicate::PySundanceCellPredicate(PyObject* functor) 
  : py_functor_(functor), evalOpCallback_(), descrCallback_()

{
  // Increment the reference count
  Py_XINCREF(py_functor_);

  // If the python object has a "evalOp" attribute, set it
  // to be the PySundanceCellPredicate computeF callback function. If not,
  // we have an error!
  if (PyObject_HasAttrString (py_functor_, "evalOp")) 
    {
      setEvalOp(PyObject_GetAttrString(py_functor_, "evalOp"));
    } 
  else
    {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
                         "PySundanceCellPredicate bound to a Python object "
                         "without a method called evalOp().");
    }

  // If the python object has a "description" attribute, bind
  // it to the description() function call. Otherwise, use
  // a raw description of this object as a description. 
  if (PyObject_HasAttrString (py_functor_, "description")) 
    {
      setDescr(PyObject_GetAttrString(py_functor_, "description"));
    }
}

PySundanceCellPredicate::~PySundanceCellPredicate() {
  // Decrement the reference count
  Py_XDECREF(py_functor_);
}


bool PySundanceCellPredicate::operator()(const Point& x) const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(evalOpCallback_.get()==0, std::runtime_error,
                     "null pointer to python evalOp() method");

  PyObject * arglist;
  PyObject * result;

  switch(x.dim())
    {
    case 1:
      arglist = Py_BuildValue("(d)", x[0]);
      break;
    case 2:
      arglist = Py_BuildValue("(dd)", x[0], x[1]);
      break;
    case 3:
      arglist = Py_BuildValue("(ddd)", x[0], x[1], x[2]);
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
                         "point dimension = " << x << " not supported");
    }
  result = PyEval_CallObject(evalOpCallback_->getFunction(), arglist);
  Py_DECREF(arglist);  // All done with argument list

  if (0 == result) {
    PyErr_Print();
    return false;
  }
  Py_DECREF(result); // All done with returned result object
  return (bool) PyObject_IsTrue(result);
}

string PySundanceCellPredicate::description() const 
{
  PyObject* result = 0 ;

  if (descrCallback_.get() != 0)
    {
      PyObject* arglist;
      
      arglist = Py_BuildValue("()");
      result = PyEval_CallObject(descrCallback_->getFunction(), arglist);
      Py_DECREF(arglist);  // All done with argument list
    }
  else
    {
      result = PyObject_Str(py_functor_);
    }
      
  if (0 == result) 
    {
      PyErr_Print();
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "zero result from python callback");
    }
  
  Py_DECREF(result); // All done with returned result object
  
  char* str = 0;
  Py_ssize_t len = 0;
  PyString_AsStringAndSize(result, &str, &len);

  return std::string(str);
}


PyObject * PySundanceCellPredicate::setEvalOp(PyObject* pyClass)
{
  PySundanceCallback* cb = new PySundanceCallback();
  evalOpCallback_ = rcp(cb);
  return evalOpCallback_->setFunction(pyClass);
}


PyObject * PySundanceCellPredicate::setDescr(PyObject* pyClass)
{
  PySundanceCallback* cb = new PySundanceCallback();
  descrCallback_ = rcp(cb);
  return descrCallback_->setFunction(pyClass);
}

