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


#include "PySundanceCallback.hpp"

#include <iostream>

using std::cout;
using std::endl;

PySundanceCallback::PySundanceCallback()
  : callback_(NULL)
{
  // Nothing here yet
}

PySundanceCallback::~PySundanceCallback()
{
  Py_XDECREF(callback_);   /* Dispose of callback */
  callback_ = 0;
}

PyObject * PySundanceCallback::setFunction( PyObject * pyMethod)
{
    PyObject *p_result = NULL;
    PyObject *p_temp   = NULL;
    PyObject *p_func   = NULL;


    assert(0 != pyMethod && "Null argument passed to setFunction()");
    
    if (PyArg_ParseTuple(pyMethod, "O", &p_temp)) {
      // Assume that if this is a tuple, the item in the tuple is a
      // PyObject that is a pointer to a Python function. This is the
      // case when this function is called from Python.
      p_func = p_temp;
    } else {
      // Otherwise we assume that this function is directly passed a
      // PyObject that is a pointer to a Python function.  This is the
      // case when this function is called from C++.
      p_func = pyMethod;
    }

    if (!PyCallable_Check(p_func)) {
      PyErr_SetString(PyExc_TypeError,
		      "Function parameter must be callable");
      cout << "PyObject passed to function is not callable" << std::endl ;
      return NULL;
    }
    Py_XINCREF(p_func);          /* Add a reference to new callback */
    Py_XDECREF(callback_);     /* Dispose of previous callback    */
    callback_ = p_func;        /* Remember new callback           */
    /* Boilerplate to return "None" */
    Py_INCREF(Py_None);
    p_result = Py_None;

    assert(0 != callback_ && "Pointer to callback not set");
    return p_result;
}
 
PyObject * PySundanceCallback::getFunction()
{
  assert (0 != callback_ && "PySundanceCallback function not yet assigned");
  return callback_;
}

const PyObject * PySundanceCallback::getFunction() const
{
  assert (0 != callback_ && "PySundanceCallback function not yet assigned");
  return callback_;
}
