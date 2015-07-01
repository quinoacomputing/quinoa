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


#include "PyTeuchos_Utils.hpp"
#include "Teuchos_RefCountPtr.hpp"

// Creates a newly allocated Teuchos parameter list from the input
// object, which must be a Python dictionary.
//
// "bool", "int", "double" and "string" are automatically recognized.
// Other types can be defined here as tuples. For example, the Python
// dictionary can be something like:
// List = {
//   "double parameter": 12.0,
//   "int parameter"   : 12,
//   "string parameter": "12"
// }
//
// \author Marzio Sala, SNL 9215
//
// \date Last modified on 08-Aug-05

Teuchos::ParameterList dict2ParameterList(PyObject* obj)
{
  int i;
  Teuchos::ParameterList List;
  if (!PyDict_Check(obj)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a dictionary");
    return List;
  }

  int size = PyDict_Size(obj);
  PyObject* Keys = PyDict_Keys(obj);
  PyObject* Values = PyDict_Values(obj);

  for (i = 0; i < size ; i++)
    {
      PyObject *s = PyList_GetItem(Keys,i);
      PyObject *t = PyList_GetItem(Values,i);

      // Get the parameter name
      if (!PyString_Check(s)) {
        PyErr_SetString(PyExc_ValueError, "Dictionary keys must be std::strings");
        return List;
      }
      std::string ParameterName = PyString_AsString(s);

      // now parse for the parameter value and type
      // This can be a "int", "double", "string", or a tuple
      // for more general types

      if (PyBool_Check(t))
        {
          if (t == Py_True)
            List.set(ParameterName, true);
          else
            List.set(ParameterName, false);
        }
      else if (PyInt_Check(t))
        {
          int ParameterValue = PyInt_AsLong(t);
          List.set(ParameterName, ParameterValue);
        }
      else if (PyFloat_Check(t))
        {
          double ParameterValue = PyFloat_AsDouble(t);
          List.set(ParameterName, ParameterValue);
        }
      else if (PyString_Check(t))
        {
          std::string ParameterValue = PyString_AsString(t);
          List.set(ParameterName, ParameterValue);
        }
      else if (PyTuple_Check(t))
        {
          if (!PyString_Check(PyTuple_GetItem(t, 0)) ||
              !PyString_Check(PyTuple_GetItem(t, 1))) {
            PyErr_SetString(PyExc_ValueError, "tuples must contain std::strings");
            return List;
          }
          std::string ParameterType = PyString_AsString(PyTuple_GetItem(t, 0));
          std::string ParameterValue = PyString_AsString(PyTuple_GetItem(t, 1));
          if (ParameterType == "bool")
            {
              if (ParameterValue == "true")
                List.set(ParameterName, true);
              else
                List.set(ParameterName, false);
            }
          else if (ParameterType == "int")
            {
              List.set(ParameterName, (int)atoi(ParameterValue.c_str()));
            }
          else if (ParameterType == "double")
            {
              List.set(ParameterName, (double)atof(ParameterValue.c_str()));
            }
          else if (ParameterType == "string")
            {
              List.set(ParameterName, std::string(ParameterValue));
            }
          else
            {
              PyErr_SetString(PyExc_ValueError, "type in tuple not recognized");
              return List;
            }
        }
      else if (PyDict_Check(t))
        {
          Teuchos::ParameterList sublist = dict2ParameterList(t);
          List.set(ParameterName, sublist);
        }
      else
        {
          PyErr_SetString(PyExc_ValueError, "Type in list not recognized");
          return List;
        }
    }

  return(List);
}

