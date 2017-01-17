// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//          PyTrilinos: Python Interfaces to Trilinos Packages
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
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
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

%define %nox_solver_docstring
"
PyTrilinos.NOX.Solver is the python interface to the Solver namespace
of the Trilinos package NOX:

    http://trilinos.sandia.gov/packages/nox

The purpose of NOX.Solver is to provide solver manager classes for
NOX.  NOX.Solver provides the following user-level classes:

    * Generic                  - Base class for solver managers
    * LineSearchBased          - Line-search-based solver manager
    * TrustRegionBased         - Trust-region-based solver manager
    * InexactTrustRegionBased  - Inexact-trust-region-based solver
                                 manager
    * TensorBased              - Tensor-based solver manager

in addition to the following factory function:

    * buildSolver              - Recommended method for creating solver
                                 managers (note that without loss of
                                 functionality, the Factory class is not
                                 currently provided).
"
%enddef

%module(package      = "PyTrilinos.NOX",
	autodoc      = "1",
	implicitconv = "1",
	docstring    = %nox_solver_docstring) Solver

%{
// PyTrilinos includes
#include "PyTrilinos_config.h"
#include "PyTrilinos_LinearProblem.hpp"

// Teuchos includes
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif
#include "PyTrilinos_Teuchos_Util.hpp"

// NOX includes
#include "NOX_StatusTest_Generic.H"
#include "NOX_StatusTest_NormWRMS.H"
#include "NOX_StatusTest_Stagnation.H"
#include "NOX_StatusTest_MaxIters.H"
#include "NOX_StatusTest_Combo.H"
#include "NOX_StatusTest_FiniteValue.H"
#include "NOX_StatusTest_NormF.H"
#include "NOX_StatusTest_NormUpdate.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Solver_LineSearchBased.H"
#include "NOX_Solver_TrustRegionBased.H"
#include "NOX_Solver_InexactTrustRegionBased.H"
#include "NOX_Solver_TensorBased.H"
#include "NOX_Solver_Factory.H"

// Local includes
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"
%}

// Configuration and optional includes
%include "PyTrilinos_config.h"
#ifdef HAVE_NOX_EPETRA
%{
#include "Epetra_SerialComm.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_SerialDistributor.h"
#include "Epetra_OffsetIndex.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Export.h"
#include "Epetra_IntVector.h"
#include "Epetra_Vector.h"
#include "Epetra_FEVector.h"
#include "Epetra_SerialDenseSVD.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialSymDenseMatrix.h"
#include "Epetra_InvOperator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVbrMatrix.h"
#include "Epetra_JadMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_MapColoring.h"
#include "NOX_Epetra_Group.H"
#include "NOX_Epetra_Vector.H"
#include "NOX_Epetra_ModelEvaluatorInterface.H"
#include "NOX_Epetra_MatrixFree.H"
#include "NOX_Epetra_FiniteDifference.H"
#include "NOX_Epetra_FiniteDifferenceColoring.H"
#include "NOX_Epetra_LinearSystem_AztecOO.H"
#include "NOX_Epetra_Interface_Preconditioner.H"
%}
#endif

// Standard exception handling
%include "exception.i"

// Include NOX documentation
%include "NOX_dox.i"

// General ignore directives
%ignore operator<<(ostream &, NOX::StatusTest::StatusType );
%ignore *::print(ostream& stream, int indent = 0) const;

// Rename directives
%rename(StatusTest_Generic) NOX::StatusTest::Generic;
%rename(StatusTest_None   ) NOX::StatusTest::None;

// Trilinos imports
%import  "Teuchos.i"
%import  "NOX.Abstract.i"
%import  "NOX.StatusTest.i"

// General exception handling
%feature("director:except")
{
  if ($error != NULL)
  {
    throw Swig::DirectorMethodException();
  }
}

%exception
{
  try
  {
    $action
    if (PyErr_Occurred()) SWIG_fail;
  }
  catch(Teuchos::Exceptions::InvalidParameterType & e)
  {
    SWIG_exception(SWIG_TypeError, e.what());
  }
  catch(Teuchos::Exceptions::InvalidParameter & e)
  {
    PyErr_SetString(PyExc_KeyError, e.what());
    SWIG_fail;
  }
  catch(Swig::DirectorException &e)
  {
    SWIG_fail;
  }
  catch(const char* e)
  {
    PyErr_SetString(PyExc_RuntimeError, e);
    SWIG_fail;
  }
  SWIG_CATCH_STDEXCEPT
  catch(...)
  {
    SWIG_exception(SWIG_UnknownError, "Unknown C++ exception");
  }
}

//////////////////////////////////
// NOX::Solver::Generic support //
//////////////////////////////////
%ignore *::getSolutionGroup;
%ignore *::getPreviousSolutionGroup;
%ignore *::getList;
%rename(getSolutionGroup           ) *::getSolutionGroupPtr;
%rename(getPreviousSolutionGroupPtr) *::getPreviousSolutionGroupPtr;
%rename(getList                    ) *::getListPtr;
%teuchos_rcp(NOX::Solver::Generic)
%include "NOX_Solver_Generic.H"

//////////////////////////////////////////
// NOX::Solver::LineSearchBased support //
//////////////////////////////////////////
%teuchos_rcp(NOX::Solver::LineSearchBased)
%include "NOX_Solver_LineSearchBased.H"

///////////////////////////////////////////
// NOX::Solver::TrustRegionBased support //
///////////////////////////////////////////
%teuchos_rcp(NOX::Solver::TrustRegionBased)
%include "NOX_Solver_TrustRegionBased.H"

//////////////////////////////////////////////////
// NOX::Solver::InexactTrustRegionBased support //
//////////////////////////////////////////////////
%teuchos_rcp(NOX::Solver::InexactTrustRegionBased)
%include "NOX_Solver_InexactTrustRegionBased.H"

//////////////////////////////////////
// NOX::Solver::TensorBased support //
//////////////////////////////////////
%teuchos_rcp(NOX::Solver::TensorBased)
%include "NOX_Solver_TensorBased.H"

//////////////////////////////////////
// NOX::Solver::buildSolver support //
//////////////////////////////////////
%rename (buildSolver) myBuildSolver;
// NOX::Solver::buildSolver in NOX_Solver_Factory.H returns a
// Teuchos::RCP<NOX::Solver::Generic>.  As far as I can tell, SWIG
// cannot properly downcast the NOX::Solver::Generic object wrapped
// within the Teuchos::RCP<> in order to, say, call its solve()
// method.  Therefore, I write my own wrapper around buildSolver()
// that does this downcasting explicitly and returns a python wrapper
// around the appropriate derived class.
%inline
{
  PyObject *
    myBuildSolver(const Teuchos::RCP< NOX::Abstract::Group     > & grp,
		  const Teuchos::RCP< NOX::StatusTest::Generic > & tests,
		  const Teuchos::RCP< Teuchos::ParameterList   > & params)
  {
    // SWIG type queries
    static swig_type_info * swig_NSLSB_ptr  =
      SWIG_TypeQuery("Teuchos::RCP< NOX::Solver::LineSearchBased > *");
    static swig_type_info * swig_NSTRB_ptr  =
      SWIG_TypeQuery("Teuchos::RCP< NOX::Solver::TrustRegionBased > *");
    static swig_type_info * swig_NSITRB_ptr =
      SWIG_TypeQuery("Teuchos::RCP< NOX::Solver::InexactTrustRegionBased > *");
    static swig_type_info * swig_NSTB_ptr   =
      SWIG_TypeQuery("Teuchos::RCP< NOX::Solver::TensorBased > *");
    // Build a NOX::Solver::Generic object via the buildSolver factory
    Teuchos::RCP< NOX::Solver::Generic > rcp_solver = NOX::Solver::buildSolver(grp, tests, params);
    // Try to downcast to a derived class
    {
      Teuchos::RCP< NOX::Solver::LineSearchBased > result =
	Teuchos::rcp_dynamic_cast< NOX::Solver::LineSearchBased >(rcp_solver);
      if (!result.is_null())
      {
	Teuchos::RCP< NOX::Solver::LineSearchBased > *smartresult =
	  new Teuchos::RCP< NOX::Solver::LineSearchBased >(result);
	return SWIG_NewPointerObj((void*)smartresult, swig_NSLSB_ptr, 1);
      }
    }
    {
      Teuchos::RCP< NOX::Solver::TrustRegionBased > result =
	Teuchos::rcp_dynamic_cast< NOX::Solver::TrustRegionBased >(rcp_solver);
      if (!result.is_null())
      {
	Teuchos::RCP< NOX::Solver::TrustRegionBased > *smartresult =
	  new Teuchos::RCP< NOX::Solver::TrustRegionBased >(result);
	return SWIG_NewPointerObj((void*)smartresult, swig_NSTRB_ptr, 1);
      }
    }
    {
      Teuchos::RCP< NOX::Solver::InexactTrustRegionBased > result =
	Teuchos::rcp_dynamic_cast< NOX::Solver::InexactTrustRegionBased >(rcp_solver);
      if (!result.is_null())
      {
	Teuchos::RCP< NOX::Solver::InexactTrustRegionBased > *smartresult =
	  new Teuchos::RCP< NOX::Solver::InexactTrustRegionBased >(result);
	return SWIG_NewPointerObj((void*)smartresult, swig_NSITRB_ptr, 1);
      }
    }
    {
      Teuchos::RCP< NOX::Solver::TensorBased > result =
	Teuchos::rcp_dynamic_cast< NOX::Solver::TensorBased >(rcp_solver);
      if (!result.is_null())
      {
	Teuchos::RCP< NOX::Solver::TensorBased > *smartresult =
	  new Teuchos::RCP< NOX::Solver::TensorBased >(result);
	return SWIG_NewPointerObj((void*)smartresult, swig_NSTB_ptr, 1);
      }
    }
    PyErr_SetString(PyExc_RuntimeError, "NOX::Solver::buildSolver returned unrecognized "
		    "derivative of NOX::Solver::Generic");
    return NULL;
  }
}

// Turn off the exception handling
%exception;
