/*! @HEADER */
/*
************************************************************************

                CTrilinos:  C interface to Trilinos
                Copyright (2009) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither the name of the Corporation nor the names of the
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Questions? Contact M. Nicole Lemaster (mnlemas@sandia.gov)

************************************************************************
*/
/*! @HEADER */

#ifdef HAVE_CTRILINOS_EXPERIMENTAL
#ifndef CNOX_INTERFACE_HPP
#define CNOX_INTERFACE_HPP

// Interface to the NLS_PetraGroup to provide for
// residual and matrix fill routines.

// ---------- Standard Includes ----------
#include <iostream>

#include "Epetra_Operator.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "LOCA_Epetra.H"
#include "LOCA_Parameter_Vector.H"

 /*!
   \brief  General NOX/LOCA Epetra Interface class for use in
   C and Fortran codes in Matrix-Free mode.

   The residual and (optional) preconditioner functions call back 
   using function pointers supplied by the user. This class is used
   in conjunction with the functions in CNOX_Driver.cpp.
 */

class  CNOX_Interface :
  public LOCA::Epetra::Interface::Required,
  public NOX::Epetra::Interface::Preconditioner,
  public Epetra_Operator
{
public:
  CNOX_Interface(int* nelems, double* statevector,
                 const LOCA::ParameterVector& pVector_,
                 const Epetra_Comm& comm_,
                 void* blackbox_res, void* blackbox_prec,
                 void (*residualFunction)(double *, double *, int, void *),
                 void (*precFunction)(double *, double *, int, double*, void *));
  ~CNOX_Interface();

  //! Compute and return F
  bool computeF(const Epetra_Vector& x, Epetra_Vector& F, FillType flag);

  //! Set a parameter in the user's code.
  void setParameters(const LOCA::ParameterVector& params);

  //! Print solution to output file
  virtual void printSolution(const Epetra_Vector& x, double conParam);

  //! Application Operator: Object that points to the user's evaluation routines.
  /*! This is used to point to the actual routines and to store
   *  auxiliary data required by the user's application for function/Jacobian
   *  evaluations that NOX does not need to know about.  This is type of
   *  passdown class design by the application code.
   */

  Teuchos::RCP<Epetra_Vector> getVector() const;

  // 1 Method for inheritance from NOX::Epetra::Interface::Preconditioner
  // Compute preconditioner \f$M\f$.
  virtual bool computePreconditioner(const Epetra_Vector& x,
                                     Epetra_Operator& Prec,
                                     Teuchos::ParameterList* p = 0);



  // 10 Methods for inheritance from Epetra_Operator
  // Only ApplyInverse is non-trivial -- first 9 satisfied here in header

    int SetUseTranspose(bool UseTranspose)
        { cerr<<"ERROR: No noxlocainterface::SetUseTranspose"<<endl; return -1;};
    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
        { cerr<<"ERROR: No noxlocainterface::Apply"<<endl; return -1;};
    double NormInf() const
        { cerr<<"ERROR: No noxlocainterface::Apply"<<endl; return 1.0;};
    const char* Label() const { return "noxlocainterface::user preconditioner";};
    bool UseTranspose() const { return false;};
    bool HasNormInf() const { return false;};
    const Epetra_Comm& Comm() const {return comm;};
    const Epetra_Map& OperatorDomainMap() const {return *globalMap;};
    const Epetra_Map& OperatorRangeMap() const {return *globalMap;};

    void resetBlackbox(void* blackbox_res_,  void* blackbox_prec_) {
       blackbox_res=blackbox_res_; blackbox_prec=blackbox_prec_; }

    //! Apply the preconditioner
    int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;


  private:

    int N;
    const Epetra_Comm& comm;
    Teuchos::RCP<Epetra_Vector> solution;
    //Teuchos::RCP<Epetra_Vector> global_solution;
    //Teuchos::RCP<Epetra_Import> global_importer;
    Teuchos::RCP<Epetra_Map> globalMap;
    LOCA::ParameterVector pVector;
    void* blackbox_res;
    void* blackbox_prec;
    void (*residualFunction)(double *, double *, int, void *);
    void (*precFunction)(double *, double *, int, double*, void *);

};

#endif

#endif //HAVE_CTRILINOS_EXPERIMENTAL
