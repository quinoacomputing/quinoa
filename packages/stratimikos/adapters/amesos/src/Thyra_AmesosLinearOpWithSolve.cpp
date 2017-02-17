/*
// @HEADER
// ***********************************************************************
//
//         Stratimikos: Thyra-based strategies for linear solvers
//                Copyright (2006) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov)
//
// ***********************************************************************
// @HEADER
*/

#include "Thyra_AmesosLinearOpWithSolve.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Epetra_MultiVector.h"
#include "Teuchos_TimeMonitor.hpp"


namespace Thyra {


// Constructors/initializers/accessors


AmesosLinearOpWithSolve::AmesosLinearOpWithSolve():
  amesosSolverTransp_(Thyra::NOTRANS),
  amesosSolverScalar_(1.0)
{}


AmesosLinearOpWithSolve::AmesosLinearOpWithSolve(
  const Teuchos::RCP<const LinearOpBase<double> > &fwdOp,
  const Teuchos::RCP<const LinearOpSourceBase<double> > &fwdOpSrc,
  const Teuchos::RCP<Epetra_LinearProblem> &epetraLP,
  const Teuchos::RCP<Amesos_BaseSolver> &amesosSolver,
  const EOpTransp amesosSolverTransp,
  const double amesosSolverScalar
  )
{
  this->initialize(fwdOp,fwdOpSrc,epetraLP,amesosSolver,
    amesosSolverTransp,amesosSolverScalar);
}


void AmesosLinearOpWithSolve::initialize(
  const Teuchos::RCP<const LinearOpBase<double> > &fwdOp,
  const Teuchos::RCP<const LinearOpSourceBase<double> > &fwdOpSrc,
  const Teuchos::RCP<Epetra_LinearProblem> &epetraLP,
  const Teuchos::RCP<Amesos_BaseSolver> &amesosSolver,
  const EOpTransp amesosSolverTransp,
  const double amesosSolverScalar
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(fwdOp.get()==NULL);
  TEUCHOS_TEST_FOR_EXCEPT(fwdOpSrc.get()==NULL);
  TEUCHOS_TEST_FOR_EXCEPT(epetraLP.get()==NULL);
  TEUCHOS_TEST_FOR_EXCEPT(amesosSolver.get()==NULL);
  TEUCHOS_TEST_FOR_EXCEPT(epetraLP->GetLHS()!=NULL);
  TEUCHOS_TEST_FOR_EXCEPT(epetraLP->GetRHS()!=NULL);
#endif
  fwdOp_ = fwdOp;
  fwdOpSrc_ = fwdOpSrc;
  epetraLP_ = epetraLP;
  amesosSolver_ = amesosSolver;
  amesosSolverTransp_ = amesosSolverTransp;
  amesosSolverScalar_ = amesosSolverScalar;
  const std::string fwdOpLabel = fwdOp_->getObjectLabel();
  if(fwdOpLabel.length())
    this->setObjectLabel( "lows("+fwdOpLabel+")" );
}


Teuchos::RCP<const LinearOpSourceBase<double> >
AmesosLinearOpWithSolve::extract_fwdOpSrc()
{
  Teuchos::RCP<const LinearOpSourceBase<double> >
    _fwdOpSrc = fwdOpSrc_;
  fwdOpSrc_ = Teuchos::null;
  return _fwdOpSrc;
}


void AmesosLinearOpWithSolve::uninitialize(
  Teuchos::RCP<const LinearOpBase<double> > *fwdOp,
  Teuchos::RCP<const LinearOpSourceBase<double> > *fwdOpSrc,
  Teuchos::RCP<Epetra_LinearProblem> *epetraLP,
  Teuchos::RCP<Amesos_BaseSolver> *amesosSolver,
  EOpTransp *amesosSolverTransp,
  double *amesosSolverScalar
  )
{

  if(fwdOp) *fwdOp = fwdOp_;
  if(fwdOpSrc) *fwdOpSrc = fwdOpSrc_;
  if(epetraLP) *epetraLP = epetraLP_;
  if(amesosSolver) *amesosSolver = amesosSolver_;
  if(amesosSolverTransp) *amesosSolverTransp = amesosSolverTransp_;
  if(amesosSolverScalar) *amesosSolverScalar = amesosSolverScalar_;

  fwdOp_ = Teuchos::null;
  fwdOpSrc_ = Teuchos::null;
  epetraLP_ = Teuchos::null;
  amesosSolver_ = Teuchos::null;
  amesosSolverTransp_ = NOTRANS;
  amesosSolverScalar_ = 0.0;

}


// Overridden from LinearOpBase


Teuchos::RCP< const VectorSpaceBase<double> >
AmesosLinearOpWithSolve::range() const
{
  return ( fwdOp_.get() ? fwdOp_->range() : Teuchos::null );
}


Teuchos::RCP< const VectorSpaceBase<double> >
AmesosLinearOpWithSolve::domain() const
{
  return  ( fwdOp_.get() ? fwdOp_->domain() : Teuchos::null );
}


Teuchos::RCP<const LinearOpBase<double> >
AmesosLinearOpWithSolve::clone() const
{
  return Teuchos::null; // Not supported yet but could be
}


// Overridden from Teuchos::Describable


std::string AmesosLinearOpWithSolve::description() const
{
  std::ostringstream oss;
  oss << Teuchos::Describable::description();
  if(!is_null(amesosSolver_)) {
    oss
      << "{fwdOp="<<fwdOp_->description()
      << ",amesosSolver="<<typeName(*amesosSolver_)<<"}";
  }
  return oss.str();
}


void AmesosLinearOpWithSolve::describe(
  Teuchos::FancyOStream &out,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  using Teuchos::OSTab;
  using Teuchos::typeName;
  using Teuchos::describe;
  switch(verbLevel) {
    case Teuchos::VERB_DEFAULT:
    case Teuchos::VERB_LOW:
      out << this->description() << std::endl;
      break;
    case Teuchos::VERB_MEDIUM:
    case Teuchos::VERB_HIGH:
    case Teuchos::VERB_EXTREME:
    {
      out
        << Teuchos::Describable::description() << "{"
        << "rangeDim=" << this->range()->dim()
        << ",domainDim="<< this->domain()->dim() << "}\n";
      OSTab tab(out);
      if(!is_null(fwdOp_)) {
        out << "fwdOp = " << describe(*fwdOp_,verbLevel);
      }
      if(!is_null(amesosSolver_)) {
        out << "amesosSolver=" << typeName(*amesosSolver_) << "\n";
      }
      break;
    }
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true); // Should never get here!
  }
}


// protected


// Overridden from LinearOpBase


bool AmesosLinearOpWithSolve::opSupportedImpl(EOpTransp M_trans) const
{
  return ::Thyra::opSupported(*fwdOp_,M_trans);
}


void AmesosLinearOpWithSolve::applyImpl(
  const EOpTransp M_trans,
  const MultiVectorBase<double> &X,
  const Ptr<MultiVectorBase<double> > &Y,
  const double alpha,
  const double beta
  ) const
{
  Thyra::apply( *fwdOp_, M_trans, X, Y, alpha, beta );
}


// Overridden from LinearOpWithSolveBase


bool AmesosLinearOpWithSolve::solveSupportsImpl(EOpTransp M_trans) const
{
  if (Thyra::real_trans(M_trans) == Thyra::NOTRANS) {
    // Assume every amesos solver supports a basic forward solve!
    return true;
  }
  // Query the amesos solver to see if it supports the transpose operation.
  // NOTE: Amesos_BaseSolver makes you change the state of the object to
  // determine if the object supports an adjoint solver.  This is a bad design
  // but I have no control over that.  This is why you see this hacked
  // oldUseTranspose variable and logic.  NOTE: This function meets the basic
  // guarantee but if setUseTransplse(...) throws, then the state of
  // UseTranspose() may be different.
  const bool oldUseTranspose = amesosSolver_->UseTranspose();
  const bool supportsAdjoint = (amesosSolver_->SetUseTranspose(true) == 0);
  amesosSolver_->SetUseTranspose(oldUseTranspose);
  return supportsAdjoint;
}


bool AmesosLinearOpWithSolve::solveSupportsSolveMeasureTypeImpl(
  EOpTransp M_trans, const SolveMeasureType& solveMeasureType
  ) const
{
  return true; // I am a direct solver so I should be able to do it all!
}


SolveStatus<double>
AmesosLinearOpWithSolve::solveImpl(
  const EOpTransp M_trans,
  const MultiVectorBase<double> &B,
  const Ptr<MultiVectorBase<double> > &X,
  const Ptr<const SolveCriteria<double> > solveCriteria
  ) const
{
  using Teuchos::rcpFromPtr;
  using Teuchos::rcpFromRef;
  using Teuchos::OSTab;

  Teuchos::Time totalTimer("");
  totalTimer.start(true);

  THYRA_FUNC_TIME_MONITOR("Stratimikos: AmesosLOWS");

  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  OSTab tab = this->getOSTab();
  if(out.get() && static_cast<int>(verbLevel) > static_cast<int>(Teuchos::VERB_NONE))
    *out << "\nSolving block system using Amesos solver "
         << typeName(*amesosSolver_) << " ...\n\n";

  //
  // Get the op(...) range and domain maps
  //
  const EOpTransp amesosOpTransp = real_trans(trans_trans(amesosSolverTransp_,M_trans));
  const Epetra_Operator *amesosOp = epetraLP_->GetOperator();
  const Epetra_Map
    &opRangeMap  = ( amesosOpTransp == NOTRANS
      ? amesosOp->OperatorRangeMap()  : amesosOp->OperatorDomainMap() ),
    &opDomainMap = ( amesosOpTransp == NOTRANS
      ? amesosOp->OperatorDomainMap() : amesosOp->OperatorRangeMap()  );

  //
  // Get Epetra_MultiVector views of B and X
  //
  Teuchos::RCP<const Epetra_MultiVector>
    epetra_B = get_Epetra_MultiVector(opRangeMap, rcpFromRef(B));
  Teuchos::RCP<Epetra_MultiVector>
    epetra_X = get_Epetra_MultiVector(opDomainMap, rcpFromPtr(X));

  //
  // Set B and X in the linear problem
  //
  epetraLP_->SetLHS(&*epetra_X);
  epetraLP_->SetRHS(const_cast<Epetra_MultiVector*>(&*epetra_B));
  // Above should be okay but cross your fingers!

  //
  // Solve the linear system
  //
  const bool oldUseTranspose = amesosSolver_->UseTranspose();
  amesosSolver_->SetUseTranspose(amesosOpTransp==TRANS);
  const int err = amesosSolver_->Solve();
  TEUCHOS_TEST_FOR_EXCEPTION( 0!=err, CatastrophicSolveFailure,
    "Error, the function Solve() on the amesos solver of type\n"
    "\'"<<typeName(*amesosSolver_)<<"\' failed with error code "<<err<<"!"
    );
  amesosSolver_->SetUseTranspose(oldUseTranspose);

  //
  // Unset B and X
  //
  epetraLP_->SetLHS(NULL);
  epetraLP_->SetRHS(NULL);
  epetra_X = Teuchos::null;
  epetra_B = Teuchos::null;

  //
  // Scale X if needed
  //
  if(amesosSolverScalar_!=1.0)
    Thyra::scale(1.0/amesosSolverScalar_, X);

  //
  // Set the solve status if requested
  //
  SolveStatus<double> solveStatus;
  solveStatus.solveStatus = SOLVE_STATUS_CONVERGED;
  solveStatus.achievedTol = SolveStatus<double>::unknownTolerance();
  solveStatus.message =
    std::string("Solver ")+typeName(*amesosSolver_)+std::string(" converged!");

  //
  // Report the overall time
  //
  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out
      << "\nTotal solve time = "<<totalTimer.totalElapsedTime()<<" sec\n";

  return solveStatus;

}


}	// end namespace Thyra
