/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_CHEBYSHEV_DEF_HPP
#define IFPACK2_CHEBYSHEV_DEF_HPP


namespace Ifpack2 {

//Definitions for the Chebyshev methods:

//==========================================================================
template<class MatrixType>
Chebyshev<MatrixType>::Chebyshev(const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A)
: A_(A),
  Comm_(A->getRowMap()->getComm()),
  Time_( Teuchos::rcp( new Teuchos::Time("Ifpack2::Chebyshev") ) ),
  PolyDegree_(1),
  EigRatio_(30.0),
  LambdaMin_(0.0),
  LambdaMax_(100.0),
  MinDiagonalValue_(0.0),
  ZeroStartingSolution_(true),
  Condest_(-1.0),
  IsInitialized_(false),
  IsComputed_(false),
  NumInitialize_(0),
  NumCompute_(0),
  NumApply_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyTime_(0.0),
  ComputeFlops_(0.0),
  ApplyFlops_(0.0),
  NumMyRows_(0),
  NumGlobalRows_(0),
  NumGlobalNonzeros_(0)
{
  TEUCHOS_TEST_FOR_EXCEPTION(A_ == Teuchos::null, std::runtime_error, 
      Teuchos::typeName(*this) << "::Chebyshev(): input matrix reference was null.");
}

//==========================================================================
template<class MatrixType>
Chebyshev<MatrixType>::~Chebyshev() {
}

//==========================================================================
template<class MatrixType>
void Chebyshev<MatrixType>::setParameters(const Teuchos::ParameterList& List) {

  Ifpack2::getParameter(List, "chebyshev: ratio eigenvalue", EigRatio_);
  Ifpack2::getParameter(List, "chebyshev: min eigenvalue", LambdaMin_);
  Ifpack2::getParameter(List, "chebyshev: max eigenvalue", LambdaMax_);
  Ifpack2::getParameter(List, "chebyshev: degree",PolyDegree_);
  Ifpack2::getParameter(List, "chebyshev: min diagonal value", MinDiagonalValue_);
  Ifpack2::getParameter(List, "chebyshev: zero starting solution", ZeroStartingSolution_);

  Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>* ID = 0;
  Ifpack2::getParameter(List, "chebyshev: operator inv diagonal", ID);

  if (ID != 0) {
    InvDiagonal_ = Teuchos::rcp( new Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(*ID) );
  }
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Teuchos::Comm<int> > & 
Chebyshev<MatrixType>::getComm() const{
  return(Comm_);
}

//==========================================================================
template<class MatrixType>
Teuchos::RCP<const Tpetra::RowMatrix<typename MatrixType::scalar_type,typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> >
Chebyshev<MatrixType>::getMatrix() const {
  return(A_);
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> >&
Chebyshev<MatrixType>::getDomainMap() const {
  return A_->getDomainMap();
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> >&
Chebyshev<MatrixType>::getRangeMap() const {
  return A_->getRangeMap();
}

//==========================================================================
template<class MatrixType>
bool Chebyshev<MatrixType>::hasTransposeApply() const {
  return true;
}

//==========================================================================
template<class MatrixType>
int Chebyshev<MatrixType>::getNumInitialize() const {
  return(NumInitialize_);
}

//==========================================================================
template<class MatrixType>
int Chebyshev<MatrixType>::getNumCompute() const {
  return(NumCompute_);
}

//==========================================================================
template<class MatrixType>
int Chebyshev<MatrixType>::getNumApply() const {
  return(NumApply_);
}

//==========================================================================
template<class MatrixType>
double Chebyshev<MatrixType>::getInitializeTime() const {
  return(InitializeTime_);
}

//==========================================================================
template<class MatrixType>
double Chebyshev<MatrixType>::getComputeTime() const {
  return(ComputeTime_);
}

//==========================================================================
template<class MatrixType>
double Chebyshev<MatrixType>::getApplyTime() const {
  return(ApplyTime_);
}

//==========================================================================
template<class MatrixType>
double Chebyshev<MatrixType>::getComputeFlops() const {
  return(ComputeFlops_);
}

//==========================================================================
template<class MatrixType>
double Chebyshev<MatrixType>::getApplyFlops() const {
  return(ApplyFlops_);
}

//==========================================================================
template<class MatrixType>
typename Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType
Chebyshev<MatrixType>::getCondEst() const {
  return(Condest_);
}

//==========================================================================
template<class MatrixType>
typename Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType
Chebyshev<MatrixType>::computeCondEst(
                     CondestType CT,
                     LocalOrdinal MaxIters, 
                     magnitudeType Tol,
                     const Teuchos::Ptr<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &matrix) {
  if (!isComputed()) // cannot compute right now
    return(-1.0);

  // always compute it. Call Condest() with no parameters to get
  // the previous estimate.
  Condest_ = Ifpack2::Condest(*this, CT, MaxIters, Tol, matrix);

  return(Condest_);
}

//==========================================================================
template<class MatrixType>
void Chebyshev<MatrixType>::apply(
           const Tpetra::MultiVector<typename MatrixType::scalar_type, typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type>& X,
                 Tpetra::MultiVector<typename MatrixType::scalar_type, typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type>& Y,
                Teuchos::ETransp mode,
                 Scalar alpha,
                 Scalar beta) const {
  TEUCHOS_TEST_FOR_EXCEPTION(!isComputed(), std::runtime_error, 
      "Ifpack2::Chebyshev::apply() ERROR, not yet computed.");

  TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
     "Ifpack2::Chebyshev::apply() ERROR: X.getNumVectors() != Y.getNumVectors().");

  Time_->start();

  // If X and Y are pointing to the same memory location,
  // we need to create an auxiliary vector, Xcopy
  Teuchos::RCP< const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Xcopy;
  if (X.getLocalMV().getValues() == Y.getLocalMV().getValues())
    Xcopy = Teuchos::rcp( new Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(X) );
  else
    Xcopy = Teuchos::rcp( &X, false );

  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > xView = Xcopy->get2dView();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > yView = Y.get2dViewNonConst();
  Teuchos::ArrayRCP<const Scalar> invdiag = InvDiagonal_->get1dView();

  size_t nVecs = Y.getNumVectors();

  //--- Do a quick solve when the matrix is identity
  if ((LambdaMin_ == 1.0) && (LambdaMax_ == LambdaMin_)) {
    if (nVecs == 1) {
      Teuchos::ArrayRCP<Scalar> y = yView[0];
      Teuchos::ArrayRCP<const Scalar> x = xView[0];
      for (size_t i = 0; i < NumMyRows_; ++i)
        y[i] = x[i]*invdiag[i];
    }
    else {
      for (size_t i = 0; i < NumMyRows_; ++i) {
        const Scalar& coeff = invdiag[i];
        for (size_t k = 0; k < nVecs; ++k)
          yView[k][i] = xView[k][i] * coeff;
      }
    } // if (nVec == 1)
    return;
  } // if ((LambdaMin_ == 1.0) && (LambdaMax_ == LambdaMin_))

  //--- initialize coefficients
  // Note that delta stores the inverse of ML_Cheby::delta
  Scalar alpha_cheby = LambdaMax_ / EigRatio_;
  Scalar beta_cheby = 1.1 * LambdaMax_;
  Scalar delta = 2.0 / (beta_cheby - alpha_cheby);
  Scalar theta = 0.5 * (beta_cheby + alpha_cheby);
  Scalar s1 = theta * delta;

  //--- Define vectors
  // In ML_Cheby, V corresponds to pAux and W to dk
  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> V(X);
  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> W(X);

  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > vView = V.get2dView();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > wView = W.get2dViewNonConst();

  Scalar one = Teuchos::ScalarTraits<Scalar>::one();

  Scalar oneOverTheta = one/theta;

  // Do the smoothing when block scaling is turned OFF
  // --- Treat the initial guess
  if (ZeroStartingSolution_ == false) {
    A_->apply(Y, V);
    // compute W = invDiag * ( X - V )/ Theta
    if (nVecs == 1) {
      Teuchos::ArrayRCP<const Scalar> x = xView[0];
      Teuchos::ArrayRCP<Scalar> w = wView[0];
      Teuchos::ArrayRCP<const Scalar> v = vView[0];
      for (size_t i = 0; i < NumMyRows_; ++i)
        w[i] = invdiag[i] * (x[i] - v[i]) * oneOverTheta;
    }
    else {
      for (size_t k = 0; k < nVecs; ++k) {
        Teuchos::ArrayRCP<Scalar> wk = wView[k];
        Teuchos::ArrayRCP<const Scalar> vk = vView[k];
        for (size_t i = 0; i < NumMyRows_; ++i) {
          Scalar coeff = invdiag[i]*oneOverTheta;
          wk[i] = (xView[k][i] - (vk[i])) * coeff;
        }
      }
    } // if (nVec == 1)
    // Update the vector Y
    Y.update(one, W, one);
  }
  else {
    // compute W = invDiag * X / Theta
    if (nVecs == 1) {
      Teuchos::ArrayRCP<const Scalar> x= xView[0];
      Teuchos::ArrayRCP<Scalar> w = wView[0];
      Teuchos::ArrayRCP<Scalar> y = yView[0];
      for (size_t i = 0; i < NumMyRows_; ++i) {
        w[i] = invdiag[i] * x[i] * oneOverTheta;
        y[i] = w[i];
      }
    }
    else {
      for (size_t k = 0; k < nVecs; ++k) {
        for (size_t i = 0; i < NumMyRows_; ++i) {
          Scalar coeff = invdiag[i]*oneOverTheta;
          wView[k][i] = xView[k][i] * coeff;
          yView[k][i] = wView[k][i];
        }
      }
    } // if (nVec == 1)
  } // if (ZeroStartingSolution_ == false)

  //--- Apply the polynomial
  Scalar rhok = 1.0/s1, rhokp1;
  Scalar dtemp1, dtemp2;
  int degreeMinusOne = PolyDegree_ - 1;
  Scalar two = 2.0;
  for (int deg = 0; deg < degreeMinusOne; ++deg) {
    A_->apply(Y, V);
    rhokp1 = one / (two *s1 - rhok);
    dtemp1 = rhokp1 * rhok;
    dtemp2 = two * rhokp1 * delta;
    rhok = rhokp1;
    // compute W = dtemp1 * W
    W.scale(dtemp1);
    // compute W = W + dtemp2 * invDiag * ( X - V )
    for (size_t k = 0; k < nVecs; ++k) {
      for (size_t i = 0; i < NumMyRows_; ++i) {
        Scalar coeff = invdiag[i]*dtemp2;
        wView[k][i] += (xView[k][i] - (vView[k][i])) * coeff;
      }
    }
    // Update the vector Y
    Y.update(one, W, one);
  } // for (deg = 0; deg < degreeMinusOne; ++deg)

  // Flops are updated in each of the following. 

  ++NumApply_;
  Time_->stop();
  ApplyTime_ += Time_->totalElapsedTime();
}

//==========================================================================
template<class MatrixType>
void Chebyshev<MatrixType>::applyMat(
           const Tpetra::MultiVector<typename MatrixType::scalar_type, typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type>& X,
                 Tpetra::MultiVector<typename MatrixType::scalar_type, typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type>& Y,
             Teuchos::ETransp mode) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(isComputed() == false, std::runtime_error,
     "Ifpack2::Chebyshev::applyMat() ERROR: isComputed() must be true prior to calling applyMat().");
  TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
     "Ifpack2::Chebyshev::applyMat() ERROR: X.getNumVectors() != Y.getNumVectors().");
  A_->apply(X, Y, mode);
}

//==========================================================================
template<class MatrixType>
void Chebyshev<MatrixType>::initialize() {
  IsInitialized_ = false;

  TEUCHOS_TEST_FOR_EXCEPTION(A_ == Teuchos::null, std::runtime_error, 
      "Ifpack2::Chebyshev::initialize ERROR: A_ == Teuchos::null");

  TEUCHOS_TEST_FOR_EXCEPTION(A_->getGlobalNumRows() != A_->getGlobalNumCols(), std::runtime_error,
     "Ifpack2::Chebyshev::initialize ERROR: only square matrices are supported");

  NumMyRows_ = A_->getNodeNumRows();
  NumGlobalRows_ = A_->getGlobalNumRows();
  NumGlobalNonzeros_ = A_->getGlobalNumEntries();

  ++NumInitialize_;
  Time_->stop();
  InitializeTime_ += Time_->totalElapsedTime();
  IsInitialized_ = true;
}

//==========================================================================
template<class MatrixType>
void Chebyshev<MatrixType>::compute()
{
  if (!isInitialized()) {
    initialize();
  }

  Time_->start(true);

  // reset values
  IsComputed_ = false;
  Condest_ = -1.0;

  TEUCHOS_TEST_FOR_EXCEPTION(PolyDegree_ <= 0, std::runtime_error,
    "Ifpack2::Chebyshev::compute() ERROR: PolyDegree_ must be at least one");
  
  if (InvDiagonal_ == Teuchos::null)
  {
    InvDiagonal_ = Teuchos::rcp( new Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(A_->getRowMap()) );

    A_->getLocalDiagCopy(*InvDiagonal_);

    // Inverse diagonal elements
    // Replace zeros with 1.0
    Teuchos::ArrayRCP<Scalar> diagvals = InvDiagonal_->get1dViewNonConst();
    for (size_t i = 0 ; i < NumMyRows_ ; ++i) {
      Scalar& diag = diagvals[i];
      if (Teuchos::ScalarTraits<Scalar>::magnitude(diag) < Teuchos::ScalarTraits<Scalar>::magnitude(MinDiagonalValue_))
        diag = MinDiagonalValue_;
      else
        diag = 1.0 / diag;
    }
  }
  // otherwise the inverse of the diagonal has been given by the user

  ComputeFlops_ += NumMyRows_;

  ++NumCompute_;
  Time_->stop();
  ComputeTime_ += Time_->totalElapsedTime();
  IsComputed_ = true;
}

//==========================================================================
template<class MatrixType>
void Chebyshev<MatrixType>::
PowerMethod(const Tpetra::Operator<typename MatrixType::scalar_type, typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type>& Operator, 
            const Tpetra::Vector<typename MatrixType::scalar_type, typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type>& InvPointDiagonal, 
            const int MaximumIterations, 
            typename MatrixType::scalar_type& lambda_max)
{
  // this is a simple power method
  lambda_max = 0.0;
  Teuchos::Array<Scalar> RQ_top(1), RQ_bottom(1);
  Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> x(Operator.getDomainMap());
  Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> y(Operator.getRangeMap());
  x.randomize();
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;
  Teuchos::Array<magnitudeType> norms(x.getNumVectors());
  x.norm2(norms());

  x.scale(1.0 / norms[0]);

  Scalar one = Teuchos::ScalarTraits<Scalar>::one();
  Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();

  for (int iter = 0; iter < MaximumIterations; ++iter)
  {
    Operator.apply(x, y);
    y.elementWiseMultiply(1.0, InvPointDiagonal, y, 0.0);
    y.dot(x, RQ_top());
    x.dot(x, RQ_bottom());
    lambda_max = RQ_top[0] / RQ_bottom[0];
    y.norm2(norms());
    TEUCHOS_TEST_FOR_EXCEPTION(norms[0] == zero, std::runtime_error, "Ifpack2::Chebyshev::PowerMethod ERROR, norm == 0");
    x.update( one / norms[0], y, zero);
  }
}

//==========================================================================
template<class MatrixType>
void Chebyshev<MatrixType>::
CG(const Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Operator, 
            const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& InvPointDiagonal, 
   const int MaximumIterations, 
   Scalar& lambda_min, Scalar& lambda_max)
{
#ifdef HAVE_IFPACK2_AZTECOO
  Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> x(Operator.getDomainMap());
  Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> y(Operator.getRangeMap());
  x.Random();
  y.putScalar(0.0);

  Tpetra::LinearProblem LP(const_cast<Tpetra::Operator*>(&Operator), &x, &y);
  AztecOO solver(LP);
  solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  solver.SetAztecOption(AZ_output, AZ_none);

  Ifpack2_DiagPreconditioner diag(Operator.OperatorDomainMap(),
                                 Operator.OperatorRangeMap(),
                                 InvPointDiagonal);
  solver.SetPrecOperator(&diag);
  solver.Iterate(MaximumIterations, 1e-10);

  const double* status = solver.GetAztecStatus();

  lambda_min = status[AZ_lambda_min];
  lambda_max = status[AZ_lambda_max];

  return(0);
#else
  throw std::runtime_error("Ifpack2::Chebyshev::CG: support for AztecOO not currently implemented.");
#endif
}

//==========================================================================
template <class MatrixType>
std::string Chebyshev<MatrixType>::description() const {
  std::ostringstream oss;
  oss << Teuchos::Describable::description();
  if (isInitialized()) {
    if (isComputed()) {
      oss << "{status = initialized, computed";
    }
    else {
      oss << "{status = initialized, not computed";
    }
  }
  else {
    oss << "{status = not initialized, not computed";
  }
  //
  oss << ", global rows = " << A_->getGlobalNumRows()
      << ", global cols = " << A_->getGlobalNumCols()
      << "}";
  return oss.str();
}

//==========================================================================
template <class MatrixType>
void Chebyshev<MatrixType>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
  using std::endl;
  using std::setw;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;
  Teuchos::EVerbosityLevel vl = verbLevel;
  if (vl == VERB_DEFAULT) vl = VERB_LOW;
  const int myImageID = Comm_->getRank();
  Teuchos::OSTab tab(out);

  Scalar MinVal, MaxVal;
  if (IsComputed_) {
    Teuchos::ArrayRCP<const Scalar> DiagView = InvDiagonal_->get1dView();
    Scalar myMinVal = DiagView[0];
    Scalar myMaxVal = DiagView[0];
    for(typename Teuchos::ArrayRCP<Scalar>::size_type i=1; i<DiagView.size(); ++i) {
      if (Teuchos::ScalarTraits<Scalar>::magnitude(myMinVal) > Teuchos::ScalarTraits<Scalar>::magnitude(DiagView[i])) myMinVal = DiagView[i];
      if (Teuchos::ScalarTraits<Scalar>::magnitude(myMaxVal) < Teuchos::ScalarTraits<Scalar>::magnitude(DiagView[i])) myMaxVal = DiagView[i];
    }
    Teuchos::reduceAll(*Comm_, Teuchos::REDUCE_MIN, 1, &myMinVal, &MinVal);
    Teuchos::reduceAll(*Comm_, Teuchos::REDUCE_MAX, 1, &myMaxVal, &MaxVal);
  }

  //    none: print nothing
  //     low: print O(1) info from node 0
  //  medium: 
  //    high: 
  // extreme: 
  if (vl != VERB_NONE && myImageID == 0) {
    out << this->description() << endl;
    out << endl;
    out << "===============================================================================" << std::endl;
    out << "Degree of polynomial      = " << PolyDegree_ << std::endl;
    if   (ZeroStartingSolution_) { out << "Using zero starting solution" << endl; }
    else                         { out << "Using input starting solution" << endl; }
    if   (Condest_ == -1.0) { out << "Condition number estimate       = N/A" << endl; }
    else                    { out << "Condition number estimate       = " << Condest_ << endl; }
    if (IsComputed_) {
      out << "Minimum value on stored inverse diagonal = " << MinVal << std::endl;
      out << "Maximum value on stored inverse diagonal = " << MaxVal << std::endl;
    }
    out << std::endl;
    out << "Phase           # calls    Total Time (s)     Total MFlops      MFlops/s       " << endl;
    out << "------------    -------    ---------------    ---------------   ---------------" << endl;
    out << setw(12) << "initialize()" << setw(5) << getNumInitialize() << "    " << setw(15) << getInitializeTime() << endl;
    out << setw(12) << "compute()" << setw(5) << getNumCompute()    << "    " << setw(15) << getComputeTime() << "    " 
        << setw(15) << getComputeFlops() << "    " 
        << setw(15) << (getComputeTime() != 0.0 ? getComputeFlops() / getComputeTime() * 1.0e-6 : 0.0) << endl;
    out << setw(12) << "apply()" << setw(5) << getNumApply()    << "    " << setw(15) << getApplyTime() << "    " 
        << setw(15) << getApplyFlops() << "    " 
        << setw(15) << (getApplyTime() != 0.0 ? getApplyFlops() / getApplyTime() * 1.0e-6 : 0.0) << endl;
    out << "===============================================================================" << std::endl;
    out << endl;
  }
}

}//namespace Ifpack2

#endif // IFPACK2_CHEBYSHEV_DEF_HPP

