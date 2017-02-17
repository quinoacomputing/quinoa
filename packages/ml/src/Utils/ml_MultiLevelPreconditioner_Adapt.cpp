/*!
 *
 *  \file ml_MultiLevelPreconditioner_Adapt.cpp
 *
 *  \brief Methods to define adaptive smoothed aggregation.
 *
 */
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#include "ml_common.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)

#include "ml_include.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"
#include "Epetra_Time.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "ml_epetra.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_viz_xyz.h"
#include "ml_agg_info.h"
#ifdef HAVE_ML_ANASAxI
#include "ml_anasazi.h"
#endif

// ================================================ ====== ==== ==== == =
// the tentative null space is in input because the user
// has to remember to allocate and fill it, and then to delete
// it after calling this method.
int ML_Epetra::MultiLevelPreconditioner::
ComputeAdaptivePreconditioner(int TentativeNullSpaceSize,
                              double* TentativeNullSpace)
{

  if ((TentativeNullSpaceSize == 0) || (TentativeNullSpace == 0))
    ML_CHK_ERR(-1);

  // ================================== //
  // get parameters from the input list //
  // ================================== //

  // maximum number of relaxation sweeps
  int MaxSweeps = List_.get("adaptive: max sweeps", 10);
  // number of std::vector to be added to the tentative null space
  int NumAdaptiveVectors = List_.get("adaptive: num vectors", 1);

  if (verbose_) {
    std::cout << PrintMsg_ << "*** Adaptive Smoother Aggregation setup ***" << std::endl;
    std::cout << PrintMsg_ << "    Maximum relaxation sweeps     = " << MaxSweeps << std::endl;
    std::cout << PrintMsg_ << "    Additional vectors to compute = " << NumAdaptiveVectors << std::endl;
  }

  // ==================================================== //
  // compute the preconditioner, set null space from user //
  // (who will have to delete std::vector TentativeNullSpace)  //
  // ==================================================== //

  double* NewNullSpace = 0;
  double* OldNullSpace = TentativeNullSpace;
  int OldNullSpaceSize = TentativeNullSpaceSize;

  // need some work otherwise matvec() with Epetra_Vbr fails.
  // Also, don't differentiate between range and domain here
  // as ML will not work if range != domain
  const Epetra_VbrMatrix* VbrA = NULL;
  VbrA = dynamic_cast<const Epetra_VbrMatrix*>(RowMatrix_);

  Epetra_Vector* LHS = 0;
  Epetra_Vector* RHS = 0;

  if (VbrA != 0) {
    LHS = new Epetra_Vector(VbrA->DomainMap());
    RHS = new Epetra_Vector(VbrA->DomainMap());
  } else {
    LHS = new Epetra_Vector(RowMatrix_->OperatorDomainMap());
    RHS = new Epetra_Vector(RowMatrix_->OperatorDomainMap());
  }

  // destroy what we may already have
  if (IsComputePreconditionerOK_ == true) {
    DestroyPreconditioner();
  }

  // build the preconditioner for the first time
  List_.set("null space: type", "pre-computed");
  List_.set("null space: dimension", OldNullSpaceSize);
  List_.set("null space: vectors", OldNullSpace);
  ComputePreconditioner();

  // ====================== //
  // add one std::vector at time //
  // ====================== //

  for (int istep = 0 ; istep < NumAdaptiveVectors ; ++istep) {

    if (verbose_) {
      std::cout << PrintMsg_ << "\tAdaptation step " << istep << std::endl;
      std::cout << PrintMsg_ << "\t---------------" << std::endl;
    }

    // ==================== //
    // look for "bad" modes //
    // ==================== //

    // note: should an error occur, ML_CHK_ERR will return,
    // and LHS and RHS will *not* be delete'd (--> memory leak).
    // Anyway, this means that something wrong happened in the code
    // and should be fixed by the user.

    LHS->Random();
    double Norm2;

    for (int i = 0 ; i < MaxSweeps ; ++i) {
      // RHS = (I - ML^{-1} A) LHS
      ML_CHK_ERR(RowMatrix_->Multiply(false,*LHS,*RHS));
      // FIXME: can do something slightly better here
      ML_CHK_ERR(ApplyInverse(*RHS,*RHS));
      ML_CHK_ERR(LHS->Update(-1.0,*RHS,1.0));
      LHS->Norm2(&Norm2);
      if (verbose_) {
        std::cout << PrintMsg_ << "\titer " << i << ", ||x||_2 = ";
        std::cout << Norm2 << std::endl;
      }
    }

    // scaling vectors
    {
      double theNormInf;
      LHS->NormInf(&theNormInf);
      LHS->Scale(1.0 / theNormInf);
    }

    // ========================================================= //
    // copy tentative and computed null space into NewNullSpace, //
    // which now becomes the standard null space                 //
    // ========================================================= //

    int NewNullSpaceSize = OldNullSpaceSize + 1;
    NewNullSpace = new double[NumMyRows() * NewNullSpaceSize];
    assert (NewNullSpace != 0);
    int itmp = OldNullSpaceSize * NumMyRows();
    for (int i = 0 ; i < itmp ; ++i) {
      NewNullSpace[i] = OldNullSpace[i];
    }

    for (int j = 0 ; j < NumMyRows() ; ++j) {
      NewNullSpace[itmp + j] = (*LHS)[j];
    }

    // =============== //
    // visualize modes //
    // =============== //

    if (List_.get("adaptive: visualize", false)) {

      double* x_coord = List_.get("viz: x-coordinates", (double*)0);
      double* y_coord = List_.get("viz: y-coordinates", (double*)0);
      double* z_coord = List_.get("viz: z-coordinates", (double*)0);
      assert (x_coord != 0);

      std::vector<double> plot_me(NumMyRows()/NumPDEEqns_);
      ML_Aggregate_Viz_Stats info;
      info.Amatrix = &(ml_->Amat[LevelID_[0]]);
      info.x = x_coord;
      info.y = y_coord;
      info.z = z_coord;
      info.Nlocal = NumMyRows() / NumPDEEqns_;
      info.Naggregates = 1;
      ML_Operator_AmalgamateAndDropWeak(&(ml_->Amat[LevelID_[0]]),
                                        NumPDEEqns_, 0.0);

      for (int ieqn = 0 ; ieqn < NumPDEEqns_ ; ++ieqn) {
        for (int j = 0 ; j < NumMyRows() ; j+=NumPDEEqns_) {
          plot_me[j / NumPDEEqns_] = (*LHS)[j + ieqn];
        }
        char FileName[80];
        sprintf(FileName,"nullspace-mode%d-eq%d.xyz", istep, ieqn);
        if (verbose_)
          std::cout << PrintMsg_ << "writing file " << FileName << "..." << std::endl;
        ML_Aggregate_VisualizeXYZ(info,FileName,
                                  ml_->comm,&plot_me[0]);
      }

      ML_Operator_UnAmalgamateAndDropWeak(&(ml_->Amat[LevelID_[0]]),
                                          NumPDEEqns_, 0.0);
    }

    // Destroy the old preconditioner
    DestroyPreconditioner();

    // ==================================================== //
    // build the new preconditioner with the new null space //
    // ==================================================== //

    List_.set("null space: type", "pre-computed");
    List_.set("null space: dimension", NewNullSpaceSize);
    List_.set("null space: vectors", NewNullSpace);

    ML_CHK_ERR(ComputePreconditioner());

    if (istep && (istep != NumAdaptiveVectors))
      delete OldNullSpace;

    OldNullSpace = NewNullSpace;
    OldNullSpaceSize = NewNullSpaceSize;

  }

  // keep trace of this pointer, it will be delete'd later
  NullSpaceToFree_ = NewNullSpace;

  delete LHS;
  delete RHS;

  return(0);

}

#endif /*ifdef HAVE_ML_EPETRA && HAVE_ML_TEUCHOS*/
#if NOT_DEFINED
    else if (AdaptType == "Anasazi") {

      // FIXME: right now it works for symmetric problems only
      double tol = List_.get("eigen-analysis: tolerance", 1e-5);

      Teuchos::ParameterList AnasaziList;
      AnasaziList.set("eigen-analysis: matrix operation", "I-ML^{-1}A");
      // not so sure about that
      AnasaziList.set("eigen-analysis: use diagonal scaling", false);
      AnasaziList.set("eigen-analysis: symmetric problem", true);
      // FIXME: what is length and block size???
      AnasaziList.set("eigen-analysis: length", 10);
      AnasaziList.set("eigen-analysis: block-size", 1);
      AnasaziList.set("eigen-analysis: tolerance", tol);
      AnasaziList.set("eigen-analysis: action", "LM");
      AnasaziList.set("eigen-analysis: restart", 1);
      AnasaziList.set("eigen-analysis: output", 10);

      // data to hold real and imag for eigenvalues and eigenvectors
      std::vector<double> RealEigenvalues(BlockSize);
      std::vector<double> ImagEigenvalues(BlockSize);

      std::vector<double> RealEigenvectors(BlockSize * NumMyRows());

      // this is the starting value -- random
      Epetra_MultiVector EigenVectors(OperatorDomainMap(),BlockSize);
      EigenVectors.Random();

      int NumRealEigenvectors = 0, NumImagEigenvectors = 0;

#ifdef HAVE_ML_ANASAxI
      // 2.- call Anasazi and store the results in eigenvectors
      ML_Anasazi::Interface(RowMatrix_,EigenVectors,&RealEigenvalues[0],
                            &ImagEigenvalues[0], AnasaziList, &RealEigenvectors[0], 0,
                            &NumRealEigenvectors, &NumImagEigenvectors, ml_);
#else
      if( Comm().MyPID() == 0 ) {
        std::cerr << ErrorMsg_ << "ML has been configure without the Anasazi interface" << std::endl
          << ErrorMsg_ << "You must add the option --enable-anasazi to use" << std::endl
          << ErrorMsg_ << "filtering and Anasazi" << std::endl;
      }
      ML_EXIT(EXIT_FAILURE);
#endif

      assert (NumRealEigenvectors != 0);

      if (verbose_) {
        std::cout << PrintMsg_ << "\t- Computed eigenvalues of I - ML^{-1}A:" << std::endl;
        for (int i = 0 ; i < BlockSize ; ++i) {
          std::cout << PrintMsg_ << "\t  z = " << std::setw(10) << RealEigenvalues[i]
            << " + i(" << std::setw(10) << ImagEigenvalues[i] << " ),  |z| = "
            << std::sqrt(RealEigenvalues[i]*RealEigenvalues[i] + ImagEigenvalues[i]*ImagEigenvalues[i]) << std::endl;
        }
        std::cout << PrintMsg_ << "\t- Using " << NumRealEigenvectors << " real and "
          << NumImagEigenvectors << " imaginary eigenvector(s)" << std::endl;
      }

      // FIXME: this is not very efficient...
      for (int i = 0 ; i < NumRealEigenvectors ; ++i) {
        for (int j = 0 ; j < NumMyRows() ; ++j)
          (*LHS)[i][j] = RealEigenvectors[j + i * NumMyRows()];
      }
    }
    else {
      std::cerr << ErrorMsg_ << "`adaptive: type' has an incorrect value" << std::endl;
      std::cerr << ErrorMsg_ << "(" << AdaptType << "). It should be:" << std::endl;
      std::cerr << ErrorMsg_ << "<relaxation> / <Anasazi>" << std::endl;
      ML_EXIT(-1);
    }
#endif
