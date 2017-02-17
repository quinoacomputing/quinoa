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

#ifndef IFPACK2_AMESOS_HPP
#define IFPACK2_AMESOS_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Tpetra_Operator.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"

class Tpetra_Map;
class Tpetra_Time;
class Tpetra_Comm;
class Amesos_BaseSolver;
class Tpetra_LinearProblem;
class Tpetra_RowMatrix;

//! Ifpack2_Amesos: a class to use Amesos' factorizations as preconditioners.
/*!
Class Ifpack2_Amesos enables the use of Amesos' factorizations as 
Ifpack2_Preconditioners.

Ifpack2_Amesos is just a bare-bone wrap to Amesos. Currently, the
only parameter required recognized by SetParameters() is
\c "amesos: solver type" (defaulted to \c
"Amesos_Klu"), which defined the Amesos solver. The Teuchos list
in input to SetParameters() is copied, then the copied list is
used to set the parameters of the Amesos object.

This class works with matrices whose communicator contains only one
process, that is, either serial matrices, or Ifpack2_LocalFilter'd matrices.

\warning The number of flops is NOT updated.

\author Michael Heroux, SNL 9214.

\date Last update Oct-04.

*/
class Ifpack2_Amesos : public Ifpack2_Preconditioner {
      
public:

  //@{ \name Constructors/Destructors.

  //! Constructor.
  Ifpack2_Amesos(Tpetra_RowMatrix* Matrix);

  //! Copy constructor.
  Ifpack2_Amesos(const Ifpack2_Amesos& rhs);

  //! Operator=.
  Ifpack2_Amesos& operator=(const Ifpack2_Amesos& rhs);

  //@{ \name Destructor.
  //! Destructor
  virtual ~Ifpack2_Amesos() {};

  //@}

  //@{ \name Attribute set methods.

   //! If set true, transpose of this operator will be applied (not implemented).
    /*! This flag allows the transpose of the given operator to be used 
     * implicitly.  
      
    \param 
	   UseTranspose_in - (In) If true, multiply by the transpose of operator, 
	   otherwise just use operator.

    \return Integer error code, set to 0 if successful.  Set to -1 if this implementation does not support transpose.
  */

  virtual int SetUseTranspose(bool UseTranspose_in);
  //@}
  
  //@{ \name Mathematical functions.

    //! Applies the matrix to an Tpetra_MultiVector.
  /*! 
    \param
    X - (In) A Tpetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param 
    Y - (Out) A Tpetra_MultiVector of dimension NumVectors containing the result.

    \return Integer error code, set to 0 if successful.
    */
    virtual int Apply(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

    //! Applies the preconditioner to X, returns the result in Y.
  /*! 
    \param 
    X - (In) A Tpetra_MultiVector of dimension NumVectors to be preconditioned.
    \param 
    Y - (Out) A Tpetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.

    \warning In order to work with AztecOO, any implementation of this method 
    must support the case where X and Y are the same object.
    */
    virtual int ApplyInverse(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

    //! Returns the infinity norm of the global matrix (not implemented)
    virtual double NormInf() const;
  //@}
  
  //@{ \name Attribute access functions

    //! Returns a character string describing the operator
    virtual const char * Label() const;

    //! Returns the current UseTranspose setting.
    virtual bool UseTranspose() const;

    //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
    virtual bool HasNormInf() const;

    //! Returns a pointer to the Tpetra_Comm communicator associated with this operator.
    virtual const Tpetra_Comm & Comm() const;

    //! Returns the Tpetra_Map object associated with the domain of this operator.
    virtual const Tpetra_Map & OperatorDomainMap() const;

    //! Returns the Tpetra_Map object associated with the range of this operator.
    virtual const Tpetra_Map & OperatorRangeMap() const;

  //@}

  //@{ \name Construction and application methods.
 
  //! Returns \c true is the preconditioner has been successfully initialized.
  virtual bool IsInitialized() const
  {
    return(IsInitialized_);
  }

  //! Initializes the preconditioners.
  /*! \return
   * 0 if successful, 1 if problems occurred.
   */
  virtual int Initialize();

  //! Returns \c true if the preconditioner has been successfully computed.
  virtual bool IsComputed() const
  {
    return(IsComputed_);
  }

  //! Computes the preconditioners.
  /*! \return
   * 0 if successful, 1 if problems occurred.
   */
  virtual int Compute();

  //! Sets all the parameters for the preconditioner.
  /*! Parameters currently supported:
   * - \c "amesos: solver type" : Specifies the solver type
   *   for Amesos. Default: \c Amesos_Klu.
   *
   * The input list will be copied, then passed to the Amesos
   * object through Amesos::SetParameters().
   */   
  virtual int SetParameters(Teuchos::ParameterList& List);

  //@}

  //@{ \name Query methods.

  //! Returns a const reference to the internally stored matrix.
  virtual const Tpetra_RowMatrix& Matrix() const
  {
    return(*Matrix_);
  }

  //! Returns the estimated condition number, computes it if necessary.
  virtual double Condest(const Ifpack2_CondestType CT = Ifpack2_Cheap,
                         const int MaxIters = 1550,
                         const double Tol = 1e-9,
			 Tpetra_RowMatrix* Matrix_in= 0);
  
  //! Returns the estimated condition number, never computes it.
  virtual double Condest() const
  {
    return(Condest_);
  }

  //! Returns the number of calls to Initialize().
  virtual int NumInitialize() const
  {
    return(NumInitialize_);
  }

  //! Returns the number of calls to Compute().
  virtual int NumCompute() const
  {
    return(NumCompute_);
  }

  //! Returns the number of calls to ApplyInverse().
  virtual int NumApplyInverse() const
  {
    return(NumApplyInverse_);
  }

  //! Returns the total time spent in Initialize().
  virtual double InitializeTime() const
  {
    return(InitializeTime_);
  }

  //! Returns the total time spent in Compute().
  virtual double ComputeTime() const
  {
    return(ComputeTime_);
  }

  //! Returns the total time spent in ApplyInverse().
  virtual double ApplyInverseTime() const
  {
    return(ApplyInverseTime_);
  }

  //! Returns the number of flops in the initialization phase.
  virtual double InitializeFlops() const
  {
    return(0.0);
  }

  //! Returns the total number of flops to computate the preconditioner.
  virtual double ComputeFlops() const
  {
    return(ComputeFlops_);
  }

  //! Returns the total number of flops to apply the preconditioner.
  virtual double ApplyInverseFlops() const
  {
    return(ApplyInverseFlops_);
  }

  // Returns a constant reference to the internally stored 
  virtual const Teuchos::ParameterList& List() const 
  {
    return(List_);
  }

  //! Prints on ostream basic information about \c this object.
  virtual std::ostream& Print(std::ostream& os) const;

  //@}

protected:
  
  //@{ \name Methods to get/set private data

  //! Sets the label.
  inline void SetLabel(const char* Label_in) 
  {
    Label_ = Label_in;
  }

  //! Sets \c IsInitialized_.
  inline void SetIsInitialized(const bool IsInitialized_in)
  {
    IsInitialized_ = IsInitialized_in;
  }

  //! Sets \c IsComputed_.
  inline void SetIsComputed(const int IsComputed_in)
  {
    IsComputed_ = IsComputed_in;
  }

  //! Sets \c NumInitialize_.
  inline void SetNumInitialize(const int NumInitialize_in)
  {
    NumInitialize_ = NumInitialize_in;
  }

  //! Sets \c NumCompute_.
  inline void SetNumCompute(const int NumCompute_in)
  {
    NumCompute_ = NumCompute_in;
  }

  //! Sets \c NumApplyInverse_.
  inline void SetNumApplyInverse(const int NumApplyInverse_in)
  {
    NumApplyInverse_ = NumApplyInverse_in;
  }

  //! Sets \c InitializeTime_.
  inline void SetInitializeTime(const double InitializeTime_in)
  {
    InitializeTime_ = InitializeTime_in;
  }

  //! Sets \c ComputeTime_.
  inline void SetComputeTime(const double ComputeTime_in)
  {
    ComputeTime_ = ComputeTime_in;
  }

  //! Sets \c ApplyInverseTime_.
  inline void SetApplyInverseTime(const double ApplyInverseTime_in)
  {
    ApplyInverseTime_ = ApplyInverseTime_in;
  }

  //! Sets \c ComputeFlops_.
  inline void SetComputeFlops(const double ComputeFlops_in)
  {
    ComputeFlops_ = ComputeFlops_in;
  }

  //! Sets \c ComputeFlops_.
  inline void SetApplyInverseFlops(const double ApplyInverseFlops_in)
  {
    ApplyInverseFlops_ = ApplyInverseFlops_in;
  }

  //! Set \c List_.
  inline void SetList(const Teuchos::ParameterList& List_in)
  {
    List_ = List_in;
  }
  //@}
  
private:

  //! Pointers to the matrix to be preconditioned.
  Teuchos::RCP<const Tpetra_RowMatrix> Matrix_;

  //! Linear problem required by Solver_.
  Teuchos::RCP<Tpetra_LinearProblem> Problem_;
  //! Amesos solver, use to apply the inverse of the local matrix.
  Teuchos::RCP<Amesos_BaseSolver> Solver_;
  //! Contains a copy of the input parameter list.
  Teuchos::ParameterList List_;

  //! Contains the label of \c this object.
  string Label_;
  //! If true, the preconditioner has been successfully initialized.
  bool IsInitialized_;
  //! If true, the preconditioner has been successfully computed.
  bool IsComputed_;
  //! If true, the preconditioner solves for the transpose of the matrix.
  bool UseTranspose_;

  //! Contains the number of successful calls to Initialize().
  int NumInitialize_;
  //! Contains the number of successful call to Compute().
  int NumCompute_;
  //! Contains the number of successful call to ApplyInverse().
  mutable int NumApplyInverse_;

  //! Contains the time for all successful calls to Initialize().
  double InitializeTime_;
  //! Contains the time for all successful calls to Compute().
  double ComputeTime_;
  //! Contains the time for all successful calls to ApplyInverse().
  mutable double ApplyInverseTime_;
  //! Time object.
  Teuchos::RCP<Tpetra_Time> Time_;

  //! Contains the number of flops for Compute().
  double ComputeFlops_;
  //! Contain sthe number of flops for ApplyInverse().
  double ApplyInverseFlops_;

  //! Contains the estimated condition number.
  double Condest_;
};

#endif // IFPACK2_AMESOS_HPP
