// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef XPETRA_EPETRAOPERATOR_HPP
#define XPETRA_EPETRAOPERATOR_HPP

#include "Xpetra_EpetraConfigDefs.hpp"

#include <Epetra_Operator.hpp>

#include "Xpetra_Map.hpp"
#include "Xpetra_EpetraMap.hpp"
#include "Xpetra_MultiVector.hpp"
#include "Xpetra_EpetraMultiVector.hpp"
#include "Xpetra_Operator.hpp"

#include "Xpetra_Utils.hpp"

namespace Xpetra {

  template<class EpetraGlobalOrdinal>
  class EpetraOperator : public Operator<double, int, EpetraGlobalOrdinal>
  {
    typedef double                                                      Scalar;
    typedef int                                                         LocalOrdinal;
    typedef EpetraGlobalOrdinal                                         GlobalOrdinal;
    typedef typename Operator<double, int, GlobalOrdinal>::node_type    Node;

  public:
    //@{

    //! The Map associated with the domain of this operator, which must be compatible with X.getMap().
    virtual Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getDomainMap() const {
      XPETRA_MONITOR("EpetraOperator::getDomainMap()");
      return toXpetra(op_->OperatorDomainMap());
    }

    //! The Map associated with the range of this operator, which must be compatible with Y.getMap().
    virtual Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getRangeMap() const {
      XPETRA_MONITOR("EpetraOperator::getRangeMap()");
      return toXpetra(op_->OperatorRangeMap());
    }

    //! \brief Computes the operator-multivector application.
    /*! Loosely, performs \f$Y = \alpha \cdot A^{\textrm{mode}} \cdot X + \beta \cdot Y\f$. However, the details of operation
        vary according to the values of \c alpha and \c beta. Specifically
        - if <tt>beta == 0</tt>, apply() <b>must</b> overwrite \c Y, so that any values in \c Y (including NaNs) are ignored.
        - if <tt>alpha == 0</tt>, apply() <b>may</b> short-circuit the operator, so that any values in \c X (including NaNs) are ignored.
     */
    virtual void
    apply (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
           MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
           Teuchos::ETransp mode = Teuchos::NO_TRANS,
           Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
           Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const {
      XPETRA_MONITOR("EpetraOperator::apply");

      XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT<GlobalOrdinal>, X, eX, "Xpetra::EpetraOperator->apply(): cannot cast input to Xpetra::EpetraMultiVectorT");
      XPETRA_DYNAMIC_CAST(      EpetraMultiVectorT<GlobalOrdinal>, Y, eY, "Xpetra::EpetraOperator->apply(): cannot cast input to Xpetra::EpetraMultiVectorT");

      TEUCHOS_TEST_FOR_EXCEPTION((mode != Teuchos::NO_TRANS) && (mode != Teuchos::TRANS), Exceptions::NotImplemented,
                                 "Xpetra::EpetraOperator->apply(): can only accept mode == NO_TRANS or mode == TRANS");
      TEUCHOS_TEST_FOR_EXCEPTION(mode == Teuchos::TRANS && !hasTransposeApply, Exceptions::RuntimeError,
                                 "Xpetra::EpetraOperator->apply(): cannot apply transpose as underlying Epetra operator does not support it");

      // Helper vector for string A*X
      RCP<Epetra_MultiVector> epY = eY.getEpetra_MultiVector();
      RCP<Epetra_MultiVector> tmp = Teuchos::rcp(new Epetra_MultiVector(*epY));
      tmp->PutScalar(0.0);

      bool curTranspose = op_->UseTranspose();
      op_->SetUseTranspose(mode == Teuchos::TRANS);
      XPETRA_ERR_CHECK(op_->Multiply(*eX.getEpetra_MultiVector(), *tmp));
      op_->setUseTranspose(curTranspose);

      // calculate alpha * A * x + beta * y
      XPETRA_ERR_CHECK(eY.getEpetra_MultiVector()->Update(alpha, *tmp, beta));
    }

    /// \brief Whether this operator supports applying the transpose or conjugate transpose.
    virtual bool hasTransposeApply() const {
      if (op_->UseTransepose()) {
        // If current setting is to use transpose, then obviously we can use it
        return true;
      }
      // We do not currently use transpose, try setting it
      int err = op_->SetUseTranspose(true);
      SetUseTranspose(false);
      return (err == 0);
    }

    //@}

    //! @name Overridden from Teuchos::Describable
    //@{

    //! A simple one-line description of this object.
    std::string description() const { XPETRA_MONITOR("EpetraOperator::description"); return "Epetra_Operator"; }

    //! Print the object with the given verbosity level to a FancyOStream.
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const {
      XPETRA_MONITOR("EpetraOperator::describe");
      out << "Epetra_Operator" << std::endl;
    }

    //@}

    //! @name Xpetra specific
    //@{

    //! TpetraOperator constructor to wrap a Tpetra::Operator object
    EpetraOperator(const Teuchos::RCP<Tpetra::Operator< Scalar, LocalOrdinal, GlobalOrdinal, Node> > &op) : op_(op) { } //TODO removed const

    //@}

  private:
    //! The Tpetra::Operator which this class wraps.
    RCP< Epetra_Operator> op_;

  }; // TpetraOperator class


} // Xpetra namespace

#endif // XPETRA_EPETRAOPERATOR_HPP
