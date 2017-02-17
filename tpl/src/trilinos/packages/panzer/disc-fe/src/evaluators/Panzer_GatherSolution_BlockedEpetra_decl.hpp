// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_GATHER_SOLUTION_BLOCKEDEPETRA_DECL_HPP
#define PANZER_EVALUATOR_GATHER_SOLUTION_BLOCKEDEPETRA_DECL_HPP

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "PanzerDiscFE_config.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_CloneableEvaluator.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

// thyra forward decl
namespace Thyra {
  template <typename> class ProductVectorBase;
}

namespace panzer {

class BlockedEpetraLinearObjContainer;

template <typename LocalOrdinalT,typename GlobalOrdinalT>
class UniqueGlobalIndexer; //forward declaration

/** \brief Gathers solution values from the Newton solution vector into
    the nodal fields of the field manager

    Currently makes an assumption that the stride is constant for dofs
    and that the nmber of dofs is equal to the size of the solution
    names vector.
*/
template<typename EvalT, typename TRAITS,typename LO,typename GO> class GatherSolution_BlockedEpetra
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Residual, TRAITS>,
    public panzer::CloneableEvaluator  {
public:
   typedef typename EvalT::ScalarT ScalarT;

   GatherSolution_BlockedEpetra(const Teuchos::ParameterList& p) {}

   virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
   { return Teuchos::rcp(new GatherSolution_BlockedEpetra<EvalT,TRAITS,LO,GO>(pl)); }

   void postRegistrationSetup(typename TRAITS::SetupData d, PHX::FieldManager<TRAITS>& vm)
   { }
   void evaluateFields(typename TRAITS::EvalData d)
   { std::cout << "unspecialized version of \"GatherSolution_BlockedEpetra::evaluateFields\" on "+PHX::typeAsString<EvalT>()+"\" should not be used!" << std::endl;
     TEUCHOS_ASSERT(false); }
};

// **************************************************************
// **************************************************************
// * Specializations
// **************************************************************
// **************************************************************


// **************************************************************
// Residual
// **************************************************************
template<typename TRAITS,typename LO,typename GO>
class GatherSolution_BlockedEpetra<panzer::Traits::Residual,TRAITS,LO,GO>
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Residual, TRAITS>,
    public panzer::CloneableEvaluator  {


public:

  GatherSolution_BlockedEpetra(const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > & indexers)
     : indexers_(indexers) {}

  GatherSolution_BlockedEpetra(const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > & indexers,
                               const Teuchos::ParameterList& p);

  void postRegistrationSetup(typename TRAITS::SetupData d,
                             PHX::FieldManager<TRAITS>& vm);

  void preEvaluate(typename TRAITS::PreEvalData d);

  void evaluateFields(typename TRAITS::EvalData d);

  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new GatherSolution_BlockedEpetra<panzer::Traits::Residual,TRAITS,LO,GO>(indexers_,pl)); }


private:

  typedef typename panzer::Traits::Residual EvalT;
  typedef typename panzer::Traits::Residual::ScalarT ScalarT;

  std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > indexers_;

  std::vector<int> indexerIds_;   // block index
  std::vector<int> subFieldIds_; // sub field numbers

  std::vector< PHX::MDField<ScalarT,Cell,NODE> > gatherFields_;

  std::vector<std::string> indexerNames_;
  bool useTimeDerivativeSolutionVector_;
  std::string globalDataKey_; // what global data does this fill?

  Teuchos::RCP<Thyra::ProductVectorBase<double> > x_;

  // Fields for storing tangent components dx/dp of solution vector x
  // These are not actually used by the residual specialization of this evaluator,
  // even if they are supplied, but it is useful to declare them as dependencies anyway
  // when saving the tangent components to the output file
  bool has_tangent_fields_;
  std::vector< std::vector< PHX::MDField<const ScalarT,Cell,NODE> > > tangentFields_;

  GatherSolution_BlockedEpetra();
};

// **************************************************************
// Tangent
// **************************************************************
template<typename TRAITS,typename LO,typename GO>
class GatherSolution_BlockedEpetra<panzer::Traits::Tangent,TRAITS,LO,GO>
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Tangent, TRAITS>,
    public panzer::CloneableEvaluator  {


public:

  GatherSolution_BlockedEpetra(const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > & indexers)
     : indexers_(indexers) {}

  GatherSolution_BlockedEpetra(const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > & indexers,
                               const Teuchos::ParameterList& p);

  void postRegistrationSetup(typename TRAITS::SetupData d,
                             PHX::FieldManager<TRAITS>& vm);

  void preEvaluate(typename TRAITS::PreEvalData d);

  void evaluateFields(typename TRAITS::EvalData d);

  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new GatherSolution_BlockedEpetra<panzer::Traits::Tangent,TRAITS,LO,GO>(indexers_,pl)); }


private:

  typedef typename panzer::Traits::Tangent EvalT;
  typedef typename panzer::Traits::Tangent::ScalarT ScalarT;

  std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > indexers_;

  std::vector<int> indexerIds_;   // block index
  std::vector<int> subFieldIds_; // sub field numbers

  std::vector< PHX::MDField<ScalarT,Cell,NODE> > gatherFields_;

  std::vector<std::string> indexerNames_;
  bool useTimeDerivativeSolutionVector_;
  std::string globalDataKey_; // what global data does this fill?

  Teuchos::RCP<Thyra::ProductVectorBase<double> > x_;

  // Fields for storing tangent components dx/dp of solution vector x
  bool has_tangent_fields_;
  std::vector< std::vector< PHX::MDField<const ScalarT,Cell,NODE> > > tangentFields_;

  GatherSolution_BlockedEpetra();
};

// **************************************************************
// Jacobian
// **************************************************************
template<typename TRAITS,typename LO,typename GO>
class GatherSolution_BlockedEpetra<panzer::Traits::Jacobian,TRAITS,LO,GO>
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Jacobian, TRAITS>,
    public panzer::CloneableEvaluator  {

public:

  GatherSolution_BlockedEpetra(const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > & indexers)
     : indexers_(indexers) {}

  GatherSolution_BlockedEpetra(const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > & indexers,
                               const Teuchos::ParameterList& p);

  void postRegistrationSetup(typename TRAITS::SetupData d,
     PHX::FieldManager<TRAITS>& vm);

  void preEvaluate(typename TRAITS::PreEvalData d);

  void evaluateFields(typename TRAITS::EvalData d);

  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new GatherSolution_BlockedEpetra<panzer::Traits::Jacobian,TRAITS,LO,GO>(indexers_,pl)); }

private:

  typedef typename panzer::Traits::Jacobian EvalT;
  typedef typename panzer::Traits::Jacobian::ScalarT ScalarT;

  std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > indexers_;

  std::vector<int> indexerIds_;   // block index
  std::vector<int> subFieldIds_; // sub field numbers

  std::vector< PHX::MDField<ScalarT,Cell,NODE> > gatherFields_;

  std::vector<std::string> indexerNames_;
  bool useTimeDerivativeSolutionVector_;
  bool disableSensitivities_;
  std::string sensitivitiesName_; // This sets which gather operations have sensitivities
  bool applySensitivities_;       // This is a local variable that is used by evaluateFields
                                  // to turn on/off a certain set of sensitivities
  std::string globalDataKey_; // what global data does this fill?
  int gatherSeedIndex_; // what gather seed in the workset to use
                        // if less than zero then use alpha or beta
                        // as appropriate

  Teuchos::RCP<Thyra::ProductVectorBase<double> > x_;

  GatherSolution_BlockedEpetra();
};

}

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
#include "Panzer_GatherSolution_BlockedEpetra_Hessian.hpp"
#endif

// **************************************************************
#endif
