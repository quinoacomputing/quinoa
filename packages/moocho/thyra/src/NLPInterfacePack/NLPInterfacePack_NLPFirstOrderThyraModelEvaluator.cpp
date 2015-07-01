// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
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

#include <assert.h>

#include <algorithm>

#include "NLPInterfacePack_NLPFirstOrderThyraModelEvaluator.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "AbstractLinAlgPack_ThyraAccessors.hpp"
#include "AbstractLinAlgPack_VectorSpaceThyra.hpp"
#include "AbstractLinAlgPack_VectorMutableThyra.hpp"
#include "AbstractLinAlgPack_MatrixOpNonsingThyra.hpp"
#include "AbstractLinAlgPack_BasisSystemComposite.hpp"
#include "AbstractLinAlgPack_VectorSpaceBlocked.hpp"
#include "AbstractLinAlgPack_VectorAuxiliaryOps.hpp"
#include "AbstractLinAlgPack_MatrixSymPosDefCholFactor.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace NLPInterfacePack {

NLPFirstOrderThyraModelEvaluator::NLPFirstOrderThyraModelEvaluator()
{}

NLPFirstOrderThyraModelEvaluator::NLPFirstOrderThyraModelEvaluator(
  const Teuchos::RCP<Thyra::ModelEvaluator<value_type> >  &model   
  ,const int                                                      p_idx 
  ,const int                                                      g_idx 
  )
{
  initialize(model,p_idx,g_idx);
}

void NLPFirstOrderThyraModelEvaluator::initialize(
  const Teuchos::RCP<Thyra::ModelEvaluator<value_type> >  &model
  ,const int                                                      p_idx
  ,const int                                                      g_idx
  )
{
  initializeBase(model,p_idx,g_idx);
}
  
// Overridden public members from NLP

void NLPFirstOrderThyraModelEvaluator::initialize(bool test_setup)
{
  if(initialized_) {
    NLPFirstOrder::initialize(test_setup);
    return;
  }
  NLPThyraModelEvaluatorBase::initialize(test_setup);
  NLPFirstOrder::initialize(test_setup);
}

void NLPFirstOrderThyraModelEvaluator::unset_quantities()
{
  NLPFirstOrder::unset_quantities();
}

// Overridden public members from NLPFirstOrder

void NLPFirstOrderThyraModelEvaluator::set_Gc(MatrixOp* Gc)
{
  NLPFirstOrder::set_Gc(Gc);
  Gc_updated_ = false;
}

const NLPFirstOrder::mat_fcty_ptr_t
NLPFirstOrderThyraModelEvaluator::factory_Gc() const
{
  return factory_Gc_;
}

const NLPFirstOrder::basis_sys_ptr_t
NLPFirstOrderThyraModelEvaluator::basis_sys() const
{
  return basis_sys_;
}

// Overridden protected members from NLPFirstOrder

void NLPFirstOrderThyraModelEvaluator::imp_calc_Gc(const Vector& x, bool newx, const FirstOrderInfo& first_order_info) const
{
  evalModel(x,newx,NULL,NULL,&first_order_info);
}

// private

void NLPFirstOrderThyraModelEvaluator::evalModel( 
  const Vector            &x
  ,bool                   newx
  ,const ZeroOrderInfo    *zero_order_info
  ,const ObjGradInfo      *obj_grad_info
  ,const FirstOrderInfo   *first_order_info
  ) const
{
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  using Teuchos::dyn_cast;
  using Teuchos::RCP;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using AbstractLinAlgPack::VectorMutableThyra;
  using AbstractLinAlgPack::MatrixOpThyra;
  using AbstractLinAlgPack::MatrixOpNonsingThyra;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::VerboseObjectTempState<MEB> VOTSME;
  typedef MEB::DerivativeMultiVector<value_type> DerivMV;
  typedef MEB::Derivative<value_type> Deriv;
  //
  // Get output and verbosity
  //
  const Teuchos::RCP<Teuchos::FancyOStream>
    out = this->getOStream();
  const Teuchos::EVerbosityLevel
    verbLevel = ( showModelEvaluatorTrace() ? this->getVerbLevel() : Teuchos::VERB_NONE );
  Teuchos::OSTab tab(out);
  VOTSME modelOutputTempState(model_,out,verbLevel);
  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out << "\nEntering MoochoPack::NLPFirstOrderThyraModelEvaluator::calc_point(...) ...\n";
  //
  // Set the input and output arguments
  //
  MEB::InArgs<value_type>  model_inArgs  = model_->createInArgs();
  MEB::OutArgs<value_type> model_outArgs = model_->createOutArgs();
  MatrixOp            *Gc = NULL;
  VectorMutable       *Gf = NULL;
  value_type          *f  = NULL;
  VectorMutable       *c  = NULL;
  preprocessBaseInOutArgs(
    x,newx,zero_order_info,obj_grad_info,first_order_info
    ,&model_inArgs,&model_outArgs,&Gc,&Gf,&f,&c
    );
  //
  MatrixOpNonsing  *C_aggr;
  MatrixOp         *N_aggr;
  if( Gc && !Gc_updated_ ) {
    BasisSystemComposite::get_C_N( Gc, &C_aggr, &N_aggr ); // Will return NULLs if Gc is not initialized
    if(C_aggr) {
      model_outArgs.set_W(
        rcp_const_cast<Thyra::LinearOpWithSolveBase<value_type> >(
          dyn_cast<MatrixOpNonsingThyra>(*C_aggr).set_uninitialized()
          ).assert_not_null()
        );
      if(p_idx_ >= 0) {
        // ToDo: This is implemented for direct sensitivities, change this for adjoint sensitivities
        model_outArgs.set_DfDp(
          p_idx_
          ,DerivMV(
            rcp_const_cast<Thyra::MultiVectorBase<value_type> >(
              rcp_dynamic_cast<const Thyra::MultiVectorBase<value_type> >(
                dyn_cast<MatrixOpThyra>(*N_aggr).set_uninitialized()
                )
              ).assert_not_null()
            ,MEB::DERIV_MV_BY_COL
            )
          );
      }
    }
    else {
      model_outArgs.set_W(model_->create_W().assert_not_null());
      if(p_idx_>=0)
        model_outArgs.set_DfDp(p_idx_,Thyra::create_DfDp_mv(*model_,p_idx_,MEB::DERIV_MV_BY_COL));
    }
    if (model_inArgs.supports(MEB::IN_ARG_alpha))
      model_inArgs.set_alpha(0.0);
    if (model_inArgs.supports(MEB::IN_ARG_beta))
      model_inArgs.set_beta(1.0);
  }
  //
  // Evaluate the model
  //
  model_->evalModel(model_inArgs,model_outArgs);
  //
  // Postprocess the output arguments
  //
  postprocessBaseOutArgs(&model_outArgs,Gf,f,c);
  //
  if( Gc && !Gc_updated_ ) {
    RCP<MatrixOpNonsing> C_ptr;
    RCP<MatrixOp>        N_ptr;
    if(!C_aggr) {
      C_ptr  = Teuchos::rcp(new MatrixOpNonsingThyra());
      C_aggr = &*C_ptr;
      if(p_idx_>=0) {
        N_ptr  = Teuchos::rcp(new MatrixOpThyra());
        N_aggr = &*N_ptr;
      }
    }
    RCP<Thyra::LinearOpWithSolveBase<value_type> >
      model_W = model_outArgs.get_W();
    model_W->setOStream(out);
    if(showModelEvaluatorTrace())
      model_W->setVerbLevel(verbLevel);
    dyn_cast<MatrixOpNonsingThyra>(*C_aggr).initialize(model_W,BLAS_Cpp::no_trans);
    // ToDo: This is implemented for direct sensitivities, change this for adjoint sensitivities
    if(p_idx_>=0)
      dyn_cast<MatrixOpThyra>(*N_aggr).initialize(model_outArgs.get_DfDp(p_idx_).getDerivativeMultiVector().getMultiVector(),BLAS_Cpp::no_trans);
    if( C_ptr.get() ) {
      BasisSystemComposite::initialize_Gc(
        this->space_x(), basis_sys_->var_dep(), basis_sys_->var_indep()
        ,this->space_c()
        ,C_ptr, N_ptr
        ,Gc
        );
    }
    Gc_updated_ = true;
  }
  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out << "\nLeaving MoochoPack::NLPFirstOrderThyraModelEvaluator::calc_point(...) ...\n";
}

}	// end namespace NLPInterfacePack
