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

#include "NLPInterfacePack_NLPDirectThyraModelEvaluator.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "AbstractLinAlgPack_ThyraAccessors.hpp"
#include "AbstractLinAlgPack_VectorSpaceThyra.hpp"
#include "AbstractLinAlgPack_VectorMutableThyra.hpp"
#include "AbstractLinAlgPack_MultiVectorMutableThyra.hpp"
#include "AbstractLinAlgPack_MatrixOpNonsingThyra.hpp"
#include "AbstractLinAlgPack_VectorSpaceBlocked.hpp"
#include "AbstractLinAlgPack_VectorAuxiliaryOps.hpp"
#include "AbstractLinAlgPack_VectorDenseEncap.hpp"
#include "AbstractLinAlgPack_BasisSystem.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace NLPInterfacePack {

NLPDirectThyraModelEvaluator::NLPDirectThyraModelEvaluator()
  :DfDp_is_const_(false)
{}

NLPDirectThyraModelEvaluator::NLPDirectThyraModelEvaluator(
  const Teuchos::RCP<Thyra::ModelEvaluator<value_type> > &model 
  ,const int p_idx 
  ,const int g_idx 
  ,const objDirecFiniteDiffCalculator_ptr_t objDirecFiniteDiffCalculator
  ,const conDirecFiniteDiffCalculator_ptr_t conDirecFiniteDiffCalculator
  )
  :DfDp_is_const_(false)
{
  initialize(
    model,p_idx,g_idx
    ,objDirecFiniteDiffCalculator,conDirecFiniteDiffCalculator
    );
}

void NLPDirectThyraModelEvaluator::initialize(
  const Teuchos::RCP<Thyra::ModelEvaluator<value_type> > &model
  ,const int p_idx
  ,const int g_idx
  ,const objDirecFiniteDiffCalculator_ptr_t objDirecFiniteDiffCalculator
  ,const conDirecFiniteDiffCalculator_ptr_t conDirecFiniteDiffCalculator
  )
{
  typedef Thyra::ModelEvaluatorBase MEB;
  if(objDirecFiniteDiffCalculator.get())
    this->set_objDirecFiniteDiffCalculator(objDirecFiniteDiffCalculator);
  if(conDirecFiniteDiffCalculator.get())
    this->set_conDirecFiniteDiffCalculator(conDirecFiniteDiffCalculator);
  initializeBase(model,p_idx,g_idx);
  Thyra::ModelEvaluatorBase::OutArgs<double> model_outArgs = model->createOutArgs();
  MEB::DerivativeProperties model_W_properties = model_outArgs.get_W_properties();
  if( p_idx >= 0 ) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      (conDirecFiniteDiffCalculator_.get()==0 && !model_outArgs.supports(MEB::OUT_ARG_DfDp,p_idx).supports(MEB::DERIV_MV_BY_COL))
      ,std::invalid_argument
      ,"Error, model must support computing DfDp("<<p_idx<<") as a"
      " column-oriented multi-vector if not using finite differences!"
      );
  }
  DfDp_is_const_  = false;
}

// Overridden public members from NLP

void NLPDirectThyraModelEvaluator::initialize(bool test_setup)
{
  if(initialized_) {
    NLPDirect::initialize(test_setup);
    return;
  }
  NLPThyraModelEvaluatorBase::initialize(test_setup);
  NLPDirect::initialize(test_setup);
}

void NLPDirectThyraModelEvaluator::unset_quantities()
{
  NLPDirect::unset_quantities();
}

// Overridden public members from NLPObjGrad

bool NLPDirectThyraModelEvaluator::supports_Gf() const
{
  if(objDirecFiniteDiffCalculator_.get())
    return false;
  return true;
  // ToDo: Change this to false when only the operator for model_DgDx is
  // supported!
}

bool NLPDirectThyraModelEvaluator::supports_Gf_prod() const
{
  return ( objDirecFiniteDiffCalculator_.get() != NULL );
}

value_type NLPDirectThyraModelEvaluator::calc_Gf_prod(
  const Vector& x, const Vector& d, bool newx
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPT(objDirecFiniteDiffCalculator_.get()==0);
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgs<value_type> basePoint = model_->createInArgs();
  MEB::InArgs<value_type> directions = model_->createInArgs();
  MEB::OutArgs<value_type> baseFunc = model_->createOutArgs();
  MEB::OutArgs<value_type> variations = model_->createOutArgs();
  set_x(x,&basePoint);
  set_x(d,&directions);
  variations.set_g(g_idx_,createMember(model_->get_g_space(g_idx_)));
  objDirecFiniteDiffCalculator_->calcVariations(
    *model_,basePoint,directions,baseFunc,variations
    );
  return Thyra::get_ele(*variations.get_g(g_idx_),0);
}

// Overridden public members from NLPDirect

Range1D NLPDirectThyraModelEvaluator::var_dep() const
{
  return basis_sys_->var_dep();
}

Range1D NLPDirectThyraModelEvaluator::var_indep() const
{
  return basis_sys_->var_indep();
}

const NLPDirect::mat_fcty_ptr_t
NLPDirectThyraModelEvaluator::factory_D() const
{
  return basis_sys_->factory_D();
}

const NLPDirect::mat_sym_nonsing_fcty_ptr_t
NLPDirectThyraModelEvaluator::factory_S() const
{
  return basis_sys_->factory_S();
}

void NLPDirectThyraModelEvaluator::calc_point(
  const Vector &x
  ,value_type *f
  ,VectorMutable *c
  ,bool recalc_c
  ,VectorMutable *Gf
  ,VectorMutable *py
  ,VectorMutable *rGf
  ,MatrixOp *GcU
  ,MatrixOp *D
  ,MatrixOp *Uz
  ) const
{
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  using Teuchos::dyn_cast;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using AbstractLinAlgPack::VectorSpaceThyra;
  using AbstractLinAlgPack::VectorMutableThyra;
  using AbstractLinAlgPack::MultiVectorMutableThyra;
  using AbstractLinAlgPack::MatrixOpThyra;
  using AbstractLinAlgPack::MatrixOpNonsingThyra;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::RCP<Thyra::VectorBase<value_type> > TVecPtr;
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
  typedef Teuchos::VerboseObjectTempState<MEB> VOTSME;
  VOTSME modelOutputTempState(model_,out,verbLevel);
  typedef Teuchos::VerboseObjectTempState<Thyra::DirectionalFiniteDiffCalculator<value_type> > VOTSDFDC;
  VOTSDFDC objDirecFiniteDiffCalculatorOutputTempState(objDirecFiniteDiffCalculator_,out,verbLevel);
  VOTSDFDC conDirecFiniteDiffCalculatorOutputTempState(conDirecFiniteDiffCalculator_,out,verbLevel);
  const bool trace = static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW);
  if(out.get() && trace)
    *out << "\nEntering MoochoPack::NLPDirectThyraModelEvaluator::calc_point(...) ...\n";
  //
  // Validate input
  //
  TEUCHOS_TEST_FOR_EXCEPT(GcU!=NULL); // Can't handle these yet!
  TEUCHOS_TEST_FOR_EXCEPT(Uz!=NULL);
  TEUCHOS_TEST_FOR_EXCEPTION(
    objDirecFiniteDiffCalculator_.get()!=NULL && Gf!=NULL, std::logic_error
    ,"Error, can not compute full gradient vector Gf when using directional finite differences!"
    );
  //
  // Compute using derivative objects that come from the underlying model.
  //
  // Set the input and output arguments
  //
  // ToDo: Disallow computation Gf and instead just compute
  // the operators for this!
  //
  MEB::InArgs<value_type> model_inArgs = model_->createInArgs();
  MEB::OutArgs<value_type> model_outArgs = model_->createOutArgs();
  NLPObjGrad::ObjGradInfo obj_grad_info;
  obj_grad_info.Gf = Gf;
  obj_grad_info.f = f;
  if(recalc_c) obj_grad_info.c = c;
  preprocessBaseInOutArgs(
    x,true,NULL,&obj_grad_info,NULL
    ,&model_inArgs,&model_outArgs,NULL,NULL,NULL,NULL
    );
  if( py || rGf || D ) {
    if(thyra_C_.get()==NULL)
      thyra_C_ = model_->create_W();
    model_outArgs.set_W(thyra_C_.assert_not_null());
  }
  bool new_thyra_N = false;
  if( rGf || D ) {
    if(thyra_N_.get()==NULL) {
      thyra_N_ = Thyra::create_DfDp_mv(*model_,p_idx_,MEB::DERIV_MV_BY_COL).getMultiVector();
      new_thyra_N = true;
    }
    if(
      !conDirecFiniteDiffCalculator_.get()
      &&
      ( new_thyra_N
        || model_outArgs.get_DfDp_properties(p_idx_).linearity!=MEB::DERIV_LINEARITY_CONST )
      )
    {
      model_outArgs.set_DfDp(p_idx_,DerivMV(thyra_N_.assert_not_null(),MEB::DERIV_MV_BY_COL));
    }
  }
  if ( !is_null(model_outArgs.get_W()) || !is_null(model_outArgs.get_W_op()) ) {
    if (model_inArgs.supports(MEB::IN_ARG_alpha))
      model_inArgs.set_alpha(0.0);
    if (model_inArgs.supports(MEB::IN_ARG_beta))
      model_inArgs.set_beta(1.0);
  }
  //
  // Evaluate the functions
  //
  model_->evalModel(model_inArgs,model_outArgs);
  //
  // Postprocess the evaluation
  //
  postprocessBaseOutArgs(&model_outArgs,Gf,f,recalc_c?c:NULL);
  // Setup solve components
  const VectorSpaceThyra *space_c;
  const VectorSpaceThyra *space_xD;
  Teuchos::RCP<const Thyra::VectorBase<value_type> > thyra_c;
  Teuchos::RCP<Thyra::VectorBase<value_type> > thyra_py;
  RCP<MatrixOp> D_used = rcp(D,false);
  RCP<Thyra::MultiVectorBase<value_type> > thyra_D;

  if( py || ( ( rGf || D ) && conDirecFiniteDiffCalculator_.get() ) ) {
    space_c = &dyn_cast<const VectorSpaceThyra>(c->space()),
      space_xD = &dyn_cast<const VectorSpaceThyra>(py->space());
    get_thyra_vector(*space_c,*c,&thyra_c);
    get_thyra_vector(*space_xD,py,&thyra_py);
  }

  if( D || rGf ) {
    if(!D) D_used = this->factory_D()->create();
    thyra_D =
      rcp_const_cast<Thyra::MultiVectorBase<value_type> >(
        rcp_dynamic_cast<const Thyra::MultiVectorBase<value_type> >(
          dyn_cast<MultiVectorMutableThyra>(*D_used).thyra_multi_vec()
          )
        );
  }

  if(
    ( rGf || D )
    &&
    !is_null(conDirecFiniteDiffCalculator_)
    &&
    ( new_thyra_N || !DfDp_is_const() )
    )
  {
    if(out.get() && trace)
      *out << "\nComputing thyra_N using directional finite differences ...\n";
    // !
    if(thyra_c.get())
      model_outArgs.set_f(rcp_const_cast<Thyra::VectorBase<value_type> >(thyra_c)); // Okay, will not be changed!
    typedef Thyra::DirectionalFiniteDiffCalculator<value_type>::SelectedDerivatives SelectedDerivatives; 
    MEB::OutArgs<value_type>
      model_fdOutArgs = conDirecFiniteDiffCalculator_->createOutArgs(
        *model_,SelectedDerivatives().supports(MEB::OUT_ARG_DfDp,0)
        );
    model_fdOutArgs.set_DfDp(p_idx_,DerivMV(thyra_N_.assert_not_null(),MEB::DERIV_MV_BY_COL));
    conDirecFiniteDiffCalculator_->calcDerivatives(
      *model_
      ,model_inArgs // The base point to compute the derivatives at
      ,model_outArgs // The base function values
      ,model_fdOutArgs // The derivatives that we want to compute
      );
  }
  // Perform solve
  if( ( D || rGf ) && py ) {
    // Solve for py and D together
    if(out.get() && trace)
      *out << "\nSolving C*[py,D] = -[c,N] simultaneously ...\n";
    const int nind = thyra_N_->domain()->dim();
    RCP<Thyra::MultiVectorBase<value_type> > 
      thyra_cN = Thyra::createMembers(thyra_N_->range(),nind+1);
    Thyra::assign(thyra_cN->col(0).ptr(), *thyra_c);
    Thyra::assign(thyra_cN->subView(Teuchos::Range1D(1,nind)).ptr(), *thyra_N_);
    RCP<Thyra::MultiVectorBase<value_type> > 
      thyra_pyD = Thyra::createMembers(thyra_D->range(),nind+1);
    Thyra::assign(thyra_pyD.ptr(), 0.0);
    Thyra::SolveStatus<value_type>
      solveStatus = Thyra::solve(*thyra_C_, Thyra::NOTRANS, *thyra_cN, thyra_pyD.ptr());
    if(out.get() && trace) {
      *out
        << "\nsolve status:\n";
      OSTab(out).o() << solveStatus;
    }
    Thyra::scale(-1.0, thyra_pyD.ptr());
    Thyra::assign(thyra_py.ptr(), *thyra_pyD->col(0));
    Thyra::assign(thyra_D.ptr(), *thyra_pyD->subView(Teuchos::Range1D(1,nind)));
  }
  else {
    // Solve for py or D
    if( py ) {
      if(out.get() && trace)
        *out << "\nSolving C*py = -c ...\n";
      Thyra::assign(thyra_py.ptr(), 0.0);
      Thyra::SolveStatus<value_type> solveStatus =
        Thyra::solve<value_type>(*thyra_C_, Thyra::NOTRANS, *thyra_c, thyra_py.ptr());
      if (nonnull(out) && trace) {
        *out << "\nsolve status:\n";
        OSTab(out).o() << solveStatus;
      }
      Thyra::Vt_S(thyra_py.ptr(), -1.0);
    }
    if( D || rGf ) {
      if(out.get() && trace)
        *out << "\nSolving C*D = -N ...\n";
      Thyra::assign(thyra_D.ptr(), 0.0);
      Thyra::SolveStatus<value_type>
        solveStatus = Thyra::solve(*thyra_C_, Thyra::NOTRANS, *thyra_N_, thyra_D.ptr());
      if(out.get() && trace) {
        *out
          << "\nsolve status:\n";
        OSTab(out).o() << solveStatus;
      }
      Thyra::scale(-1.0, thyra_D.ptr());
    }
  }
  if(thyra_py.get()) {
    free_thyra_vector(*space_c,*c,&thyra_c);
    commit_thyra_vector(*space_xD,py,&thyra_py);
  }
  // Compute reduced gradient
  if(rGf) {
    // rGf = D' * Gf_xD + Gf_xI
    const Range1D
      var_dep = basis_sys_->var_dep(),
      var_indep = basis_sys_->var_indep();
    if(objDirecFiniteDiffCalculator_.get()) {
      if(out.get() && trace)
        *out << "\nComputing rGf using directional finite differences ...\n";
      //
      // Compute each component of the reduced gradient as
      //
      // rGf(i) = fd_product( model.g, [ D(:,i); e(i) ] )
      //
      AbstractLinAlgPack::VectorDenseMutableEncap d_rGf(*rGf);
      TVecPtr e_i = createMember(model_->get_p_space(p_idx_));
      Thyra::assign(e_i.ptr(), 0.0);
      MEB::InArgs<value_type> dir = model_->createInArgs();
      dir.set_p(p_idx_,e_i);
      MEB::OutArgs<value_type> bfunc = model_->createOutArgs();
      if(f) {
        bfunc.set_g(g_idx_,createMember(model_->get_g_space(g_idx_)));
        Thyra::set_ele(0, *f, bfunc.get_g(g_idx_).ptr());
      }
      MEB::OutArgs<value_type> var = model_->createOutArgs();
      var.set_g(g_idx_,createMember(model_->get_g_space(g_idx_)));
      const int np = var_indep.size();
      for( int i = 0; i < np; ++i ) {
        set_ele(i, 1.0, e_i.ptr());
        dir.set_x(thyra_D->col(i));
        objDirecFiniteDiffCalculator_->calcVariations(
          *model_,model_inArgs,dir,bfunc,var
          );
        dir.set_x(Teuchos::null);
        Thyra::set_ele(i, 0.0, e_i.ptr());
        d_rGf()[i] = Thyra::get_ele(*var.get_g(g_idx_),0);
      }
    }
    else {
      LinAlgOpPack::V_MtV( rGf, *D_used, BLAS_Cpp::trans, *Gf->sub_view(var_dep) );
      LinAlgOpPack::Vp_V( rGf, *Gf->sub_view(var_indep) );
      // ToDo: Just compute the operators associated with Gf and not Gf directly!
    }
  }
  // * ToDo: Add specialized algorithm for computing D using an inexact Jacobian
  // * ToDo: Add in logic for inexact solves
  if(out.get() && trace)
    *out << "\nLeaving MoochoPack::NLPDirectThyraModelEvaluator::calc_point(...) ...\n";
}

void NLPDirectThyraModelEvaluator::calc_semi_newton_step(
  const Vector &x
  ,VectorMutable *c
  ,bool recalc_c
  ,VectorMutable *py
  ) const
{
  if(recalc_c) {
    //
    // Recompute c
    //
    TEUCHOS_TEST_FOR_EXCEPT(true);
  }
  // Compute py = - inv(C)*c
  TEUCHOS_TEST_FOR_EXCEPT(true);
}

}	// end namespace NLPInterfacePack
