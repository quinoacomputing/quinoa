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

//
// Note: I am using a BasisSystem here just to make it easy to generate
// var_dep() and var_indep() and for not much else.
//

#include "NLPInterfacePack_NLPThyraModelEvaluatorBase.hpp"
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
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace NLPInterfacePack {
  
// Overridden public members from NLP

void NLPThyraModelEvaluatorBase::initialize(bool test_setup)
{
  x_guess_bounds_updated_ = false;
  updateInitialGuessAndBounds();
  if(initialized_) {
    NLPObjGrad::initialize(test_setup);
    return;
  }
  //TEUCHOS_TEST_FOR_EXCEPT(true); // Todo: push the variables in bounds!
  num_bounded_x_ = AbstractLinAlgPack::num_bounded(*xl_,*xu_,NLP::infinite_bound());
  NLPObjGrad::initialize(test_setup);
  initialized_ = true;
}

bool NLPThyraModelEvaluatorBase::is_initialized() const
{
  return initialized_;
}

NLP::vec_space_ptr_t
NLPThyraModelEvaluatorBase::space_x() const
{
  return space_x_;
}

NLP::vec_space_ptr_t
NLPThyraModelEvaluatorBase::space_c() const
{
  return space_c_;
}

size_type NLPThyraModelEvaluatorBase::num_bounded_x() const
{
  return num_bounded_x_;
}

void NLPThyraModelEvaluatorBase::force_xinit_in_bounds(bool force_xinit_in_bounds)
{
  force_xinit_in_bounds_ = force_xinit_in_bounds;
}

bool NLPThyraModelEvaluatorBase::force_xinit_in_bounds() const
{
  return force_xinit_in_bounds_;
}

const Vector& NLPThyraModelEvaluatorBase::xinit() const
{
  updateInitialGuessAndBounds();
  return *xinit_;
}

const Vector& NLPThyraModelEvaluatorBase::xl() const
{
  updateInitialGuessAndBounds();
  return *xl_;
}

const Vector& NLPThyraModelEvaluatorBase::xu() const
{
  updateInitialGuessAndBounds();
  return *xu_;
}

value_type NLPThyraModelEvaluatorBase::max_var_bounds_viol() const
{
  return 1e-5; // I have no idea?
}

void NLPThyraModelEvaluatorBase::set_f(value_type* f)
{
  NLP::set_f(f);
  f_updated_ = false;
}

void NLPThyraModelEvaluatorBase::set_c(VectorMutable* c)
{
  NLP::set_c(c);
  c_updated_ = false;
}

void NLPThyraModelEvaluatorBase::unset_quantities()
{
  NLP::unset_quantities();
}

void NLPThyraModelEvaluatorBase::scale_f( value_type scale_f )
{
  obj_scale_ = scale_f;
}

value_type NLPThyraModelEvaluatorBase::scale_f() const
{
  return obj_scale_;
}

void NLPThyraModelEvaluatorBase::report_final_solution(
  const Vector&    x
  ,const Vector*   lambda
  ,const Vector*   nu
  ,bool            optimal
  )
{
  using Teuchos::dyn_cast;
  using Teuchos::RCP;
  typedef Thyra::ModelEvaluatorBase MEB;
  using AbstractLinAlgPack::VectorMutableThyra;
  MEB::InArgs<value_type> model_finalPoint = model_->createInArgs();
  if( basis_sys_.get() ) {
    const Range1D
      var_dep   = basis_sys_->var_dep(),
      var_indep = basis_sys_->var_indep();
    RCP<const Vector> xD = x.sub_view(var_dep), xI;
    if(p_idx_>=0) xI = x.sub_view(var_indep);
    model_finalPoint.set_x(dyn_cast<const VectorMutableThyra>(*xD).thyra_vec().assert_not_null());
    if(p_idx_ >= 0)
      model_finalPoint.set_p(p_idx_,dyn_cast<const VectorMutableThyra>(*xI).thyra_vec().assert_not_null());
    else if( model_finalPoint.Np() >= 1 )
      model_finalPoint.set_p(0,model_->getNominalValues().get_p(0)); // Assume we will find in it zero!
  }
  else { // no dependent vars
    TEUCHOS_TEST_FOR_EXCEPT(p_idx_<0);
    model_finalPoint.set_p(p_idx_,dyn_cast<const VectorMutableThyra>(x).thyra_vec().assert_not_null());
  }
  model_->reportFinalPoint(model_finalPoint,optimal);
}

// Overridden public members from NLPObjGrad

void NLPThyraModelEvaluatorBase::set_Gf(VectorMutable* Gf)
{
  NLPObjGrad::set_Gf(Gf);
  Gf_updated_ = false;
}

// Overridden protected members from NLP

void NLPThyraModelEvaluatorBase::imp_calc_f(
  const Vector& x, bool newx
  ,const ZeroOrderInfo& zero_order_info
  ) const
{
  evalModel(x,newx,&zero_order_info,NULL);
}

void NLPThyraModelEvaluatorBase::imp_calc_c(
  const Vector& x, bool newx
  ,const ZeroOrderInfo& zero_order_info
  ) const
{
  evalModel(x,newx,&zero_order_info,NULL);
}

// Overridden protected members from NLPObjGrad

void NLPThyraModelEvaluatorBase::imp_calc_Gf(
  const Vector& x, bool newx
  ,const ObjGradInfo& obj_grad_info
  ) const
{
  evalModel(x,newx,NULL,&obj_grad_info);
}

// Protected functions to be used by subclasses

NLPThyraModelEvaluatorBase::NLPThyraModelEvaluatorBase()
  :showModelEvaluatorTrace_(false),initialized_(false)
  ,obj_scale_(1.0),has_bounds_(false)
  ,force_xinit_in_bounds_(true),num_bounded_x_(0)
  ,x_guess_bounds_updated_(false)
{}

void NLPThyraModelEvaluatorBase::initializeBase(
  const Teuchos::RCP<Thyra::ModelEvaluator<value_type> >  &model,
  const int p_idx,
  const int g_idx
  )
{

  using Teuchos::dyn_cast;
  using AbstractLinAlgPack::VectorSpaceThyra;
  using AbstractLinAlgPack::VectorMutableThyra;
  using AbstractLinAlgPack::MatrixOpNonsingThyra;
  typedef ::Thyra::ModelEvaluatorBase MEB;

  initialized_ = false;
  x_guess_bounds_updated_ = false;
  model_g_updated_ = model_Dg_updated_ = f_updated_ = c_updated_ = Gf_updated_ = Gc_updated_ = false;

  const char msg_err[] = "NLPThyraModelEvaluatorBase::initialize(...): Errror!";
  Thyra::ModelEvaluatorBase::OutArgs<value_type> model_outArgs = model->createOutArgs();
  TEUCHOS_TEST_FOR_EXCEPTION( model.get() == NULL, std::invalid_argument, msg_err );
  TEUCHOS_TEST_FOR_EXCEPTION( p_idx >= 0 && ( p_idx > model_outArgs.Np()-1 ), std::invalid_argument, msg_err );
  TEUCHOS_TEST_FOR_EXCEPTION( g_idx >= 0 && ( g_idx > model_outArgs.Ng()-1 ), std::invalid_argument, msg_err );
  //
  Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> > model_space_x(model->get_x_space());
  const bool no_model_x = (model_space_x.get() == NULL);
  Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> > model_space_f(model->get_f_space());
  const bool no_model_f = (model_space_f.get() == NULL);
  //
  if( !no_model_f ) {
    TEUCHOS_TEST_FOR_EXCEPTION( !model_outArgs.supports(MEB::OUT_ARG_W), std::invalid_argument, msg_err );
    MEB::DerivativeProperties model_W_properties = model_outArgs.get_W_properties();
    TEUCHOS_TEST_FOR_EXCEPTION( model_W_properties.supportsAdjoint==false, std::invalid_argument, msg_err );
    TEUCHOS_TEST_FOR_EXCEPTION( model_W_properties.rank==MEB::DERIV_RANK_DEFICIENT, std::invalid_argument, msg_err );
    /*
    if(p_idx >= 0 ) {
      TEUCHOS_TEST_FOR_EXCEPTION( model_outArgs.supports(MEB::OUT_ARG_DfDp,p_idx).none(), std::invalid_argument, msg_err );
      if(g_idx >= 0) {
        TEUCHOS_TEST_FOR_EXCEPTION( model_outArgs.supports(MEB::OUT_ARG_DgDp,g_idx,p_idx).none(), std::invalid_argument, msg_err );
      }
    }
    if(g_idx >= 0) {
      TEUCHOS_TEST_FOR_EXCEPTION( model_outArgs.supports(MEB::OUT_ARG_DgDx,g_idx).none(), std::invalid_argument, msg_err );
    }
    */
  }

  model_ = model;
  p_idx_ = p_idx;
  g_idx_ = g_idx;

  if(p_idx >= 0 ) {
    DfDp_supports_op_ = model_outArgs.supports(MEB::OUT_ARG_DfDp,p_idx).supports(MEB::DERIV_LINEAR_OP);
    DfDp_supports_mv_ = model_outArgs.supports(MEB::OUT_ARG_DfDp,p_idx).supports(MEB::DERIV_MV_BY_COL);
  }
  else {
    DfDp_supports_op_ = false;
    DfDp_supports_mv_ = false;
  }

  VectorSpace::space_ptr_t space_xI;
  if(p_idx >= 0)
    space_xI = Teuchos::rcp(new VectorSpaceThyra(model_->get_p_space(p_idx)));
  VectorSpace::space_ptr_t space_xD;
  //
  if(!no_model_x) {
    space_xD = Teuchos::rcp(new VectorSpaceThyra(model_space_x));
    if (p_idx >= 0)  {
      VectorSpace::space_ptr_t spaces_xD_xI[] = { space_xD, space_xI };
      space_x_ = Teuchos::rcp(new VectorSpaceBlocked(spaces_xD_xI,2));
    }
    else {
      space_x_ = space_xD;
    }
  }
  else {
    space_x_ = space_xI;
  } 
  TEUCHOS_TEST_FOR_EXCEPT(!space_x_.get());

  if(!no_model_f)
    space_c_ = Teuchos::rcp(new VectorSpaceThyra(model_space_f));
  else
    space_c_ = Teuchos::null;

  xinit_ = space_x_->create_member();  *xinit_ = 0.0;
  xl_    = space_x_->create_member();  *xl_    = -NLP::infinite_bound();
  xu_    = space_x_->create_member();  *xu_    = +NLP::infinite_bound();

  if(!no_model_f) {

    factory_Gc_ = BasisSystemComposite::factory_Gc();

    basis_sys_ = Teuchos::rcp(
      new BasisSystemComposite(
        space_x_
        ,space_c_
        ,Teuchos::rcp(new Teuchos::AbstractFactoryStd<MatrixOpNonsing,MatrixOpNonsingThyra>())          // factory_C
        ,Teuchos::rcp(new Teuchos::AbstractFactoryStd<MatrixSymOp,MatrixSymPosDefCholFactor>())         // factory_transDtD
        ,Teuchos::rcp(new Teuchos::AbstractFactoryStd<MatrixSymOpNonsing,MatrixSymPosDefCholFactor>())  // factory_S
        )
      );

  }
  else {

    factory_Gc_ = Teuchos::null;

    basis_sys_ = Teuchos::null;

  }
  
  if(g_idx >= 0) {
    model_g_ = createMember(model_->get_g_space(g_idx));
  }

}

void NLPThyraModelEvaluatorBase::updateInitialGuessAndBounds() const
{

  using Teuchos::dyn_cast;
  using AbstractLinAlgPack::VectorSpaceThyra;
  using AbstractLinAlgPack::VectorMutableThyra;
  using AbstractLinAlgPack::MatrixOpNonsingThyra;
  typedef ::Thyra::ModelEvaluatorBase MEB;

  if (x_guess_bounds_updated_)
    return;

  Thyra::ModelEvaluatorBase::OutArgs<value_type>
    model_outArgs = model_->createOutArgs();
  Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> >
    model_space_x = model_->get_x_space();
  const bool
    no_model_x = (model_space_x.get() == NULL);
  Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> >
    model_space_f = model_->get_f_space();
  const bool
    no_model_f = (model_space_f.get() == NULL);

  Thyra::ModelEvaluatorBase::InArgs<value_type>
    model_initialGuess = model_->getNominalValues(),
    model_lowerBounds = model_->getLowerBounds(),
    model_upperBounds = model_->getUpperBounds();

  if(!no_model_x) {
    VectorSpace::vec_mut_ptr_t xinit_D = xinit_->sub_view(basis_sys_->var_dep());
    copy_from_model_x( model_initialGuess.get_x().get(), &*xinit_D );
    VectorSpace::vec_mut_ptr_t xl_D = xl_->sub_view(basis_sys_->var_dep());
    copy_from_model_x( model_lowerBounds.get_x().get(), &*xl_D );
    VectorSpace::vec_mut_ptr_t xu_D = xu_->sub_view(basis_sys_->var_dep());
    copy_from_model_x( model_upperBounds.get_x().get(), &*xu_D );
  }

  if(p_idx_ >= 0) {
    Range1D var_indep = ( basis_sys_.get() ? basis_sys_->var_indep() : Range1D() );
    VectorSpace::vec_mut_ptr_t xinit_I = xinit_->sub_view(var_indep);
    copy_from_model_p( model_initialGuess.get_p(p_idx_).get(), &*xinit_I );
    VectorSpace::vec_mut_ptr_t xl_I = xl_->sub_view(var_indep);
    copy_from_model_p( model_lowerBounds.get_p(p_idx_).get(), &*xl_I );
    VectorSpace::vec_mut_ptr_t xu_I = xu_->sub_view(var_indep);
    copy_from_model_p( model_upperBounds.get_p(p_idx_).get(), &*xu_I );
  }

  x_guess_bounds_updated_ = true;

}

void NLPThyraModelEvaluatorBase::assert_is_initialized() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !is_initialized(), NLP::UnInitialized
    ,"NLPThyraModelEvaluatorBase::assert_is_initialized() : Error, "
    "NLPThyraModelEvaluatorBase::initialize() has not been called yet."
    );
}

void NLPThyraModelEvaluatorBase::copy_from_model_x( const Thyra::VectorBase<value_type>* model_x, VectorMutable* x_D ) const
{
  if(!model_x) return;
  *x_D = AbstractLinAlgPack::VectorMutableThyra(Teuchos::rcp(const_cast<Thyra::VectorBase<value_type>*>(model_x),false));
}

void NLPThyraModelEvaluatorBase::copy_from_model_p( const Thyra::VectorBase<value_type> *model_p, VectorMutable* x_I ) const
{
  if(!model_p) return;
  *x_I = AbstractLinAlgPack::VectorMutableThyra(Teuchos::rcp(const_cast<Thyra::VectorBase<value_type>*>(model_p),false));
}

void NLPThyraModelEvaluatorBase::set_x(
  const Vector                                      &x
  ,Thyra::ModelEvaluatorBase::InArgs<value_type>    *model_inArgs_inout
  ) const
{
  using Teuchos::dyn_cast;
  using Teuchos::RCP;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using AbstractLinAlgPack::VectorMutableThyra;
  typedef Thyra::ModelEvaluatorBase MEB;
  //
  // Set the input arguments
  //
  MEB::InArgs<value_type> &model_inArgs = *model_inArgs_inout;
  if( basis_sys_.get() ) {
    const Range1D
      var_dep   = basis_sys_->var_dep(),
      var_indep = basis_sys_->var_indep();
    RCP<const Vector> xD = x.sub_view(var_dep), xI;
    if(p_idx_>=0) xI = x.sub_view(var_indep);
    model_inArgs.set_x(dyn_cast<const VectorMutableThyra>(*xD).thyra_vec().assert_not_null());
    if(p_idx_ >= 0)
      model_inArgs.set_p(p_idx_,dyn_cast<const VectorMutableThyra>(*xI).thyra_vec().assert_not_null());
  }
  else { // no dependent vars
    TEUCHOS_TEST_FOR_EXCEPT(p_idx_<0);
    model_inArgs.set_p(p_idx_,dyn_cast<const VectorMutableThyra>(x).thyra_vec().assert_not_null());
  }
}

void NLPThyraModelEvaluatorBase::preprocessBaseInOutArgs(
  const Vector                                      &x
  ,bool                                             newx
  ,const ZeroOrderInfo                              *zero_order_info
  ,const ObjGradInfo                                *obj_grad_info
  ,const NLPFirstOrder::FirstOrderInfo              *first_order_info
  ,Thyra::ModelEvaluatorBase::InArgs<value_type>    *model_inArgs_inout
  ,Thyra::ModelEvaluatorBase::OutArgs<value_type>   *model_outArgs_inout
  ,MatrixOp*                                        *Gc_out
  ,VectorMutable*                                   *Gf_out
  ,value_type*                                      *f_out
  ,VectorMutable*                                   *c_out
  ) const
{
  using Teuchos::dyn_cast;
  using Teuchos::RCP;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using AbstractLinAlgPack::VectorMutableThyra;
  using AbstractLinAlgPack::MatrixOpThyra;
  using AbstractLinAlgPack::MatrixOpNonsingThyra;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef MEB::DerivativeMultiVector<value_type> DerivMV;
  typedef MEB::Derivative<value_type> Deriv;
  //
  if(newx) model_g_updated_ = model_Dg_updated_ = f_updated_ = c_updated_ = Gf_updated_ = Gc_updated_ = false;
  //
  // Set the input arguments
  //
  MEB::InArgs<value_type> &model_inArgs = *model_inArgs_inout;
  set_x(x,&model_inArgs);
  //
  // Set the output arguments
  //
  MatrixOp            *Gc = NULL;
  VectorMutable       *Gf = NULL;
  value_type          *f  = NULL;
  VectorMutable       *c  = NULL;
  //
  if(zero_order_info) {
    f = zero_order_info->f;
    c = zero_order_info->c;
  }
  else if(obj_grad_info) {
    Gf = obj_grad_info->Gf;
    f  = obj_grad_info->f;
    c  = obj_grad_info->c;
  }
  else if(first_order_info) {
    Gc = first_order_info->Gc;
    Gf = first_order_info->Gf;
    f  = first_order_info->f;
    c  = first_order_info->c;
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPT(true); // Should never be called!
  }
  //
  MEB::OutArgs<value_type> &model_outArgs = *model_outArgs_inout;
  if( f && (g_idx_>=0) && !f_updated_ ) {
    model_outArgs.set_g(g_idx_,model_g_.assert_not_null()); // ToDo: Make more general!
  }
  if( c && !c_updated_ ) {
    Teuchos::RCP<Thyra::VectorBase<value_type> > thyra_c;
    get_thyra_vector(*space_c_,c,&thyra_c);
    model_outArgs.set_f(thyra_c.assert_not_null());
  }
  if( Gf && !Gf_updated_ ) {
    if(g_idx_>=0) {
      if(p_idx_>=0) {
        const Range1D
          var_dep   = ( basis_sys_.get() ? basis_sys_->var_dep() : Range1D::Invalid ),
          var_indep = ( basis_sys_.get() ? basis_sys_->var_indep() : Range1D() );
        if( var_dep.size() ) {
          model_outArgs.set_DgDx(
            g_idx_
            ,DerivMV(
              rcp_const_cast<Thyra::VectorBase<value_type> >(
                dyn_cast<VectorMutableThyra>(*Gf->sub_view(var_dep)).thyra_vec()
                ).assert_not_null()
              ,MEB::DERIV_TRANS_MV_BY_ROW
              )
            );
        }
        model_outArgs.set_DgDp(
          g_idx_,p_idx_
          ,DerivMV(
            rcp_const_cast<Thyra::VectorBase<value_type> >(
              dyn_cast<VectorMutableThyra>(*Gf->sub_view(var_indep)).thyra_vec()
              ).assert_not_null()
            ,MEB::DERIV_TRANS_MV_BY_ROW
            )
          );
      }
    }
  }
  if(Gc_out) *Gc_out = Gc;
  if(Gf_out) *Gf_out = Gf;
  if(f_out)  *f_out  = f;
  if(c_out)  *c_out  = c;
}

void NLPThyraModelEvaluatorBase::postprocessBaseOutArgs(
  Thyra::ModelEvaluatorBase::OutArgs<value_type>        *model_outArgs_inout
  ,VectorMutable                                        *Gf
  ,value_type                                           *f
  ,VectorMutable                                        *c
  ) const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::OutArgs<value_type> &model_outArgs = *model_outArgs_inout;
  if( f && !f_updated_ ) {
    if(g_idx_>=0) {
      *f = obj_scale_ * ::Thyra::get_ele(*model_g_,0);
    }
    else {
      *f = 0.0;
    }
    f_updated_ = true;
  }
  if( c && !c_updated_ ) {
    Teuchos::RCP<Thyra::VectorBase<value_type> >
      thyra_c = model_outArgs.get_f();
    commit_thyra_vector(*space_c_,c,&thyra_c);
    model_outArgs.set_f(Teuchos::null);
    c_updated_ = true;
  }
  if( Gf && !Gf_updated_ ) {
    if(g_idx_>=0) {
      const Range1D
        var_dep   = ( basis_sys_.get() ? basis_sys_->var_dep() : Range1D::Invalid ),
        var_indep = ( basis_sys_.get() ? basis_sys_->var_indep() : Range1D() );
      if (obj_scale_ != 1.0 )
        Vt_S( Gf, obj_scale_ );
      if(var_dep.size())
        Gf->sub_view(var_dep)->has_changed();
      if(var_indep.size())
        Gf->sub_view(var_indep)->has_changed();
    }
    else {
      *Gf = 0.0;
    }
    Gf_updated_ = true;
  }
}

// private

void NLPThyraModelEvaluatorBase::evalModel( 
  const Vector            &x
  ,bool                   newx
  ,const ZeroOrderInfo    *zero_order_info
  ,const ObjGradInfo      *obj_grad_info
  ) const
{
  typedef Thyra::ModelEvaluatorBase MEB;
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
  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out << "\nEntering MoochoPack::NLPThyraModelEvaluatorBase::evalModel(...) ...\n";
  //
  // Set the input and output arguments
  //
  MEB::InArgs<value_type>  model_inArgs  = model_->createInArgs();
  MEB::OutArgs<value_type> model_outArgs = model_->createOutArgs();
  VectorMutable       *Gf = NULL;
  value_type          *f  = NULL;
  VectorMutable       *c  = NULL;
  preprocessBaseInOutArgs(
    x,newx,zero_order_info,obj_grad_info,NULL
    ,&model_inArgs,&model_outArgs,NULL,&Gf,&f,&c
    );
  //
  // Evaluate the model
  //
  model_->evalModel(model_inArgs,model_outArgs);
  //
  // Postprocess the output arguments
  //
  postprocessBaseOutArgs(&model_outArgs,Gf,f,c);
  //
  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out << "\nLeaving MoochoPack::NLPThyraModelEvaluatorBase::evalModel(...) ...\n";
}

}	// end namespace NLPInterfacePack
