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

#include <sstream>
#include <typeinfo>

#include "MoochoPack_NLPAlgoState.hpp"
#include "ConstrainedOptPack_MeritFuncNLP.hpp"
#include "AbstractLinAlgPack_MatrixSymOp.hpp"
#include "AbstractLinAlgPack_MatrixOpNonsing.hpp"
#include "Teuchos_dyn_cast.hpp"

#include "IterationPack_IterQuantityAccess.hpp"
#include "IterationPack_cast_iq.hpp"
#include "IterationPack_IterQuantityAccessContiguous.hpp"

// NLPAlgoState iteration quantities names

// Iteration Info
const std::string MoochoPack::num_basis_name		= "num_basis";
// NLP Problem Info 
const std::string MoochoPack::x_name				= "x";
const std::string MoochoPack::f_name				= "f";
const std::string MoochoPack::Gf_name				= "Gf";
const std::string MoochoPack::HL_name				= "HL";
const std::string MoochoPack::c_name				= "c";
const std::string MoochoPack::Gc_name				= "Gc";
// Constraint Gradient Null Space / Range Space Decomposition Info
const std::string MoochoPack::Y_name				= "Y";
const std::string MoochoPack::Z_name				= "Z";
const std::string MoochoPack::R_name				= "R";
const std::string MoochoPack::Uy_name				= "Uy";
const std::string MoochoPack::Uz_name				= "Uz";
// Search Direction Info
const std::string MoochoPack::py_name				= "py";
const std::string MoochoPack::Ypy_name				= "Ypy";
const std::string MoochoPack::pz_name				= "pz";
const std::string MoochoPack::Zpz_name				= "Zpz";
const std::string MoochoPack::d_name				= "d";
// Reduced QP Subproblem Info
const std::string MoochoPack::rGf_name				= "rGf";
const std::string MoochoPack::rHL_name				= "rHL";
const std::string MoochoPack::w_name				= "w";
const std::string MoochoPack::zeta_name			= "zeta";
const std::string MoochoPack::qp_grad_name			= "qp_grad";
const std::string MoochoPack::eta_name				= "eta";
// Global Convergence Info
const std::string MoochoPack::alpha_name			= "alpha";
const std::string MoochoPack::merit_func_nlp_name	= "merit_func_nlp";
const std::string MoochoPack::mu_name				= "mu";
const std::string MoochoPack::phi_name				= "phi";
// KKT Info
const std::string MoochoPack::opt_kkt_err_name		= "opt_kkt_err";
const std::string MoochoPack::feas_kkt_err_name	= "feas_kkt_err";
const std::string MoochoPack::comp_kkt_err_name	= "comp_kkt_err";
const std::string MoochoPack::GL_name				= "GL";
const std::string MoochoPack::rGL_name				= "rGL";
const std::string MoochoPack::lambda_name			= "lambda";
const std::string MoochoPack::nu_name				= "nu";

namespace MoochoPack {

// Constructors / initializers

void NLPAlgoState::set_space_range (const vec_space_ptr_t& space_range )
{
  space_range_ = space_range;
  update_vector_factories(VST_SPACE_RANGE,space_range);
}

void NLPAlgoState::set_space_null (const vec_space_ptr_t& space_null )
{
  space_null_ = space_null;
  update_vector_factories(VST_SPACE_NULL,space_null);
}

NLPAlgoState::NLPAlgoState(
  const decomp_sys_ptr_t& decomp_sys
  ,const vec_space_ptr_t& space_x
  ,const vec_space_ptr_t& space_c
  ,const vec_space_ptr_t& space_range
  ,const vec_space_ptr_t& space_null
  )
  :decomp_sys_(decomp_sys)
  ,space_x_(space_x)
  ,space_c_(space_c)
  ,space_range_(space_range)
  ,space_null_(space_null)
{}

// Iteration Info

STATE_INDEX_IQ_DEF(  NLPAlgoState,              num_basis, num_basis_name          )

// NLP Problem Info

STATE_VECTOR_IQ_DEF( NLPAlgoState,              x,         x_name,  get_space_x(), VST_SPACE_X  )
STATE_SCALAR_IQ_DEF( NLPAlgoState,              f,         f_name                               )
STATE_IQ_DEF(        NLPAlgoState, MatrixSymOp, HL,        HL_name                              )
STATE_VECTOR_IQ_DEF( NLPAlgoState,              Gf,        Gf_name, get_space_x(), VST_SPACE_X  )
STATE_VECTOR_IQ_DEF( NLPAlgoState,              c,         c_name,  get_space_c(), VST_SPACE_C  )
STATE_IQ_DEF(        NLPAlgoState, MatrixOp,    Gc,        Gc_name                              )

// Constraint Gradient Null Space / Range Space Decomposition Info

STATE_IQ_DEF(        NLPAlgoState, MatrixOp,        Y,  Y_name                  )
STATE_IQ_DEF(        NLPAlgoState, MatrixOp,        Z,  Z_name                  )
STATE_IQ_DEF(        NLPAlgoState, MatrixOpNonsing, R,  R_name                  )
STATE_IQ_DEF(        NLPAlgoState, MatrixOp,        Uy, Uy_name                 )
STATE_IQ_DEF(        NLPAlgoState, MatrixOp,        Uz, Uz_name                 )

// Search Direction Info

STATE_VECTOR_IQ_DEF( NLPAlgoState,                  py,  py_name,   get_space_range(), VST_SPACE_RANGE )
STATE_VECTOR_IQ_DEF( NLPAlgoState,                  Ypy, Ypy_name,  get_space_x(),     VST_SPACE_X     )
STATE_VECTOR_IQ_DEF( NLPAlgoState,                  pz,  pz_name,   get_space_null(),  VST_SPACE_NULL  )
STATE_VECTOR_IQ_DEF( NLPAlgoState,                  Zpz, Zpz_name,  get_space_x(),     VST_SPACE_X     )
STATE_VECTOR_IQ_DEF( NLPAlgoState,                  d,   d_name,    get_space_x(),     VST_SPACE_X     )

// QP Subproblem Info

STATE_VECTOR_IQ_DEF( NLPAlgoState,                  rGf,     rGf_name,      get_space_null(), VST_SPACE_NULL )
STATE_IQ_DEF(        NLPAlgoState, MatrixSymOp,     rHL,     rHL_name                                        )
STATE_VECTOR_IQ_DEF( NLPAlgoState,                  w,       w_name,        get_space_null(), VST_SPACE_NULL ) 
STATE_SCALAR_IQ_DEF( NLPAlgoState,                  zeta,    zeta_name                                       )
STATE_VECTOR_IQ_DEF( NLPAlgoState,                  qp_grad, qp_grad_name,  get_space_null(), VST_SPACE_NULL )
STATE_SCALAR_IQ_DEF( NLPAlgoState,                  eta,     eta_name                                        )

// Global Convergence Info

STATE_SCALAR_IQ_DEF( NLPAlgoState,                  alpha,          alpha_name          )
STATE_IQ_DEF(        NLPAlgoState, MeritFuncNLP,    merit_func_nlp, merit_func_nlp_name )
STATE_SCALAR_IQ_DEF( NLPAlgoState,                  mu,             mu_name             )
STATE_SCALAR_IQ_DEF( NLPAlgoState,                  phi,            phi_name            )

// KKT Info

STATE_SCALAR_IQ_DEF( NLPAlgoState,                  opt_kkt_err,    opt_kkt_err_name                                    )
STATE_SCALAR_IQ_DEF( NLPAlgoState,                  feas_kkt_err,   feas_kkt_err_name                                   )
STATE_SCALAR_IQ_DEF( NLPAlgoState,                  comp_kkt_err,   comp_kkt_err_name                                   )
STATE_VECTOR_IQ_DEF( NLPAlgoState,                  GL,             GL_name,           get_space_x(),    VST_SPACE_X    )
STATE_VECTOR_IQ_DEF( NLPAlgoState,                  rGL,            rGL_name,          get_space_null(), VST_SPACE_NULL )
STATE_VECTOR_IQ_DEF( NLPAlgoState,                  lambda,         lambda_name,       get_space_c(),    VST_SPACE_C    )
STATE_VECTOR_IQ_DEF( NLPAlgoState,                  nu,             nu_name,           get_space_x(),    VST_SPACE_X    )

// protected

void NLPAlgoState::update_iq_id(
  const std::string&                iq_name
  ,iq_id_encap*                     iq_id
  ) const
{
  namespace rcp = MemMngPack;
  if(iq_id->iq_id == DOES_NOT_EXIST)
    iq_id->iq_id = this->get_iter_quant_id(iq_name);
  TEUCHOS_TEST_FOR_EXCEPTION(
    iq_id->iq_id == DOES_NOT_EXIST, DoesNotExist
    ,"NLPAlgoState::update_iq_id(iq_name,iq_id) : Error, "
    " The iteration quantity with name \'" << iq_name <<
    "\' does not exist!" );
}

void NLPAlgoState::update_index_type_iq_id(
  const std::string&                iq_name
  ,iq_id_encap*                     iq_id
  )
{
  namespace rcp = MemMngPack;
  if(iq_id->iq_id == DOES_NOT_EXIST) {
    iq_id_type
      _iq_id = this->get_iter_quant_id(iq_name);
    if(_iq_id == DOES_NOT_EXIST) {
      iq_id->iq_id = this->set_iter_quant(
        iq_name
        ,Teuchos::rcp(
          new IterQuantityAccessContiguous<index_type>(
            1
            ,iq_name
#ifdef _MIPS_CXX
            ,Teuchos::RCP<Teuchos::AbstractFactoryStd<index_type,index_type> >(
              new Teuchos::AbstractFactoryStd<index_type,index_type>())
#endif
            )
          )
        );
    }
    else {
      iq_id->iq_id = _iq_id;
    }
  }
}

void NLPAlgoState::update_value_type_iq_id(
  const std::string&                iq_name
  ,iq_id_encap*                     iq_id
  )
{
  namespace rcp = MemMngPack;
  if(iq_id->iq_id == DOES_NOT_EXIST) {
    iq_id_type
      _iq_id = this->get_iter_quant_id(iq_name);
    if(_iq_id == DOES_NOT_EXIST) {
      iq_id->iq_id = this->set_iter_quant(
        iq_name
        ,Teuchos::rcp(
          new IterQuantityAccessContiguous<value_type>(
            1
            ,iq_name
#ifdef _MIPS_CXX
            ,Teuchos::RCP<Teuchos::AbstractFactoryStd<value_type,value_type> >(
              new Teuchos::AbstractFactoryStd<value_type,value_type>())
#endif
            )
          )
        );
    }
    else {
      iq_id->iq_id = _iq_id;
    }
  }
}

void NLPAlgoState::update_vector_iq_id(
  const std::string&                iq_name
  ,const VectorSpace::space_ptr_t&  vec_space
  ,EVecSpaceType                    vec_space_type
  ,iq_id_encap*                     iq_id
  )
{
  namespace rcp = MemMngPack;
  if(iq_id->iq_id == DOES_NOT_EXIST) {
    iq_id_type
      _iq_id = this->get_iter_quant_id(iq_name);
    if(_iq_id == DOES_NOT_EXIST) {
      iq_id->iq_id = this->set_iter_quant(
        iq_name
        ,Teuchos::rcp(
          new IterQuantityAccessContiguous<VectorMutable>(
            1
            ,iq_name
            ,vec_space
            )
          )
        );
    }
    else {
      iq_id->iq_id = _iq_id;
    }
    // Record the list of vectors for a given vector space. 
    vector_iqs_lists_[vec_space_type].push_back(iq_id->iq_id);
  }
}

// private

void NLPAlgoState::update_vector_factories(
  EVecSpaceType             vec_space_type
  ,const vec_space_ptr_t&   vec_space
  )
{
  using Teuchos::dyn_cast;
  iq_vector_list_t  &iq_vector_list = vector_iqs_lists_[vec_space_type];
  for( iq_vector_list_t::const_iterator iq_itr = iq_vector_list.begin(); iq_itr != iq_vector_list.end(); ++iq_itr )
    dyn_cast<IterQuantityAccessContiguous<VectorMutable> >(this->iter_quant(*iq_itr)).set_factory(vec_space);
}

}	// end namespace MoochoPack
