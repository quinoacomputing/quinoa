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
#include "MoochoPack_IpState.hpp"
#include "AbstractLinAlgPack_MatrixSymDiagStd.hpp"

// Iteration Quantity Strings
const std::string MoochoPack::barrier_parameter_name = "barrier_parameter";
const std::string MoochoPack::barrier_obj_name = "barrier_obj";
const std::string MoochoPack::grad_barrier_obj_name = "grad_barrier_obj";
const std::string MoochoPack::e_tol_name = "e_tol";
const std::string MoochoPack::comp_err_mu_name = "comp_err_mu";
const std::string MoochoPack::Vu_name = "Vu";
const std::string MoochoPack::Vl_name = "Vl";
const std::string MoochoPack::invXu_name = "invXu";
const std::string MoochoPack::invXl_name = "invXl";
const std::string MoochoPack::rHB_name = "rHB";
const std::string MoochoPack::B_name = "B";
const std::string MoochoPack::Sigma_name = "Sigma";
const std::string MoochoPack::w_sigma_name = "w_sigma";
const std::string MoochoPack::dvl_name = "dvl";
const std::string MoochoPack::dvu_name = "dvu";
const std::string MoochoPack::alpha_vl_name = "alpha_vl";
const std::string MoochoPack::alpha_vu_name = "alpha_vu";

namespace MoochoPack {

IpState::IpState(
  const decomp_sys_ptr_t& decomp_sys
  ,const vec_space_ptr_t& space_x
  ,const vec_space_ptr_t& space_c
  ,const vec_space_ptr_t& space_range
  ,const vec_space_ptr_t& space_null
  )
  :
  NLPAlgoState(decomp_sys, space_x, space_c, space_range, space_null)
  {
  }

IpState::~IpState()
  {
  }

///********** Iteration Quantities **************

STATE_SCALAR_IQ_DEF(IpState, barrier_parameter, barrier_parameter_name)

STATE_SCALAR_IQ_DEF(IpState, barrier_obj, barrier_obj_name)

STATE_VECTOR_IQ_DEF(IpState, grad_barrier_obj, grad_barrier_obj_name, get_space_x(), VST_SPACE_X)

STATE_SCALAR_IQ_DEF(IpState, e_tol, e_tol_name)

STATE_SCALAR_IQ_DEF(IpState, comp_err_mu, comp_err_mu_name)

STATE_IQ_DEF(IpState, MatrixSymDiagStd, Vu, Vu_name)

STATE_IQ_DEF(IpState, MatrixSymDiagStd, Vl, Vl_name)

STATE_IQ_DEF(IpState, MatrixSymDiagStd, invXu, invXu_name)

STATE_IQ_DEF(IpState, MatrixSymDiagStd, invXl, invXl_name)

STATE_IQ_DEF(IpState, MatrixSymOp, rHB, rHB_name)

STATE_IQ_DEF(IpState, MatrixSymOp, B, B_name)

STATE_IQ_DEF(IpState, MatrixSymDiagStd, Sigma, Sigma_name)

STATE_VECTOR_IQ_DEF(IpState, w_sigma, w_sigma_name, get_space_null(), VST_SPACE_NULL )  
STATE_VECTOR_IQ_DEF(IpState, dvl, dvl_name, get_space_x(), VST_SPACE_X)
STATE_VECTOR_IQ_DEF(IpState, dvu, dvu_name, get_space_x(), VST_SPACE_X)

STATE_SCALAR_IQ_DEF(IpState, alpha_vl, alpha_vl_name)
STATE_SCALAR_IQ_DEF(IpState, alpha_vu, alpha_vu_name)
} // end namespace MoochoPack
