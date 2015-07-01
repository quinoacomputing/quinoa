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

#if !defined IP_STATE_H
#define IP_STATE_H

#include "MoochoPack_NLPAlgoState.hpp"

namespace MoochoPack {

// Iteration Quantity Strings
extern const std::string barrier_parameter_name;
extern const std::string barrier_obj_name;
extern const std::string grad_barrier_obj_name;
extern const std::string e_tol_name;
extern const std::string comp_err_mu_name;
extern const std::string Vu_name;
extern const std::string Vl_name;
extern const std::string invXu_name;
extern const std::string invXl_name;
extern const std::string rHB_name;
extern const std::string B_name;
extern const std::string Sigma_name;
extern const std::string w_sigma_name;
extern const std::string dvl_name;
extern const std::string dvu_name;
extern const std::string alpha_vl_name;
extern const std::string alpha_vu_name;


class IpState 
  : public MoochoPack::NLPAlgoState
  {

  public:
    ///********** Iteration Quantities **************

    /// mu: barrier parameter
    STATE_SCALAR_IQ_DECL(barrier_parameter)

    /// barrier_obj: objective value with 
    //   barrier term included
    STATE_SCALAR_IQ_DECL(barrier_obj)

    /// grad_barrier_obj: gradient of the objective
    //   with barrier term included
    STATE_VECTOR_IQ_DECL(grad_barrier_obj)

    /// e_tol: current error tolerance for inner loop
    STATE_SCALAR_IQ_DECL(e_tol)

    /// comp_err_mu: perturbed complementarity error for barrier sub problem
    STATE_SCALAR_IQ_DECL(comp_err_mu)

    /// Vu - diagonal matrix of upper bound multipliers
    STATE_IQ_DECL(MatrixSymDiagStd, Vu)

    /// Vl - diagonal matrix of lower bound multipliers
    STATE_IQ_DECL(MatrixSymDiagStd, Vl)

    /// invXu - (Xu)^-1 - matrix of 1/(xu-x) diagonal
    STATE_IQ_DECL(MatrixSymDiagStd, invXu)

    /// invXl - (Xl)^-1 - matrix of 1/(x-xl) diagonal
    STATE_IQ_DECL(MatrixSymDiagStd, invXl)

    /// rHB - reduced Hessian of the barrier term (Z_Sigma_Z)
    STATE_IQ_DECL(MatrixSymOp, rHB)

    /// B - overall reduced 'Hessian' (Z_W_Z+Z_Sigma_Z)
    STATE_IQ_DECL(MatrixSymOp, B)

    /// Full space Sigma (invXl*Vl-invXu*Vu)
    STATE_IQ_DECL(MatrixSymDiagStd, Sigma)

    /// w_sigma:  crossterm correction for sigma (Z' * Sigma * Y * py)
    STATE_VECTOR_IQ_DECL(w_sigma) 

    /// dvl:  Search direction for lower bound multipliers ( n x 1 )
    STATE_VECTOR_IQ_DECL(dvl)

    /// dvu:  Search direction for upper bound multipliers ( n x 1 )
    STATE_VECTOR_IQ_DECL(dvu)

    /// alpha_vl: step size for vl
    STATE_SCALAR_IQ_DECL(alpha_vl)

    /// alpha_vl: step size for vu
    STATE_SCALAR_IQ_DECL(alpha_vu)

    /** \brief Construct
     *
     * 
     */
    IpState(
      const decomp_sys_ptr_t& decomp_sys   = Teuchos::null
      ,const vec_space_ptr_t& space_x      = Teuchos::null
      ,const vec_space_ptr_t& space_c      = Teuchos::null
      ,const vec_space_ptr_t& space_range  = Teuchos::null
      ,const vec_space_ptr_t& space_null   = Teuchos::null
      );

    virtual ~IpState();

  }; // end class IpState

} // end namespace MoochoPack



#endif // if !defined IP_STATE_H
