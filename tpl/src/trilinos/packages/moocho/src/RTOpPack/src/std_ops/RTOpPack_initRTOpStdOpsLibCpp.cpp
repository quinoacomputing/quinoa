#if 0

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

// ////////////////////////////////////////////////
// initRTOpStdOpsLibCpp.cpp

#include "RTOpPack_initRTOpStdOpsLibCpp.hpp"
#include "RTOpPack_RTOpCPostMod.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"
#include "RTOp_ROp_combined_nu_comp_err.h"
#include "RTOp_ROp_max_rel_step.h"
#include "RTOp_TOp_ele_wise_sqrt.h"
#include "RTOp_ROp_comp_err_with_mu.h"
#include "RTOp_ROp_max_step.h"
#include "RTOp_TOp_force_in_bounds.h"
#include "RTOp_ROp_dot_prod.h"
#include "RTOp_ROp_norms.h"
#include "RTOp_TOp_inv_of_difference.h"
#include "RTOp_ROp_find_nan_inf.h"
#include "RTOp_ROp_num_bounded.h"
#include "RTOp_TOp_max_abs_vec_scalar.h"
#include "RTOp_ROp_fraction_to_boundary.h"
#include "RTOp_ROp_num_nonzeros.h"
#include "RTOp_TOp_max_vec_scalar.h"
#include "RTOp_ROp_fraction_to_zero_boundary.h"
#include "RTOp_ROp_sum.h"
#include "RTOp_ROp_sum_abs.h"
#include "RTOp_TOp_multiplier_step.h"
#include "RTOp_ROp_get_ele.h"
#include "RTOp_TOp_add_scalar.h"
#include "RTOp_TOp_random_vector.h"
#include "RTOp_ROp_get_sub_vector.h"
#include "RTOp_TOp_assign_scalar.h"
#include "RTOp_TOp_scale_vector.h"
#include "RTOp_ROp_log_bound_barrier.h"
#include "RTOp_TOp_assign_vectors.h"
#include "RTOp_TOp_set_ele.h"
#include "RTOp_ROp_max_abs_ele.h"
#include "RTOp_TOp_axpy.h"
#include "RTOp_TOp_set_sub_vector.h"
#include "RTOp_ROp_max.h"
#include "RTOp_TOp_correct_multipliers.h"
#include "RTOp_TOp_sign.h"
#include "RTOp_ROp_max_inequ_viol.h"
#include "RTOp_TOp_ele_wise_divide.h"
#include "RTOp_ROp_max_near_feas_step.h"
#include "RTOp_TOp_ele_wise_prod.h"

namespace {

void add_op_factory(RTOpPack::RTOpServer<RTOp_value_type> *op_server, const RTOp_RTOp_vtbl_t &op_vtbl )
{
  namespace mmp = MemMngPack;
  namespace rtop = RTOpPack;
  typedef Teuchos::AbstractFactoryStd<rtop::RTOp,rtop::RTOpC,rtop::RTOpCPostMod>  OpFactory_t;
  op_server->add_op_factory(
    Teuchos::rcp(
      new OpFactory_t(
        rtop::RTOpCPostMod(&op_vtbl)
        ,Teuchos::AllocatorNew<rtop::RTOpC>() // RAB: 2003/10/03: MipsPro needs this *$&%^#&!!!
        )
      )
    );
}

} // namespace

void RTOpPack::initRTOpStdOpsLibCpp( RTOpPack::RTOpServer<RTOp_value_type> *op_server )
{
  namespace mmp = MemMngPack;
  typedef Teuchos::AbstractFactoryStd<RTOp,RTOpC,RTOpCPostMod>  OpFactory_t;
  
  add_op_factory(op_server,RTOp_ROp_combined_nu_comp_err_vtbl);
  add_op_factory(op_server,RTOp_ROp_combined_nu_comp_err_one_only_vtbl);
  add_op_factory(op_server,RTOp_ROp_comp_err_with_mu_vtbl);
  add_op_factory(op_server,RTOp_ROp_dot_prod_vtbl);
  add_op_factory(op_server,RTOp_ROp_find_nan_inf_vtbl);
  add_op_factory(op_server,RTOp_ROp_fraction_to_boundary_vtbl);
  add_op_factory(op_server,RTOp_ROp_fraction_to_zero_boundary_vtbl);
  add_op_factory(op_server,RTOp_ROp_get_ele_vtbl);
  add_op_factory(op_server,RTOp_ROp_get_sub_vector_vtbl);
  add_op_factory(op_server,RTOp_ROp_log_bound_barrier_vtbl);
  add_op_factory(op_server,RTOp_ROp_max_abs_ele_vtbl);
  add_op_factory(op_server,RTOp_ROp_max_vtbl);
  add_op_factory(op_server,RTOp_ROp_max_inequ_viol_vtbl);
  add_op_factory(op_server,RTOp_ROp_max_near_feas_step_vtbl);
  add_op_factory(op_server,RTOp_ROp_max_rel_step_vtbl);
  add_op_factory(op_server,RTOp_ROp_max_step_vtbl);
  add_op_factory(op_server,RTOp_ROp_norm_1_vtbl);
  add_op_factory(op_server,RTOp_ROp_norm_2_vtbl);
  add_op_factory(op_server,RTOp_ROp_norm_inf_vtbl);
  add_op_factory(op_server,RTOp_ROp_num_bounded_vtbl);
  add_op_factory(op_server,RTOp_ROp_num_nonzeros_vtbl);
  add_op_factory(op_server,RTOp_ROp_sum_vtbl);
  add_op_factory(op_server,RTOp_ROp_sum_abs_vtbl);
  add_op_factory(op_server,RTOp_TOp_add_scalar_vtbl);
  add_op_factory(op_server,RTOp_TOp_assign_scalar_vtbl);
  add_op_factory(op_server,RTOp_TOp_assign_vectors_vtbl);
  add_op_factory(op_server,RTOp_TOp_axpy_vtbl);
  add_op_factory(op_server,RTOp_TOp_Correct_Multipliers_vtbl);
  add_op_factory(op_server,RTOp_TOp_ele_wise_divide_vtbl);
  add_op_factory(op_server,RTOp_TOp_ele_wise_prod_vtbl);
  add_op_factory(op_server,RTOp_TOp_ele_wise_sqrt_vtbl);
  add_op_factory(op_server,RTOp_TOp_force_in_bounds_vtbl);
  add_op_factory(op_server,RTOp_TOp_force_in_bounds_buffer_vtbl);
  add_op_factory(op_server,RTOp_TOp_inv_of_difference_vtbl);
  add_op_factory(op_server,RTOp_TOp_max_abs_vec_scalar_vtbl);
  add_op_factory(op_server,RTOp_TOp_max_vec_scalar_vtbl);
  add_op_factory(op_server,RTOp_TOp_multiplier_step_vtbl);
  add_op_factory(op_server,RTOp_TOp_random_vector_vtbl);
  add_op_factory(op_server,RTOp_TOp_scale_vector_vtbl);
  add_op_factory(op_server,RTOp_TOp_set_ele_vtbl);
  add_op_factory(op_server,RTOp_TOp_set_sub_vector_vtbl);
  add_op_factory(op_server,RTOp_TOp_sign_vtbl);
}

#endif // 0
