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

#include <ostream>
#include <iomanip>

#include "AbstractLinAlgPack_assert_print_nan_inf.hpp"
#include "AbstractLinAlgPack_Vector.hpp"
#include "RTOp_ROp_find_nan_inf.h"
#include "RTOpPack_RTOpC.hpp"
#include "check_nan_inf.h"
#include "Teuchos_Assert.hpp"

namespace {

// Find a NaN or Inf element!
static RTOpPack::RTOpC                               find_nan_inf_op;
static Teuchos::RCP<RTOpPack::ReductTarget>  find_nan_inf_targ;

class init_rtop_server_t {
public:
  init_rtop_server_t() {
    TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_ROp_find_nan_inf_construct(&find_nan_inf_op.op() ));
    find_nan_inf_targ = find_nan_inf_op.reduct_obj_create();
  }
};

init_rtop_server_t  init_rtop_server;

} // end namespace

bool AbstractLinAlgPack::assert_print_nan_inf( const value_type& val, const char name[]
  , bool throw_excpt, std::ostream* out )
{
  if( RTOp_is_nan_inf(val) ) {
    std::ostringstream omsg;
    omsg
      << "The scalar \"" << name
      << "\" = " << val << " is not a valid bounded number";
    if(out)
      *out << omsg.str() << std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(
      throw_excpt,NaNInfException
      ,"assert_print_nan_inf(...) : Error, " << omsg.str() );
    return false;
  }
  return true;
}

bool AbstractLinAlgPack::assert_print_nan_inf(
  const Vector& v, const char name[]
  ,bool throw_excpt, std::ostream* out
  )
{
  find_nan_inf_op.reduct_obj_reinit(find_nan_inf_targ.ptr());
  const Vector* vecs[1] = { &v };
  apply_op(find_nan_inf_op,1,vecs,0,NULL,&*find_nan_inf_targ);
  RTOp_ROp_find_nan_inf_reduct_obj_t
    ele =RTOp_ROp_find_nan_inf_val(find_nan_inf_op(*find_nan_inf_targ));
  if(out && ele.i) {
    *out
      << "The vector \"" << name << "\" has the first following NaN or Inf element\n"
      << name << "(" << ele.i << ") = " << ele.v0_i << std::endl;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(
    ele.i && throw_excpt, NaNInfException
    ,"assert_print_nan_inf(...) : Error, the vector named "
    << name << " has at least one element which is NaN or Inf" );

  return ele.i == 0;
}
