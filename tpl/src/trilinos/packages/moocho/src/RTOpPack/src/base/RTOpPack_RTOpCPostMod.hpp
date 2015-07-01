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

// /////////////////////////////////////////////////////////////////////////
// RTOpPack_RTOpCPostMod.hpp

#ifndef RTOPPACK_RTOP_C_POST_MOD_HPP
#define RTOPPACK_RTOP_C_POST_MOD_HPP

#include "RTOpPack_RTOpC.hpp"

namespace RTOpPack {

class RTOpCPostMod {
public:

  /** \brief . */
  RTOpCPostMod( const RTOp_RTOp_vtbl_t *vtbl ) : vtbl_(vtbl)
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(
        !(vtbl && vtbl->obj_data_vtbl && vtbl->obj_data_vtbl->obj_create)
        ,std::logic_error, "Error!"
        );
#endif			
    }
  /** \brief . */
  void initialize(RTOpC *op) const
    {
      op->op().vtbl = vtbl_;
      op->op().vtbl->obj_data_vtbl->obj_create(NULL,NULL,&op->op().obj_data);
    }
  
private:
  
  const RTOp_RTOp_vtbl_t *vtbl_;
  
  RTOpCPostMod(); // Not defined and not to be called.

};

} // namespace RTOpPack

#endif // RTOPPACK_RTOP_C_POST_MOD_HPP
