// @HEADER
// ***********************************************************************
// 
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//                Copyright (2006) Sandia Corporation
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

#include "RTOpPack_SPMD_apply_op_decl.hpp"
#include "RTOpPack_SPMD_apply_op_def.hpp"


#ifdef RTOPPACK_ENABLE_SHOW_DUMP
// Keep for backward compatibility but it is a no-op now!
bool RTOpPack::show_spmd_apply_op_dump = false;
#endif // RTOPPACK_ENABLE_SHOW_DUMP


Teuchos::RCP<Teuchos::FancyOStream>& RTOpPack::spmdApplyOpDumpOut()
{
  static Teuchos::RCP<Teuchos::FancyOStream> dumpOut;
  return dumpOut;
}


void RTOpPack::set_SPMD_apply_op_dump_out(const RCP<FancyOStream> &dumpOut)
{
  spmdApplyOpDumpOut() = dumpOut;
}


#ifdef HAVE_RTOP_EXPLICIT_INSTANTIATION


#include "Teuchos_ExplicitInstantiationHelpers.hpp"


namespace RTOpPack {


TEUCHOS_MACRO_TEMPLATE_INSTANT_SCALAR_TYPES(
  RTOPPACK_SPMD_APPLY_OP_INSTANT_SCALAR)


} // namespace RTOpPack


#endif // HAVE_TEUCHOS_EXCPLICIT_INSTANTIATION
