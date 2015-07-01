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

#include <math.h>
#include <limits>
#include <ostream>

#include "AbstractLinAlgPack_BFGS_helpers.hpp"
#include "AbstractLinAlgPack_Vector.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"

bool AbstractLinAlgPack::BFGS_sTy_suff_p_d(
  const Vector    &s
  ,const Vector   &y
  ,const value_type     *sTy_in
  ,std::ostream         *out
  ,const char           func_name[]
  )
{
  const value_type
    sTy          = sTy_in ? *sTy_in : AbstractLinAlgPack::dot(s,y),
    nrm_s        = s.norm_2(),
    nrm_y        = y.norm_2(),
    sqrt_macheps = ::sqrt(std::numeric_limits<value_type>::epsilon()),
    min_sTy      = sqrt_macheps * nrm_s * nrm_y;
  // Skip update if: s'*y < sqrt(macheps)*||s||2*||y||2 (Dennis and Schnabel, A9.4.2)
  const bool
    sufficiently_p_d = sTy > min_sTy;
  if( !sufficiently_p_d && out ) {
    if(func_name)
      *out << func_name << " : ";
    *out
      << "Error, s'*y = " << sTy << " < sqrt(mach_eps) * ||s||2 * ||y||2 = "
      << sqrt_macheps << " * " << nrm_s << " * " << nrm_y << " = " << min_sTy
      << "\nTherefore the BFGS update is illdefined!\n";
  }
  return sufficiently_p_d;
}
