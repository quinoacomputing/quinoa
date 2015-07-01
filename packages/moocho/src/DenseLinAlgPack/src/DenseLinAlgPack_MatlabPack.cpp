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

#include <limits>
#include <iomanip>
#include <ostream>

#include "DenseLinAlgPack_MatlabPack.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"

namespace DenseLinAlgPack {

std::ostream& MatlabPack::out( std::ostream& o, const char* name, const DVectorSlice& vs
  , BLAS_Cpp::Transp trans )
{
  int p = o.precision();
  o.precision( std::numeric_limits<value_type>::digits10 + 3 );
  try {
    o << name << " =  [ ";
    for( DVectorSlice::const_iterator itr = vs.begin(); itr != vs.end(); ++itr )
      o << ' ' << *itr << ";";
    o << "];" << (trans == BLAS_Cpp::no_trans ? ' ' : '\'' ) << std::endl;
  }
  catch(...) {
    o.precision(p);
    throw;
  }
  o.precision(p);
  return o;
}

std::ostream& MatlabPack::out( std::ostream& o, const char* name, const DMatrixSlice& gms
  , BLAS_Cpp::Transp trans )
{
  int p = o.precision();
  o.precision( std::numeric_limits<value_type>::digits10 + 3 );
  try {
    o << name << " =  [\n";
    for( size_type i = 1; i <= gms.rows(); ++i ) {
      const DVectorSlice vs = gms.row(i);
      for( DVectorSlice::const_iterator itr = vs.begin(); itr != vs.end(); ++itr )
        o << *itr << ", ";
      o << ";\n";
    }
    o << "];" << (trans == BLAS_Cpp::no_trans ? ' ' : '\'' ) << std::endl;
  }
  catch(...) {
    o.precision(p);
    throw;
  }
  o.precision(p);
  return o;
}

}	// end namespace DenseLinAlgPack
