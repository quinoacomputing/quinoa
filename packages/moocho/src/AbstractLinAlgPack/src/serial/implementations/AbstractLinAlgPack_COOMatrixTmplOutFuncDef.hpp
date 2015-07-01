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

#ifndef COO_MATRIX_TMPL_OUT_FUNC_DEF_H
#define COO_MATRIX_TMPL_OUT_FUNC_DEF_H

#include <ostream>
#include <iomanip>

#include "AbstractLinAlgPack_COOMatrixTmplOutFuncDecl.hpp"

namespace AbstractLinAlgPack {

template <class T_COOM>
std::ostream& output_COOM(std::ostream& os, const T_COOM& coom
  , SparseLinAlgPackIO::fmtflags extra_flags)
{
  int w = os.width(0) - 1; // get the set width (minus 1 since a space is inserted)

  if(    !(extra_flags & SparseLinAlgPackIO::ignore_dim_bit)
    || !(extra_flags & SparseLinAlgPackIO::ignore_nz_bit) )
  {

    os << std::setw(0) << std::left;

    if( !(extra_flags & SparseLinAlgPackIO::ignore_dim_bit) )
      os << coom.rows() << ' ' << coom.cols() << ' ';

    if( !(extra_flags & SparseLinAlgPackIO::ignore_nz_bit) )
      os << coom.nz();

    os << std::endl << std::right;
  }
  
  if(!coom.nz()) return os;	// no elements to output
  
  typename T_COOM::difference_type
    row_offset	= coom.row_offset(),
    col_offset	= coom.col_offset();
  for(typename T_COOM::const_iterator itr = coom.begin(); itr != coom.end();++itr) {
    os	<< " " << std::setw(w) << itr->value()
      << ':' << itr->row_i() + row_offset
      << ':' << itr->col_j() + col_offset;
  }

  if( !(extra_flags & SparseLinAlgPackIO::no_insert_newlines_bit) )
    os << std::endl;

  return os;
}

}	// end namespace AbstractLinAlgPack 

#endif	// COO_MATRIX_TMPL_OUT_FUNC_DEF_H
