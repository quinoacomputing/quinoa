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

#ifndef COO_MATRIX_TMPL_OUT_FUNC_DECL_H
#define COO_MATRIX_TMPL_OUT_FUNC_DECL_H

#include "SparseLinAlgPackIOBasic.hpp"

namespace AbstractLinAlgPack {

/** \brief @name COO matrix output stream function.
  *
  * This is a functions that are used to output a templated COO matrix object
  * to a char based output stream.  The COOMatrixTemplateInterface specification
  * is used.
  *
  * The output format is diferent depending on whether the
  * bits #SparseLinAlgPackIO::ignore_dim_bit#,  #LinAlgPackIO::no_insert_newlines_bit#
  * , and #SparseLinAlgPackIO::ignore_nz_bit# are set.
  * The default output format is:
  *
  *	#rows  cols  nz#\\
  * #	val1:i1:j1	val2:i2:j2 ... valnz:inz:jnz#\\
  *
  * Each of the elements (val:i:j) are are put into columns according to the width set in the 
  * output stream #os# and other formating commands when it is called.  Even if the set
  * width is 0 or less than the number of char's for the element a space ' ' will be inserted
  * between them.  The elements are formated according to the format in the stream #os#.
  *
  * If #exta_flags & SparseLinAlgPackIO::ignore_dim_bit# == 0# then #rows# and #cols# will not
  * be output.
  *
  * If #exta_flags & SparseLinAlgPackIO::ignore_nz_bit# == 0# then #nz# will not
  * be output.
  *
  * If #exta_flags & LinAlgPackIO::no_insert_newlines_bit == 0# then a newline charactor
  * will not be inserted after the last element, otherwise (by default) it will be.
  *
  * If #coom.nz() == 0# then no elements will be output.
  *
  * If any of the output operations fails then a #std::ios_base::failure# exception is thrown. 
  */
template <class T_COOM>
std::ostream& output_COOM(std::ostream& os, const T_COOM& coom
  , SparseLinAlgPackIO::fmtflags extra_flags);

}	// end namespace AbstractLinAlgPack

#endif	// COO_MATRIX_TMPL_OUT_FUNC_DECL_H
