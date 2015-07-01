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
// This is a function that reads a sparse matrix in form an imput stream.
//

#ifndef SPARSECOOREADMATRIX_H
#define SPARSECOOREADMATRIX_H

#include <istream>
#include <valarray>

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

/** \brief Read in a Coordinate Matrix from a C++ input stream and store it in valarrays.
  *
  * The format for the imput is:
  *
  * #m  n  nz#\\
  * #a1:i1:j1  a2:i2:j2 .... anz:inz:jnz#\\
  *
  * In the above format, each non-zero element is given as a three item pair:
  * value of the non-zero element, row indice (1-based) of the non-zero element,
  * and the column indice (1-based) of the non-zero element.  There must
  * be no spaces between the numbers and the \':\' charachter and there must
  * be at least one whitespace character between elements.
  *
  * The vectors #a#, #ivect#, and #jvect# are resized to #nz#.
  * If any problem is found in the input an InputException will be throw that
  * will include a descriptive message about the error.
  *
  * @param	m		number of rows of the sparse matrix
  * @param	n		number of columns of sparse matrix
  * @param	nz		number of non-zero elements of the sparse matrix
  * @param	a		vector holding the non-zero elements
  * @param	ivect	vector holding the row indices (1-based)
  * @param	jvect	vector holding the column indices (1-based)
  */
void read_coo_into_valarrays(std::istream& istrm, size_type& m, size_type& n, size_type& nz
  , std::valarray<value_type>& a, std::valarray<indice_type>& ivect
  , std::valarray<indice_type>& jvect);

}	// end namespace AbstractLinAlgPack 

#endif // SPARSECOOREADMATRIX_H
