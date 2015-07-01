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

#include <stdlib.h>

#include "AbstractLinAlgPack_SparseCOOReadMatrix.hpp"

// Throw an exception if the char is not ':'
namespace {
inline void assert_sep_char(char c) {
  if(c != ':')
    throw AbstractLinAlgPack::InputException("Sparse COO matrix input stream error:  The seperator between the element, row indice and column indice must be a \':\'");
}
inline void assert_eof(std::istream& istrm) {
  if(istrm.eof())
    throw AbstractLinAlgPack::InputException("Sparse COO matrix input stream error:  Premature end to the input file.");
}
}

void AbstractLinAlgPack::read_coo_into_valarrays(std::istream& istrm, size_type& m, size_type& n, size_type& nz
  , std::valarray<value_type>& a, std::valarray<indice_type>& ivect
  , std::valarray<indice_type>& jvect)
{
  // read in dimensions and resize
  istrm >> m;		assert_eof(istrm);
  istrm >> n;		assert_eof(istrm);
  istrm >> nz;
  a.resize(nz);
  ivect.resize(nz);
  jvect.resize(nz);

  // Read in the non-zero elements
  value_type	*p_a =			&a[0],
        *p_a_last =		p_a + nz;
  indice_type	*p_ivect =		&ivect[0],
        *p_jvect =		&jvect[0];

  for(; p_a != p_a_last; ++p_a, ++p_ivect, ++p_jvect) {
    const int bs = 50;
    char num[bs];
    char c;
    assert_eof(istrm);
    istrm.get(num, bs-1, ':');  assert_eof(istrm);  *p_a = ::atof(num);	// Read in ak 
    istrm.get(c); assert_eof(istrm); assert_sep_char(c);				// Read in ':'
    istrm.get(num, bs-1, ':'); assert_eof(istrm); *p_ivect = ::atoi(num);// Read in ik 
    istrm.get(c); assert_eof(istrm); assert_sep_char(c);				// Read in ':'
    istrm >> *p_jvect;										// Read in jk
  }
}

#endif // 0
