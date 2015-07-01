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

#ifndef VECTOR_IN_FUNC_H
#define VECTOR_IN_FUNC_H

#include "DenseLinAlgPack_IOBasic.hpp"

namespace DenseLinAlgPack {

/* * @name DVector/DVectorSlice input stream functions.
  *
  * These are functions that are used to read a DVector
  * or DVectorSlice object from a formated input stream.
  *
  * The input format is diferent depending on the on whether the
  * bit #LinAlgPackIO::ignore_dim_bit# is set.  If #exta_flags & LinAlgPackIO::ignore_dim_bit != 0#
  * then the input format is:
  *
  * Case 1\\
  *	#n#\\
  * #v(1) v(2) v(3) ... v(n)#\\
  *
  * If #exta_flags & LinAlgPackIO::ignore_dim_bit == 0# then the input format is:
  *
  * Case 2\\
  * #v(1) v(2) v(3) ... v(v.size())#\\
  *
  * The numbers of the input must be seperated by white space and be valid
  * C numeric constants.  For example, the input format for the vector {1.1, 2.2, 3.3}
  * for case 1 is:
  *
  * #3#\\
  * #1.1	2.2		3.3#\\
  *
  * And for case 2 is:
  *
  * #1.1	2.2		3.3#\\
  *
  * It is permisible for the dimension #n# in case 1 to be 0.  In this case there will be
  * no elements.  So to input an empty vector you would use:
  *
  * #0#\\
  *
  * If any of the input operations fails then a LinAlgPackIO::InputException exception
  * is thrown.  In other words if #is.fail()# or #is.eof()# is true
  * before all of the elements have been read in then the exception is thrown.
  * Also if the stream becomes corrupted (#is.bad() == true#) then a #std::ios_base::failure#
  * exception is thrown. 
  */
// @{

/** \brief . */
/* * DVector input stream function.
  *
  * Inputs a DVector object from an input stream.  If #exta_flags & LinAlgPackIO::ignore_dim_bit != 0#
  * then #v# is resized to #n# given in the file.  If #exta_flags & LinAlgPackIO::ignore_dim_bit == 0#
  * then the number of elements read in depends on the current size of #v#.
  *
  */
std::istream& input(std::istream& is, DVector* v, LinAlgPackIO::fmtflags extra_flags);

/** \brief . */
/* * DVectorSlice input stream function.
  *
  * Inputs a DVectorSlice object from an input stream.  If #exta_flags & LinAlgPackIO::ignore_dim_bit != 0#
  * then the size (!= 0) of #vs# is compared to the #n# given in the file and if they are not equal
  * then a #LinAlgPackIO::InputException# is thrown.  If #vs# is unsized then it is resized to #n#.
  * If #exta_flags & LinAlgPackIO::ignore_dim_bit == 0# then the number of elements read in depends
  * on the current size of #vs#.
  */
std::istream& input(std::istream& is, DVectorSlice* vs, LinAlgPackIO::fmtflags extra_flags);

// @}

}	// end namespace DenseLinAlgPack

#endif	// VECTOR_IN_FUNC_H
