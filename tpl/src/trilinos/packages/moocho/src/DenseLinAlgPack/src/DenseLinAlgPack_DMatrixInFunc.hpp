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

#ifndef GENMATRIX_IN_FUNC_H
#define GENMATRIX_IN_FUNC_H

#include "DenseLinAlgPack_IOBasic.hpp"

namespace DenseLinAlgPack {

/* * @name DMatrix/DMatrixSlice input stream functions.
  *
  * These are functions that are used to read a DMatrix or DMatrixSlice object in from a
  * formated input stream.
  *
  * The input format is diferent depending on the on whether the bit
  * #LinAlgPackIO::ignore_dim_bit# is set.
  * If #exta_flags & LinAlgPackIO::ignore_dim_bit != 0# then the input format is:
  *
  * Case 1\\
  *	#m  n#\\
  * #gm(1,1) gm(1,2) gm(1,3) ... gm(1,n)#\\
  * #gm(2,1) gm(2,2) gm(2,3) ... gm(2,n)#\\
  * #   .       .       .           .#\\
  * #gm(m,1) gm(m,2) gm(m,3) ... gm(m,n)#\\
  *
  * If #exta_flags & LinAlgPackIO::ignore_dim_bit == 0# then the input format is:
  *
  * Case 2\\
  * #gm(1,1) gm(1,2) gm(1,3) ... gm(1,n)#\\
  * #gm(2,1) gm(2,2) gm(2,3) ... gm(2,n)#\\
  * #   .       .       .           .#\\
  * #gm(m,1) gm(m,2) gm(m,3) ... gm(m,n)#\\
  *
  * The numbers of the input must be seperated by white space (the line breaks are
  * for looks only) and be valid C language numeric constants.
  *
  * In addition, comment lines may be inserted between rows of the input matrix.
  * These comment lines take the form of Fortan comment lines in that they start on
  * new lines (after a '\n' char) with a '*' char and end at the end of a line
  * ('\n' terminated).  After the elements for a row are read in the function
  * #eat_comment_lines(is,'*');# is called.
  *
  * For example, the input format for the matrix {1.1 1.2; 2.1 2.2}
  * with comments for case 1 is:
  *
  * #2  2#\\
  * #* This is the first row#\\
  * #1.1	1.2#\\
  * #* This is the second row#\\
  * #2.1	2.1#\\
  *
  * And for case 2 is:
  *
  * #* This is the first row#\\
  * #1.1	1.2#\\
  * #* This is the second row#\\
  * #2.1	2.1#\\
  *
  * It is permisible for the dimension #m# and #n# in case 1 to be 0.
  * In this case there will be no elements.  So to input an empty matrix you would use:
  *
  * #0  0#\\
  *
  * If one of the dimenstions is zero but not the other then a #std::length_error# 
  * exception will be thrown.
  * If any of the input operations fails then a LinAlgPackIO::InputException exception
  * is thrown.  In other words if #is.fail()# or #is.eof()# is true
  * before all of the elements have been read in then the exception is thrown.
  * Also if the stream becomes corrupted (#is.bad() == true#) then a #std::ios_base::failure#
  * exception is thrown. 
  */
// @{

/** \brief . */
/* * DMatrix input stream function.
  *
  * Inputs a DMatrix object from an input stream.
  * If #exta_flags & LinAlgPackIO::ignore_dim_bit != 0# then #gm# is resized to #m# x #n#
  * given in the file.  If #exta_flags & LinAlgPackIO::ignore_dim_bit == 0#
  * then the number of elements read in depends on the current dimension of #gm#.
  */
std::istream& input(std::istream& is, DMatrix* gm, LinAlgPackIO::fmtflags extra_flags);

/** \brief . */
/* * DMatrixSlice input stream function.
  *
  * Inputs a DMatrixSlice object from an input stream.
  * If #exta_flags & LinAlgPackIO::ignore_dim_bit != 0# then the dimension (sized) of #gms#
  * is compared to the #m# and #n# given in the file and if they are not equal
  * then a #LinAlgPackIO::InputException# is thrown.  If #gms# is unsized then it is resized
  * to #m# x #n#.  If #exta_flags & LinAlgPackIO::ignore_dim_bit == 0# then the number of
  * elements read in depends on the current size of #gms#.
  */
std::istream& input(std::istream& is, DMatrixSlice* gms, LinAlgPackIO::fmtflags extra_flags);

// @}

}	// end namespace DenseLinAlgPack

#endif	// GENMATRIX_IN_FUNC_H
