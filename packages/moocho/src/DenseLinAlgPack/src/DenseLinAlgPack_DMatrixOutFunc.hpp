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

#ifndef GENMATRIX_OUT_FUNC_H
#define GENMATRIX_OUT_FUNC_H

#include "DenseLinAlgPack_IOBasic.hpp"

namespace DenseLinAlgPack {

/* * @name DMatrixSlice output stream function.
  *
  * This is a function that are used to output a DMatrixSlice object
  * to a char based output stream.  This format is ment to be both
  * human and machine readable.  In fact the \Ref{input} function
  * can be used to read in the output produced by this function.
  * This function is ment to be the machinary for an optional output stream
  * operator \Ref{operator<<}.
  *
  * The output format is diferent depending on the on whether the
  * bits #LinAlgPackIO::ignore_dim_bit# and #LinAlgPackIO::no_insert_newlines_bit# are set.
  * The default output format is (exta_flags == 0):
  *
  * #gms.rows()  gms.cols()#\\
  * #gms(1,1)           gms(1,2)          gms(1,3)          ...      gms(1,gms.cols())#\\
  * #gms(2,1)           gms(2,2)          gms(2,3)          ...      gms(2,gms.cols())#\\
  * #   .                   .                .                               .#\\
  * #gms(gms.rows(),1) gms(gms.rows(),2)  gms(gms.rows(),3)          gms(gms.rows(),gms.cols())#\\
  *
  * If #extra_flags & LinAlgPackIO::ignore_dim_bit != 0# then the dimensions of the
  * matrix are not output and if #extra_flags & LinAlgPackIO::no_insert_newlines_bit != 0#
  * then newline char ('\n') are not inserted after each row.
  *
  * Each of the elements are are put into columns according to the width set in the 
  * output stream #os# when it is called and other formating commands.  Even if the set
  * width is 0 or less than
  * the number of char's for the element a space ' ' will be inserted between them.
  * The elements are formated according to the format in the stream #os#.
  *
  * If #gms.rows() == 0# then no elements will be output and if
  * #extra_flags & LinAlgPackIO::ignore_dim_bit != 0# and only the operation
  * #os << gms.rows() << ' ' << gms.cols() << endl;# will be performed.
  *
  * If any of the output operations fails then a #std::ios_base::failure# exception is thrown. 
  */
std::ostream& output(std::ostream& os, const DMatrixSlice& gms, LinAlgPackIO::fmtflags extra_flags);

}	// end namespace DenseLinAlgPack

#endif	// GENMATRIX_OUT_FUNC_H
