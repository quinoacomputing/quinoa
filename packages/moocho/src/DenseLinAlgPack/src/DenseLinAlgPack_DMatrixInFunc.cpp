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

#include <sstream>

#include "InputStreamHelperPack_EatInputComment.hpp"
#include "DenseLinAlgPack_DMatrixInFunc.hpp"
#include "DenseLinAlgPack_DVectorInFunc.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"

namespace {	// Local inplementation
std::istream& input_gms(std::istream& is, DenseLinAlgPack::DMatrixSlice* gms, const char func[]);
}

std::istream& DenseLinAlgPack::input(std::istream& is, DMatrix* gm, LinAlgPackIO::fmtflags extra_flags) {
  if( !(extra_flags & LinAlgPackIO::ignore_dim_bit) ) {
    size_type m, n;
    is >> m >> n;
    if(is.fail())
      throw LinAlgPackIO::InputException( "DenseLinAlgPack::input() {DMatrix}: "
        "Input operation of matrix dimension failed.  Check that the constant n "
        "is a valid integer." );
    if(is.bad())
      throw std::ios_base::failure( "DenseLinAlgPack::input() {DMatrix}: "
        "Input operation failed because the stream became currupted." );
    gm->resize(m,n);
  }
  DMatrixSlice gms = (*gm)();
  return input_gms(is,&gms,"DenseLinAlgPack::input() {DMatrix}");
}

std::istream& DenseLinAlgPack::input(std::istream& is, DMatrixSlice* gms, LinAlgPackIO::fmtflags extra_flags) {
  if( !(extra_flags & LinAlgPackIO::ignore_dim_bit) ) {
    size_type m, n;
    is >> m >> n;
    if(is.fail())
      throw LinAlgPackIO::InputException( "DenseLinAlgPack::input() {DMatrixSlice}: "
        "Input operation of matrix dimension failed.  Check that the constant n "
        " is a valid integer.");
    if(is.bad())
      throw std::ios_base::failure( "DenseLinAlgPack::input() {DMatrixSlice}: "
        "Input operation failed because the stream became currupted." );
    DenseLinAlgPack::assert_gms_lhs(*gms,m,n);
  }
  return input_gms( is, gms, "DenseLinAlgPack::input() {DMatrixSlice}" );
}

// //////////////////////
// Local implementation

namespace {

// Read in a specified number of elements into a DMatrixSlice object.
// The dim of gms is not checked.  If an element input operation fails or the end of the file
// is reached before all of the elements are read in then a LinAlgPackIO::InputException is thrown.
// If the stream becomes currupted durring the input then a std::ios_base::failure exception
// is thrown.  The state of the input steam remains the same on return accept for the char's
// that have been extracted.
std::istream& input_gms(std::istream& is, DenseLinAlgPack::DMatrixSlice* gms, const char func[]) {
  using std::ios_base;
  using DenseLinAlgPack::size_type;
  using DenseLinAlgPack::DVectorSlice;
  if(!gms->rows()) return is;	// If we are inputting an unsized matrix then there are no elements
                // to extract so just return.
  ios_base::iostate old_state = is.exceptions();		// save the old state
  is.exceptions(ios_base::badbit | ios_base::failbit | ios_base::eofbit);
  try {
    // Read in the rows
    for(size_type i = 1; i <= gms->rows(); ++i) {
      InputStreamHelperPack::eat_comment_lines(is,'*');
      DVectorSlice gms_row_i = gms->row(i);
      DenseLinAlgPack::input( is, &gms_row_i
        , (DenseLinAlgPack::LinAlgPackIO::fmtflags)(DenseLinAlgPack::LinAlgPackIO::ignore_dim_bit) );
    }
  }
  catch(...) {
    is.exceptions(old_state);
    throw;
  }
  is.exceptions(old_state);
  return is;
}

}	// end namespace
