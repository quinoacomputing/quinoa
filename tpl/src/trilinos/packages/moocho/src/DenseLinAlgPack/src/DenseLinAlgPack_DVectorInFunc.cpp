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

#include "DenseLinAlgPack_DVectorInFunc.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_AssertOp.hpp"

namespace {	// Local implementation
std::istream& input_vs(std::istream& is, DenseLinAlgPack::DVectorSlice* vs, const char func[]);
}

std::istream& DenseLinAlgPack::input(std::istream& is, DVector* v, LinAlgPackIO::fmtflags extra_flags) {
  if( !(extra_flags & LinAlgPackIO::ignore_dim_bit) ) {
    size_type n;
    is >> n;
    if(is.fail())
      throw LinAlgPackIO::InputException("DenseLinAlgPack::input() {DVector}:  Input operation of vector dimension failed.  Check that the constant n is a valid integer.");
    if(is.bad())
      throw std::ios_base::failure("DenseLinAlgPack::input() {DVector}: Input operation failed because the stream became currupted.");
    v->resize(n);
  }
  return input_vs(is,&(*v)(),"DenseLinAlgPack::input() {DVector}");
}

std::istream& DenseLinAlgPack::input(std::istream& is, DVectorSlice* vs, LinAlgPackIO::fmtflags extra_flags) {
  if( !(extra_flags & LinAlgPackIO::ignore_dim_bit) ) {
    size_type n;
    is >> n;
    if(is.fail())
      throw LinAlgPackIO::InputException("DenseLinAlgPack::input() {DVectorSlice}:  Input operation of vector dimension failed.  Check that the constant n is a valid integer.");
    if(is.bad())
      throw std::ios_base::failure("DenseLinAlgPack::input() {DVectorSlice}: Input operation failed because the stream became currupted.");
    DenseLinAlgPack::Vp_V_assert_sizes( vs->dim(), n );
  }
  return input_vs(is,vs,"DenseLinAlgPack::input() {DVectorSlice}");
}


// ///////////////////////////////////
// Local implementation

namespace {

// Read in a specified number of elements into a DVectorSlice object.
// The dim of vs is not checked.  If an element input operation fails or the end of the file
// is reached before all of the elements are read in then a LinAlgPackIO::InputException is thrown.
// If the stream becomes currupted durring the input then a std::ios_base::failure exception
// is thrown.  The state of the input steam remains the same on return accept for the char's
// that have been extracted.
std::istream& input_vs(std::istream& is, DenseLinAlgPack::DVectorSlice* vs, const char func[]) {
  using std::ios_base;
  using DenseLinAlgPack::DVectorSlice;
  if(!vs->dim()) return is;	// If there are no elements to read in just return
  ios_base::iostate old_state = is.exceptions();		// save the old state
  is.exceptions(ios_base::badbit | ios_base::failbit);
  try {
    // Read in the elements
    for(DVectorSlice::iterator itr = vs->begin(); itr != vs->end(); ++itr)
      is >> *itr;
  }
  catch(std::ios_base::failure& excpt) {
    is.exceptions(old_state);
    if(is.bad()) throw;	// The stream was bad so rethrow the exception
    if(is.fail()) {
      std::ostringstream os;
      os << func << ":  An vector element input failed.  Check that the vector element is a valid C number.  "
         << excpt.what();
      throw DenseLinAlgPack::LinAlgPackIO::InputException(os.str());			
    }
    if(is.eof()) {
      std::ostringstream os;
      os << func << ":  DVector input failed.  The end of the file was found before all of the elements where read in.  "
         << excpt.what();;
      throw DenseLinAlgPack::LinAlgPackIO::InputException(os.str());			
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
