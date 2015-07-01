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

#ifndef LINALGPACK_OUT_FORMAT_H
#define LINALGPACK_OUT_FORMAT_H

#include "DenseLinAlgPack_IOFormat.hpp"

namespace DenseLinAlgPack {

/** \brief . */
/* * Output stream operator function for const_bound_format objects.
  *
  * This template function performs the following tasks.
  * \begin{enumeration}
  * <li> Saves the formating state of #os#
  * <li> Sets the formating state of #os# to that stored in #bf.f()#
  * <li> Calls #output(os, bf.obj(), bf.extra_flags().flags())#
  * <li> Resets the streams formating state to its original.
  * \end{enumeration}
  *
  * The original formating state of #os# is preserved even if an exception
  * is thrown. 
  */
template<class T>
std::ostream& operator<<(std::ostream& os, const LinAlgPackIO::const_bound_format<T>& bf)
{
  using LinAlgPackIO::ios_format_memento;
  ios_format_memento old_format = ios_format_memento::save_format(os);
  try {
    bf.f().set_format(os);
    output( os, bf.obj(), bf.f().extra_flags().flags() );
  }
  catch(...) {
    old_format.set_format(os);
    throw;
  }
  old_format.set_format(os);
  return os;
}

/// Force a type conversion from #bound_format<T># to #const_bound_format<T># to call #operator<<#().
template<class T>
inline std::ostream& operator<<(std::ostream& os, const LinAlgPackIO::bound_format<T>& bf) {
  return operator<<( os, LinAlgPackIO::const_bound_format<T>( bf.f(), bf.obj() ) );
}

}	// end namespace DenseLinAlgPack

#endif // LINALGPACK_OUT_FORMAT_H
