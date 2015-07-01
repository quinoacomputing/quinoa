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

#include "RTOpPack_print_sub_vector.hpp"
#include "Teuchos_FancyOStream.hpp"

std::ostream& RTOpPack::output(
  std::ostream& o_arg, const SubVector& v
  ,bool print_dim , bool newline
  )
{
  Teuchos::RCP<Teuchos::FancyOStream> o = Teuchos::getFancyOStream(Teuchos::rcp(&o_arg,false));
  //Teuchos::OSTab tab(o);
  int w = o->width(0) - 1; // get the set width (minus 1 since a space is inserted)
  if( print_dim )
    *o << std::setw(0) << std::left << v.subDim() << std::endl << std::right;
  //*o << std::setiosflags(std::ios::left) << std::setw(0) << v.subDim() 
  //	  << std::endl << std::setiosflags(std::ios::right);
  // RAB: 20030916: ToDo: Fix the above by hacking std::left and std::right in config header!
  const RTOp_value_type  *v_val        = v.values();
  const ptrdiff_t        v_val_s       = v.stride();
  for( RTOp_index_type i = 1; i <= v.subDim(); ++i, v_val+=v_val_s ) {
    // insert a space to be sure there is white space
    // inbetween adjacent elements.
    *o << " " << std::setw(w) << (*v_val) << ":" << i + v.globalOffset();
  }
  if(newline) *o << std::endl;
  return o_arg;
}
