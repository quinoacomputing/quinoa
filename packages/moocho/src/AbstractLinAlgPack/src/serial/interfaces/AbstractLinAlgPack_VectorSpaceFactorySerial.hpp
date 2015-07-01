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

#ifndef VECTOR_SPACE_FACTORY_SERIAL_H
#define VECTOR_SPACE_FACTORY_SERIAL_H

#include "AbstractLinAlgPack_Types.hpp"
#include "AbstractLinAlgPack_VectorSpaceFactory.hpp"

namespace AbstractLinAlgPack {

/** \brief Abstract interface for objects that can create vector spaces of a specified dimension.
 *
 * ToDo: Finish documentation!
 */
class VectorSpaceFactorySerial
  : public AbstractLinAlgPack::VectorSpaceFactory // doxygen needs full name
{
public:

  /** \brief . */
  VectorSpaceFactorySerial( const inner_prod_ptr_t& inner_prod = Teuchos::null );

  /** @name Overridden from VectorSpaceFactory */
  //@{

  /** \brief . */
  space_ptr_t create_vec_spc(index_type dim) const;

  //@}
  
}; // end class VectorSpaceFactorySerial

} // end namespace AbstractLinAlgPack

#endif  // VECTOR_SPACE_FACTORY_SERIAL_H
