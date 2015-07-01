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

#ifndef ABSTRACT_LIN_ALG_PACK_BASIS_SYSTEM_FACTORY_H
#define ABSTRACT_LIN_ALG_PACK_BASIS_SYSTEM_FACTORY_H

#include "AbstractLinAlgPack_Types.hpp"
#include "Teuchos_AbstractFactory.hpp"
#include "Teuchos_RCP.hpp"

namespace OptionsFromStreamPack {
  class OptionsFromStream;
}

namespace AbstractLinAlgPack {

/** \brief Interface for a factory object that will create <tt>BasisSystem</tt> objects.
 *
 * 
 */
class BasisSystemFactory : public Teuchos::AbstractFactory<BasisSystem>
{
public:

  /** @name Public types */
  //@{

  /** \brief . */
  typedef Teuchos::RCP<
    const OptionsFromStreamPack::OptionsFromStream>             options_ptr_t;

  //@}

  /** \brief . */
  virtual ~BasisSystemFactory() {}

  /** \brief Set the options that will be used to determine what basis system will be returned
   * from <tt>this->create()</tt>.
   *
   * Note that it is allowed for the client to alter <tt>*options.get()</tt> after
   * this method is called so <tt>this</tt> had better read the options inside of the
   * <tt>this->create()</tt> method.
   */
  virtual void set_options( const options_ptr_t& options ) = 0;

  /** \brief Get the <tt>OptionsFromStream</tt> object being used to extract the options from.
   */
  virtual const options_ptr_t& get_options() const = 0;

}; // end class BasisSystemFactory

}  // end namespace AbstractLinAlgPack

#endif // ABSTRACT_LIN_ALG_PACK_BASIS_SYSTEM_FACTORY_H
