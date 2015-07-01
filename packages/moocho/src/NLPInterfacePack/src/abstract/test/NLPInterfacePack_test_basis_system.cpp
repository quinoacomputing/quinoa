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

#include "NLPInterfacePack_test_basis_system.hpp"
#include "NLPInterfacePack_NLPFirstOrder.hpp"
#include "AbstractLinAlgPack_BasisSystem.hpp"
#include "AbstractLinAlgPack_BasisSystemTester.hpp"
#include "AbstractLinAlgPack_BasisSystemTesterSetOptions.hpp"
#include "AbstractLinAlgPack_MatrixOpNonsing.hpp"

bool NLPInterfacePack::test_basis_system(
   NLPFirstOrder                                 *nlp
  ,BasisSystem                                  *basis_sys
  ,OptionsFromStreamPack::OptionsFromStream     *options
  ,std::ostream                                 *out
  )
{
  namespace mmp = MemMngPack;

  const index_type
    n = nlp->n(),
    m = nlp->m();

  // Create the matrices Gc and Gh
  NLPFirstOrder::mat_fcty_ptr_t::element_type::obj_ptr_t
    Gc = ( m  ? nlp->factory_Gc()->create() : Teuchos::null );
  
  // Compute the matrices at xinit
  const Vector
    &xo = nlp->xinit();
  if(m)
    nlp->set_Gc(Gc.get());
  if(m)
    nlp->calc_Gc(xo);

  // Create the matrices C and D
  BasisSystem::mat_nonsing_fcty_ptr_t::element_type::obj_ptr_t
    C = ( m ? basis_sys->factory_C()->create() : Teuchos::null);
  BasisSystem::mat_fcty_ptr_t::element_type::obj_ptr_t
    D = ( m && n > m && basis_sys->factory_C().get() ? basis_sys->factory_C()->create() : Teuchos::null);
  BasisSystem::mat_fcty_ptr_t::element_type::obj_ptr_t
    GcUP = ( m && n > m && basis_sys->factory_GcUP().get()  ? basis_sys->factory_GcUP()->create() : Teuchos::null);

  // Initialize C and D with basis_sys
  basis_sys->update_basis(
    *Gc
    ,C.get()
    ,D.get()
    ,GcUP.get()
    );

  // Test the basis and basis system objects.
  BasisSystemTester
    basis_sys_tester;
  if(options) {
    BasisSystemTesterSetOptions
      opt_setter(&basis_sys_tester);
    opt_setter.set_options(*options);
  }
  const bool result = basis_sys_tester.test_basis_system(
    *basis_sys
    ,Gc.get()
    ,C.get()
    ,NULL    // Create the N matrix internally
    ,D.get()
    ,GcUP.get()
    ,out
    );

  return result;
}
