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

#include "ConstrainedOptPack_DecompositionSystemVarReductPermStd.hpp"
#include "ConstrainedOptPack_DecompositionSystemVarReductImp.hpp"
#include "AbstractLinAlgPack_MatrixOpNonsing.hpp"
#include "AbstractLinAlgPack_BasisSystemPerm.hpp"
#include "AbstractLinAlgPack_PermutationOut.hpp"
#include "AbstractLinAlgPack_MatrixOpOut.hpp"
#include "Teuchos_Assert.hpp"

namespace ConstrainedOptPack {

// Constructors / initializers

DecompositionSystemVarReductPermStd::DecompositionSystemVarReductPermStd(
  const decomp_sys_imp_ptr_t&        decomp_sys_imp
  ,const basis_sys_ptr_t&            basis_sys
  ,bool                              basis_selected
  ,EExplicitImplicit                 D_imp
  ,EExplicitImplicit                 Uz_imp
  )
{
  this->initialize(decomp_sys_imp,basis_sys,basis_selected,D_imp,Uz_imp);
}

void DecompositionSystemVarReductPermStd::initialize(
  const decomp_sys_imp_ptr_t&        decomp_sys_imp
  ,const basis_sys_ptr_t&            basis_sys
  ,bool                              basis_selected
  ,EExplicitImplicit                 D_imp
  ,EExplicitImplicit                 Uz_imp
  )
{
  decomp_sys_imp_ = decomp_sys_imp;
  basis_sys_      = basis_sys;
  basis_selected_ = basis_selected;
  this->D_imp(D_imp);
  this->Uz_imp(Uz_imp);
}

// Overridden from DecompositionSystem

size_type DecompositionSystemVarReductPermStd::n() const
{
  return decomp_sys_imp()->n();
}

size_type DecompositionSystemVarReductPermStd::m() const
{
  return decomp_sys_imp()->m();
}

size_type DecompositionSystemVarReductPermStd::r() const
{
  return decomp_sys_imp()->r();
}

Range1D DecompositionSystemVarReductPermStd::equ_decomp() const
{
  return decomp_sys_imp()->equ_decomp();
}

Range1D DecompositionSystemVarReductPermStd::equ_undecomp() const
{
  return decomp_sys_imp()->equ_undecomp();
}

const VectorSpace::space_ptr_t
DecompositionSystemVarReductPermStd::space_range() const
{
  return decomp_sys_imp()->space_range();
}

const VectorSpace::space_ptr_t
DecompositionSystemVarReductPermStd::space_null() const
{
  return decomp_sys_imp()->space_null();
}

const DecompositionSystem::mat_fcty_ptr_t
DecompositionSystemVarReductPermStd::factory_Z() const
{
  return decomp_sys_imp()->factory_Z();
}

const DecompositionSystem::mat_fcty_ptr_t
DecompositionSystemVarReductPermStd::factory_Y() const
{
  return decomp_sys_imp()->factory_Y();
}

const DecompositionSystem::mat_nonsing_fcty_ptr_t
DecompositionSystemVarReductPermStd::factory_R() const
{
  mat_nonsing_fcty_ptr_t factory_R = decomp_sys_imp()->factory_R();
  if( factory_R.get() != NULL )
    return factory_R;
  // Else assume that R will just be the basis matrix (coordinate decomposition!)
  return basis_sys_->factory_C();
}

const DecompositionSystem::mat_fcty_ptr_t
DecompositionSystemVarReductPermStd::factory_Uz() const
{
  return decomp_sys_imp()->factory_Uz();
}

const DecompositionSystem::mat_fcty_ptr_t
DecompositionSystemVarReductPermStd::factory_Uy() const
{
  return decomp_sys_imp()->factory_Uy();
}

void DecompositionSystemVarReductPermStd::update_decomp(
  std::ostream          *out
  ,EOutputLevel         olevel
  ,ERunTests            test_what
  ,const MatrixOp       &Gc
  ,MatrixOp             *Z
  ,MatrixOp             *Y
  ,MatrixOpNonsing      *R
  ,MatrixOp             *Uz
  ,MatrixOp             *Uy
  ,EMatRelations        mat_rel
  ) const
{
  assert_basis_selected();
  decomp_sys_imp()->update_decomp(
    out,olevel,test_what,Gc,Z,Y
    ,R,Uz,Uy,mat_rel
    );
}

void DecompositionSystemVarReductPermStd::print_update_decomp(
  std::ostream& out, const std::string& L ) const
{
  // ToDo: Print basis permutation stuff also?
  decomp_sys_imp()->print_update_decomp(out,L);
}

// Overridden from DecompositionSystemVarReduct

Range1D DecompositionSystemVarReductPermStd::var_indep() const
{
  return basis_sys_.get() ? basis_sys_->var_indep() : Range1D::Invalid;
}

Range1D DecompositionSystemVarReductPermStd::var_dep() const
{
  return basis_sys_.get() ? basis_sys_->var_dep() : Range1D::Invalid;
}

// @name Overridden from DecompositionSystemVarReductPerm

const DecompositionSystemVarReductPerm::perm_fcty_ptr_t
DecompositionSystemVarReductPermStd::factory_P_var() const
{
  return basis_sys()->factory_P_var();
}

const DecompositionSystemVarReductPerm::perm_fcty_ptr_t
DecompositionSystemVarReductPermStd::factory_P_equ() const
{
  return basis_sys()->factory_P_equ();
}

bool DecompositionSystemVarReductPermStd::has_basis() const
{
  return basis_selected_;
}

void DecompositionSystemVarReductPermStd::set_decomp(
  std::ostream          *out
  ,EOutputLevel         olevel
  ,ERunTests            test_what
  ,const Permutation    &P_var
  ,const Range1D        &var_dep
  ,const Permutation    *P_equ
  ,const Range1D        *equ_decomp
  ,const MatrixOp       &Gc
  ,MatrixOp             *Z
  ,MatrixOp             *Y
  ,MatrixOpNonsing      *R
  ,MatrixOp             *Uz
  ,MatrixOp             *Uy
  ,EMatRelations        mat_rel
  )
{
  // Forward these setting on to the implementation.
  decomp_sys_imp_->D_imp(  this->D_imp()  );
  decomp_sys_imp_->Uz_imp( this->Uz_imp() );
  // Get smart pointers to the basis matrix and the direct sensistivity matrices
  // and remove references to these matrix objects from the other decomposition
  // matrices by uninitializing them.
  Teuchos::RCP<MatrixOpNonsing>  C_ptr;
  Teuchos::RCP<MatrixOp>         D_ptr;
  decomp_sys_imp_->get_basis_matrices(
    out, olevel, test_what
    ,Z, Y, R, Uz, Uy
    ,&C_ptr
    ,&D_ptr // May return D_ptr.get() == NULL if not explicit chosen
    );
  // Tell the basis system object to set this basis
  try {
    basis_sys_->set_basis(
      P_var, var_dep
      ,P_equ, equ_decomp
      ,Gc
      ,C_ptr.get()
      ,D_ptr.get() // May be NULL
      ,this->Uz_imp() == MAT_IMP_EXPLICIT ? Uz : NULL
      ,(mat_rel == MATRICES_INDEP_IMPS
        ? BasisSystem::MATRICES_INDEP_IMPS : BasisSystem::MATRICES_ALLOW_DEP_IMPS )
      ,out
      );
  }
  catch( const BasisSystem::SingularBasis& except ) {
    if(out && olevel >= PRINT_BASIC_INFO)
      *out << "Passed in basis is singular, throwing SingularDecomposition: "
         << except.what() << std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, SingularDecomposition
      ,"DecompositionSystemVarReductPermStd::set_decomp(...): Passed in basis selection "
      "gave a singular basis matrix! : " << except.what() );
  }
  // If we get here the passed in basis selection is nonsingular and the basis matrices
  // are updated.  Now give them back to the decomp_sys_imp object and update the rest
  // of the decomposition matrices.
  const size_type
    m  = Gc.cols(),
    r  = C_ptr->rows();
  decomp_sys_imp_->set_basis_matrices(
    out, olevel, test_what
    ,C_ptr
    ,D_ptr // D_ptr.get() may be NULL
    ,r > m ? Uz : NULL
    ,basis_sys_ // Always reset
    );
  C_ptr = Teuchos::null;
  D_ptr = Teuchos::null;
  decomp_sys_imp()->update_decomp(
    out,olevel,test_what,Gc,Z,Y,R
    ,r > m ? Uz : NULL
    ,r > m ? Uy : NULL
    ,mat_rel
    );
  // We have a basis!
  basis_selected_ = true;
}

void DecompositionSystemVarReductPermStd::select_decomp(
  std::ostream          *out
  ,EOutputLevel         olevel
  ,ERunTests            test_what
  ,const Vector         *nu
  ,MatrixOp             *Gc
  ,Permutation          *P_var
  ,Range1D              *var_dep
  ,Permutation          *P_equ
  ,Range1D              *equ_decomp
  ,MatrixOp             *Z
  ,MatrixOp             *Y
  ,MatrixOpNonsing      *R
  ,MatrixOp             *Uz
  ,MatrixOp             *Uy
  ,EMatRelations        mat_rel
  )
{
  // Forward these setting on to the implementation.
  decomp_sys_imp_->D_imp(  this->D_imp()  );
  decomp_sys_imp_->Uz_imp( this->Uz_imp() );
  // Get smart pointers to the basis matrix and the direct sensistivity matrices
  // and remove references to these matrix objects from the other decomposition
  // matrices by uninitializing them.
  Teuchos::RCP<MatrixOpNonsing>  C_ptr;
  Teuchos::RCP<MatrixOp>             D_ptr;
  //const bool unintialized_basis = decomp_sys_imp_->basis_sys()->var_dep().size() == 0;
  decomp_sys_imp_->get_basis_matrices(
    out, olevel, test_what
    ,Z, Y, R, Uz, Uy
    ,&C_ptr
    ,&D_ptr // May return D_ptr.get() == NULL if not explicit chosen
    );
  // Ask the basis system object to select a basis
  basis_sys_->select_basis(
    nu
    ,Gc
    ,P_var, var_dep
    ,P_equ, equ_decomp
    ,C_ptr.get()
    ,D_ptr.get() // May be NULL
    ,this->Uz_imp() == MAT_IMP_EXPLICIT ? Uz : NULL
    ,(mat_rel == MATRICES_INDEP_IMPS
      ? BasisSystem::MATRICES_INDEP_IMPS : BasisSystem::MATRICES_ALLOW_DEP_IMPS )
    ,out
    );

  if( out && (int)olevel >= (int)PRINT_BASIC_INFO ) {
    const Range1D var_indep = basis_sys_->var_indep(), equ_undecomp = basis_sys_->equ_undecomp();
    *out
      << "\nSelected a new basis\n"
      << "\nbs.var_dep()            = ["<<var_dep->lbound()<<","<<var_dep->ubound()<<"]"
      << "\nds.var_indep()          = ["<<var_indep.lbound()<<","<<var_indep.ubound()<<"]"
      << "\nds.equ_decomp()         = ["<<equ_decomp->lbound()<<","<<equ_decomp->ubound()<<"]"
      << "\nds.equ_undecomp()       = ["<<equ_undecomp.lbound()<<","<<equ_undecomp.ubound()<<"]"
      << std::endl;
  }
  if( out && (int)olevel >= (int)PRINT_VECTORS ) {
    *out
      << "\nP_var =\n" << *P_var
      << "\nP_equ =\n" << *P_equ
      ;
  }
  if( out && (int)olevel >= (int)PRINT_EVERY_THING ) {
    *out
      << "\nGc =\n" << *Gc;
  }

  // If we get here a nonsinguar basis selection has been made and the basis matrices
  // are updated.  Now give them back to the decomp_sys_imp object and update the rest
  // of the decomposition matrices.
  const size_type
    //n  = Gc->rows(),
    m  = Gc->cols(),
    r  = C_ptr->rows();
  decomp_sys_imp_->set_basis_matrices(
    out, olevel, test_what
    ,C_ptr
    ,D_ptr // D_ptr.get() may be NULL
    ,r > m ? Uz : NULL
    ,basis_sys_ // Always reset
    );
  C_ptr = Teuchos::null;
  D_ptr = Teuchos::null;
  decomp_sys_imp()->update_decomp(
    out,olevel,test_what,*Gc,Z,Y,R
    ,r > m ? Uz : NULL
    ,r > m ? Uy : NULL
    ,mat_rel
    );
  // We have a basis!
  basis_selected_ = true;
}

// private

void DecompositionSystemVarReductPermStd::assert_basis_selected() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !basis_selected_, std::logic_error
    ,"DecompositionSystemVarReductPermStd::assert_basis_selected(): Error, "
    "the methods set_decomp() or select_decomp() must be called first!" );
}

}	// end namespace ConstrainedOptPack
