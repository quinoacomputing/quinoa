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

#include <stdexcept>
#include <string>
#include <typeinfo>

#include "AbstractLinAlgPack_AssertOp.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "Teuchos_Assert.hpp"

// boilerplate code

namespace {

struct dump_vec_spaces {
public:
  dump_vec_spaces(
    const AbstractLinAlgPack::VectorSpace& _vec_space1, const char _vec_space1_name[]
    ,const AbstractLinAlgPack::VectorSpace& _vec_space2, const char _vec_space2_name[]
    )
    :vec_space1(_vec_space1),vec_space1_name(_vec_space1_name)
    ,vec_space2(_vec_space2),vec_space2_name(_vec_space2_name)
    {}
  const AbstractLinAlgPack::VectorSpace &vec_space1;
  const char                            *vec_space1_name;
  const AbstractLinAlgPack::VectorSpace &vec_space2;
  const char                            *vec_space2_name;
}; // end dum_vec_spaces

// Notice!!!!!!!  Place a breakpoint in following function in order to halt the
// program just before an exception is thrown!

std::ostream& operator<<( std::ostream& o, const dump_vec_spaces& d )
{
  o << "Error, " << d.vec_space1_name << " at address " << &d.vec_space1
    << " of type \'" << typeName(d.vec_space1)
    << "\' with dimension " << d.vec_space1_name << ".dim() = " << d.vec_space1.dim()
    << " is not compatible with "
    << d.vec_space2_name  << " at address " << &d.vec_space2
    << " of type \'" << typeName(d.vec_space2)
    << "\' with dimension " << d.vec_space2_name << ".dim() = " << d.vec_space2.dim();
  return o;
}

enum EM_VS { SPACE_COLS, SPACE_ROWS };

const AbstractLinAlgPack::VectorSpace& op(
  const AbstractLinAlgPack::MatrixOp&     M
  ,BLAS_Cpp::Transp                           M_trans
  ,EM_VS                                      M_VS
  )
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  if(M_trans == no_trans && M_VS == SPACE_COLS)
    return M.space_cols();
  if(M_trans == trans && M_VS == SPACE_COLS)
    return M.space_rows();
  if(M_trans == no_trans && M_VS == SPACE_ROWS)
    return M.space_rows();
  // M_trans == trans && M_VS == SPACE_ROWS
  return M.space_cols();
}

} // end namespace

#define ASSERT_LHS_ARG(FUNC_NAME,LHS_ARG)                                              \
  TEUCHOS_TEST_FOR_EXCEPTION(                                                                   \
    (LHS_ARG) == NULL, std::invalid_argument                                       \
    ,FUNC_NAME << " : Error!"                                                      \
    );

// Notice!!!!!!!  Setting a breakpoint a the line number that is printed by this macro
// and then trying to set the condition !is_compatible does not work (at least not
// in gdb).

#define ASSERT_VEC_SPACES_NAMES(FUNC_NAME,VS1,VS1_NAME,VS2,VS2_NAME)                   \
{                                                                                      \
  const bool is_compatible = (VS1).is_compatible(VS2);                               \
  TEUCHOS_TEST_FOR_EXCEPTION(                                                                   \
    !is_compatible, VectorSpace::IncompatibleVectorSpaces                          \
    ,FUNC_NAME << " : "	<< dump_vec_spaces(VS1,VS1_NAME,VS2,VS2_NAME)              \
    )                                                                              \
}

#define ASSERT_VEC_SPACES(FUNC_NAME,VS1,VS2)                                              \
ASSERT_VEC_SPACES_NAMES(FUNC_NAME,VS1,#VS1,VS2,#VS2)

#define ASSERT_MAT_VEC_SPACES(FUNC_NAME,M,M_T,M_VS,VS)                                    \
{                                                                                         \
  std::ostringstream M_VS_name;                                                         \
  M_VS_name << "(" #M << ( M_T == BLAS_Cpp::no_trans ? "" : "'" ) << ")"                \
         << "." << ( M_VS == SPACE_COLS ? "space_cols()" : "space_rows()" );        \
  ASSERT_VEC_SPACES_NAMES(                                                              \
    FUNC_NAME                                                                         \
    ,op(M,M_T,M_VS),M_VS_name.str().c_str()                                           \
    ,VS,#VS                                                                           \
    )                                                                                 \
}

#define ASSERT_MAT_MAT_SPACES(FUNC_NAME,M1,M1_T,M1_VS,M2,M2_T,M2_VS)                      \
{                                                                                         \
  std::ostringstream M1_VS_name, M2_VS_name;                                            \
  M1_VS_name << "(" #M1 << ( M1_T == BLAS_Cpp::no_trans ? "" : "'" ) << ")" \
         << "." << ( M1_VS == SPACE_COLS ? "space_cols()" : "space_rows()" );       \
  M2_VS_name << "(" #M2 << ( M2_T == BLAS_Cpp::no_trans ? "" : "'" ) << ")" \
         << "." << ( M2_VS == SPACE_COLS ? "space_cols()" : "space_rows()" );       \
  ASSERT_VEC_SPACES_NAMES(                                                              \
    FUNC_NAME                                                                         \
    ,op(M1,M1_T,M1_VS),M1_VS_name.str().c_str()                                       \
    ,op(M2,M2_T,M2_VS),M2_VS_name.str().c_str()                                       \
    )                                                                                 \
}

// function definitions

#ifdef ABSTRACTLINALGPACK_ASSERT_COMPATIBILITY

void AbstractLinAlgPack::Vp_V_assert_compatibility(VectorMutable* v_lhs, const Vector& v_rhs)
{
  const char func_name[] = "Vp_V_assert_compatibility(v_lhs,v_rhs)";
  ASSERT_LHS_ARG(func_name,v_lhs)
  ASSERT_VEC_SPACES("Vp_V_assert_compatibility(v_lhs,v_rhs)",v_lhs->space(),v_rhs.space());
}

void AbstractLinAlgPack::Vp_V_assert_compatibility(VectorMutable* v_lhs, const SpVectorSlice& sv_rhs)
{
  // ToDo: Check compatibility!
} 

void AbstractLinAlgPack::VopV_assert_compatibility(const Vector& v_rhs1, const Vector&  v_rhs2)
{
  const char func_name[] = "VopV_assert_compatibility(v_rhs1,v_rhs2)";
  ASSERT_VEC_SPACES(func_name,v_rhs1.space(),v_rhs2.space());
}


void AbstractLinAlgPack::VopV_assert_compatibility(const Vector& v_rhs1, const SpVectorSlice& sv_rhs2)
{
  // ToDo: Check compatibility!
} 

void AbstractLinAlgPack::VopV_assert_compatibility(const SpVectorSlice& sv_rhs1, const Vector& v_rhs2)
{
  // ToDo: Check compatibility!
} 

void AbstractLinAlgPack::Mp_M_assert_compatibility(
  MatrixOp* m_lhs, BLAS_Cpp::Transp trans_lhs
  ,const MatrixOp& m_rhs, BLAS_Cpp::Transp trans_rhs )
{
  const char func_name[] = "Mp_M_assert_compatibility(m_lhs,trans_lhs,m_rhs,trans_rhs)";
  ASSERT_LHS_ARG(func_name,m_lhs)
  ASSERT_MAT_MAT_SPACES(func_name,(*m_lhs),trans_lhs,SPACE_COLS,m_rhs,trans_rhs,SPACE_COLS)
  ASSERT_MAT_MAT_SPACES(func_name,(*m_lhs),trans_lhs,SPACE_ROWS,m_rhs,trans_rhs,SPACE_ROWS)
}

void AbstractLinAlgPack::MopM_assert_compatibility(
  const MatrixOp& m_rhs1, BLAS_Cpp::Transp trans_rhs1
  ,const MatrixOp& m_rhs2, BLAS_Cpp::Transp trans_rhs2 )
{
  const char func_name[] = "MopM_assert_compatibility(m_rhs1,trans_rhs1,m_rhs2,trans_rhs2)";
  ASSERT_MAT_MAT_SPACES(func_name,m_rhs1,trans_rhs1,SPACE_COLS,m_rhs2,trans_rhs2,SPACE_COLS)
  ASSERT_MAT_MAT_SPACES(func_name,m_rhs1,trans_rhs1,SPACE_ROWS,m_rhs2,trans_rhs2,SPACE_ROWS)
}

void AbstractLinAlgPack::MtV_assert_compatibility(
  const MatrixOp& m_rhs1, BLAS_Cpp::Transp trans_rhs1, const Vector& v_rhs2 )
{
  const char func_name[] = "MtV_assert_compatibility(m_rhs1,trans_rhs1,v_rhs2)";
  ASSERT_MAT_VEC_SPACES(func_name,m_rhs1,trans_rhs1,SPACE_ROWS,v_rhs2.space())
}

void AbstractLinAlgPack::MtV_assert_compatibility(
  const MatrixOp& m_rhs1, BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2 )
{
  // ToDo: Check compatibility!
} 

void AbstractLinAlgPack::Vp_MtV_assert_compatibility(
  VectorMutable* v_lhs
  ,const MatrixOp& m_rhs1, BLAS_Cpp::Transp trans_rhs1, const Vector& v_rhs2 )
{
  const char func_name[] = "Vp_MtV_assert_compatibility(v_lhs,m_rhs1,trans_rhs1,v_rhs2)";
  ASSERT_LHS_ARG(func_name,v_lhs)
  ASSERT_MAT_VEC_SPACES(func_name,m_rhs1,trans_rhs1,SPACE_COLS,v_lhs->space())
  ASSERT_MAT_VEC_SPACES(func_name,m_rhs1,trans_rhs1,SPACE_ROWS,v_rhs2.space())
}

void AbstractLinAlgPack::Vp_MtV_assert_compatibility(
  VectorMutable* v_lhs
  ,const MatrixOp& m_rhs1, BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2 )
{
  // ToDo: Check compatibility!
} 

void AbstractLinAlgPack::MtM_assert_compatibility(
  const MatrixOp& m_rhs1, BLAS_Cpp::Transp trans_rhs1
  ,const MatrixOp& m_rhs2, BLAS_Cpp::Transp trans_rhs2 )
{
  const char func_name[] = "MtM_assert_compatibility(m_rhs1,trans_rhs1,m_rhs2,trans_rhs2)";
  ASSERT_MAT_MAT_SPACES(func_name,m_rhs1,trans_rhs1,SPACE_COLS,m_rhs2,trans_rhs2,SPACE_ROWS)
  ASSERT_MAT_MAT_SPACES(func_name,m_rhs1,trans_rhs1,SPACE_ROWS,m_rhs2,trans_rhs2,SPACE_COLS)
}

void AbstractLinAlgPack::Mp_MtM_assert_compatibility(
  MatrixOp* m_lhs, BLAS_Cpp::Transp trans_lhs
  ,const MatrixOp& m_rhs1, BLAS_Cpp::Transp trans_rhs1
  ,const MatrixOp& m_rhs2, BLAS_Cpp::Transp trans_rhs2 )
{
  const char func_name[] = "Mp_MtM_assert_compatibility(m_lhs,trans_lhsm_rhs1,trans_rhs1,m_rhs2,trans_rhs2)";
  ASSERT_LHS_ARG(func_name,m_lhs)
  ASSERT_MAT_MAT_SPACES(func_name,(*m_lhs),trans_lhs,SPACE_COLS,m_rhs1,trans_rhs1,SPACE_COLS)
  ASSERT_MAT_MAT_SPACES(func_name,(*m_lhs),trans_lhs,SPACE_ROWS,m_rhs2,trans_rhs2,SPACE_ROWS)
  ASSERT_MAT_MAT_SPACES(func_name,m_rhs1,trans_rhs1,SPACE_ROWS,m_rhs2,trans_rhs2,SPACE_COLS)
}

#endif // ABSTRACTLINALGPACK_ASSERT_COMPATIBILITY
