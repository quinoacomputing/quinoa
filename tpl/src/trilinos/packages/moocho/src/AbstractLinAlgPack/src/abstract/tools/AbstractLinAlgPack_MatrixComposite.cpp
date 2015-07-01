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

#include "AbstractLinAlgPack_MatrixComposite.hpp"
#include "AbstractLinAlgPack_SpVectorClass.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_VectorMutableBlocked.hpp"
#include "AbstractLinAlgPack_AssertOp.hpp"
//#include "AbstractLinAlgPack/src/AbstractLinAlgPack_GenPermMatrixSliceOp.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_Assert.hpp"
#include "ProfileHackPack_profile_hack.hpp"

namespace {

// Get an element from a vector (two versions)

inline
AbstractLinAlgPack::value_type
get_element( const AbstractLinAlgPack::Vector& v, AbstractLinAlgPack::index_type i )
{
   return v.get_ele(i);
}

inline
AbstractLinAlgPack::value_type
get_element( const AbstractLinAlgPack::SpVectorSlice& v, AbstractLinAlgPack::index_type i )
{
  const AbstractLinAlgPack::SpVectorSlice::element_type
    *ele = v.lookup_element(i);
  return ele != NULL ? ele->value() : 0.0;
}

// Get a view of a vector (two versions)

inline
Teuchos::RCP<const AbstractLinAlgPack::Vector>
get_view(
  const AbstractLinAlgPack::Vector& v
  ,AbstractLinAlgPack::index_type l
  ,AbstractLinAlgPack::index_type u
  )
{
   return v.sub_view(l,u);
}

inline
Teuchos::RCP<const AbstractLinAlgPack::SpVectorSlice>
get_view(
  const AbstractLinAlgPack::SpVectorSlice& v
  ,AbstractLinAlgPack::index_type l
  ,AbstractLinAlgPack::index_type u
  )
{
  return Teuchos::rcp( new AbstractLinAlgPack::SpVectorSlice( v(l,u) ) );
}

// Perform a matrix vector multiplication

template<class V>
void Vp_StMtV_imp(
  AbstractLinAlgPack::VectorMutable* y, AbstractLinAlgPack::value_type a
  ,AbstractLinAlgPack::size_type M_rows, AbstractLinAlgPack::size_type M_cols
  ,const AbstractLinAlgPack::MatrixComposite::matrix_list_t& mat_list
  ,const AbstractLinAlgPack::MatrixComposite::vector_list_t& vec_list
  ,BLAS_Cpp::Transp M_trans
  ,const V& x, AbstractLinAlgPack::value_type b
  )
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  using BLAS_Cpp::rows;
  using BLAS_Cpp::cols;
  using AbstractLinAlgPack::dot;  // We should not have to do this but some compiles &%!#$
  typedef AbstractLinAlgPack::MatrixComposite::matrix_list_t  mat_list_t;
  typedef AbstractLinAlgPack::MatrixComposite::vector_list_t  vec_list_t;

    AbstractLinAlgPack::Vt_S( y, b );  // Will take care of b == 0.0

  if( vec_list.size() ) {
    for( vec_list_t::const_iterator itr = vec_list.begin(); itr != vec_list.end(); ++itr ) {
      const BLAS_Cpp::Transp
        op_v_trans = ( M_trans == itr->v_trans_ ? no_trans : trans );
      const AbstractLinAlgPack::index_type
        r = ( M_trans == no_trans ? itr->r_l_ : itr->c_l_ ),
        c = ( M_trans == no_trans ? itr->c_l_ : itr->r_l_ );
      if( itr->rng_G_.full_range() ) { // op(v)
        if( op_v_trans == no_trans ) {
          //
          //         [ y1 ]         [              ]   [ x1 ]
          // r:r+n-1 [ y2 ] +=  a * [   beta * v   ] * [ x2 ] c:c 
          //         [ y3 ]         [              ]   [ x3 ]
          //                            \______/        
          //                              c:c
          // =>
          //
          // y(r:r+n-1) += (a * beta * x(c)) * v
          // 
          AbstractLinAlgPack::Vp_StV( y->sub_view(r,r+itr->v_->dim()-1).get(), a * itr->beta_ * get_element(x,c), *itr->v_ );
        }
        else {
          //
          //     [ y1 ]        [                ]   [ x1 ]
          // r:r [ y2 ] += a * [     beta*v'    ] * [ x2 ] c:c+n-1
          //     [ y3 ]        [                ]   [ x3 ]
          //                         \_____/  
          //                         c:c+n-1
          // =>
          //
          // y(r) += a * beta * v'*x(c,c+n-1)
          //
//					y->set_ele( r, y->get_ele(r) + a * itr->beta_ * dot( *itr->v_, *get_view(x,c,c+itr->v_->dim()-1) ) );
          TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement the above method in VectorStdOps for Vector,SpVectorSlice!
        }
      }
      else { // op(op(G)*v) or op(v(rng_G))
        TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement when needed!
      }
    }
  }
  if( mat_list.size() ) {
    for( mat_list_t::const_iterator itr = mat_list.begin(); itr != mat_list.end(); ++itr ) {
      const AbstractLinAlgPack::index_type
        rl = rows(itr->r_l_,itr->c_l_,M_trans),
        ru = rows(itr->r_u_,itr->c_u_,M_trans),
        cl = cols(itr->r_l_,itr->c_l_,M_trans),
        cu = cols(itr->r_u_,itr->c_u_,M_trans);
      const BLAS_Cpp::Transp
        op_P_trans = ( M_trans == itr->P_trans_ ? no_trans : trans ),
        op_A_trans = ( M_trans == itr->A_trans_ ? no_trans : trans ),
        op_Q_trans = ( M_trans == itr->Q_trans_ ? no_trans : trans );
      if( itr->rng_P_.full_range() && itr->rng_Q_.full_range() ) { // op(A)
        //
        //       [ y1 ]        [                        ]   [ x1 ]
        // rl:ru [ y2 ] += a * [    alpha * op(op(A))   ] * [ x2 ] cl:cu
        //       [ y3 ]        [                        ]   [ x3 ]
        //                          \_______________/
        //                               cl:cu
        // =>
        //
        // y(rl:ru) += a * alpha * op(op(A)) * x(cl:cu)
        //
        AbstractLinAlgPack::Vp_StMtV( y->sub_view(rl,ru).get(), a * itr->alpha_, *itr->A_, op_A_trans, *get_view(x,cl,cu) );
      }
      else {
        if( itr->A_ == NULL ) { // op(P)
          TEUCHOS_TEST_FOR_EXCEPT( !(  itr->P_.get() && !itr->P_->is_identity()  ) );
          //
          //       [ y1 ]        [                        ]   [ x1 ]
          // rl:ru [ y2 ] += a * [    alpha * op(op(P))   ] * [ x2 ] cl:cu
          //       [ y3 ]        [                        ]   [ x3 ]
          //                          \_______________/
          //                               cl:cu
          // =>
          //
          // y(rl:ru) += a * alpha * op(op(P)) * x(cl:cu)
          //
// 					AbstractLinAlgPack::Vp_StMtV( y->sub_view(rl,ru).get(), a * itr->alpha_, itr->P_, op_P_trans, *get_view(x,cl,cu) );
          TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement the above method properly!
        }
        else { // op(P)*op(A)*op(Q)  [or some simplification]
          //
          //       [ y1 ]        [                                ]   [ x1 ]
          // rl:ru [ y2 ] += a * [    alpha * op(P)*op(A)*op(Q)   ] * [ x2 ] cl:cu
          //       [ y3 ]        [                                ]   [ x3 ]
          //                          \______________________/
          //                                    cl:cu
          // =>
          //
          // y(rl:ru) += a * alpha * op(P)*op(A)*op(Q) * x(cl:cu)
          //
          if( !itr->rng_P_.full_range() && !itr->rng_Q_.full_range() ) { // op(A)(rng_Q,rng_Q)
            AbstractLinAlgPack::Vp_StMtV(
              y->sub_view(rl,ru).get(), a * itr->alpha_
              ,*itr->A_->sub_view(
                itr->A_trans_==no_trans  ? itr->rng_P_ : itr->rng_Q_
                ,itr->A_trans_==no_trans ? itr->rng_Q_ : itr->rng_P_ )
              ,op_A_trans
              ,*get_view(x,cl,cu)
              );
          }
          else {
            TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement when needed!
          }
        }
      }
    }
  }
}

// Define a comparison object for comparing SubVectorEntry or SubMatrixEntry
// objects for storting by row or by column
template<class T>
class CompSubEntry {
public:
  enum EByRowCol { BY_ROW, BY_COL };
  CompSubEntry(enum EByRowCol by_row_col)
    : by_row_col_(by_row_col)
    {}
  bool operator()( const T& e1, const T& e2 )
    {
      return
        ( by_row_col_ == BY_ROW
          ? e1.r_l_ < e2.r_l_
          : e1.c_l_ < e2.c_l_
          );
    }
private:
  EByRowCol  by_row_col_;
  CompSubEntry(); // not defined and not to be called
}; // end class CompSubVectorEntry


} // end namespace

namespace AbstractLinAlgPack {

MatrixComposite::MatrixComposite( size_type rows, size_type cols )
{
  reinitialize(rows,cols);
}

void MatrixComposite::reinitialize( size_type rows, size_type cols )
{
  namespace rcp = MemMngPack;
  
  fully_constructed_ = true;
  rows_ = rows;
  cols_ = cols;
  if(matrix_list_.size()) {
    matrix_list_.erase(matrix_list_.begin(),matrix_list_.end());
  }
  if(vector_list_.size()) {
    vector_list_.erase(vector_list_.begin(),vector_list_.end());
  }
  space_rows_  = Teuchos::null;
  space_cols_  = Teuchos::null;
}

void MatrixComposite::add_vector(
  size_type                      row_offset
  ,size_type                     col_offset
  ,value_type                    beta
  ,const GenPermMatrixSlice      *G
  ,const release_resource_ptr_t  &G_release
  ,BLAS_Cpp::Transp              G_trans
  ,const Vector            *v
  ,const release_resource_ptr_t  &v_release
  ,BLAS_Cpp::Transp              v_trans
  )
{
  fully_constructed_ = false;
  TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Finish!
}

void MatrixComposite::add_vector(
  size_type                      row_offset
  ,size_type                     col_offset
  ,value_type                    beta
  ,const Range1D                 &rng_G
  ,const Vector            *v
  ,const release_resource_ptr_t  &v_release
  ,BLAS_Cpp::Transp              v_trans
  )
{
  fully_constructed_ = false;
  TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Finish!
}

void MatrixComposite::add_vector(
  size_type                      row_offset
  ,size_type                     col_offset
  ,value_type                    beta
  ,const Vector            *v
  ,const release_resource_ptr_t  &v_release
  ,BLAS_Cpp::Transp              v_trans
  )
{
  namespace rcp = MemMngPack;

  TEUCHOS_TEST_FOR_EXCEPT( !(  beta != 0.0  ) );
  TEUCHOS_TEST_FOR_EXCEPT( !(  v != NULL  ) );
  fully_constructed_ = false;
  if( v_trans == BLAS_Cpp::no_trans ) {
    TEUCHOS_TEST_FOR_EXCEPT( !(  row_offset + v->dim() <= rows_  ) );
    TEUCHOS_TEST_FOR_EXCEPT( !(  col_offset + 1 <= cols_  ) );
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPT( !(  row_offset + 1 <= rows_  ) );
    TEUCHOS_TEST_FOR_EXCEPT( !(  col_offset + v->dim() <= cols_  ) );
  }

  vector_list_.push_back(
    SubVectorEntry(
      row_offset+1,col_offset+1,beta,Range1D()
      ,Teuchos::null,Teuchos::null,BLAS_Cpp::no_trans
      ,v,v_release,v_trans ) );
}

void MatrixComposite::remove_vector( vector_list_t::iterator itr )
{
  fully_constructed_ = false;
  vector_list_.erase(itr);
}

void MatrixComposite::add_matrix(
  size_type                      row_offset
  ,size_type                     col_offset
  ,value_type                    alpha
  ,const GenPermMatrixSlice      *P
  ,const release_resource_ptr_t  &P_release
  ,BLAS_Cpp::Transp              P_trans
  ,const MatrixOp            *A
  ,const release_resource_ptr_t  &A_release
  ,BLAS_Cpp::Transp              A_trans
  ,const GenPermMatrixSlice      *Q
  ,const release_resource_ptr_t  &Q_release
  ,BLAS_Cpp::Transp              Q_trans
  )
{
  fully_constructed_ = false;
  TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Finish!
}

void MatrixComposite::add_matrix(
  size_type                      row_offset
  ,size_type                     col_offset
  ,value_type                    alpha
  ,const Range1D                 &rng_P
  ,const MatrixOp            *A
  ,const release_resource_ptr_t  &A_release
  ,BLAS_Cpp::Transp              A_trans
  ,const GenPermMatrixSlice      *Q
  ,const release_resource_ptr_t  &Q_release
  ,BLAS_Cpp::Transp              Q_trans
  )
{
  fully_constructed_ = false;
  TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Finish!
}

void MatrixComposite::add_matrix(
  size_type                      row_offset
  ,size_type                     col_offset
  ,value_type                    alpha
  ,const GenPermMatrixSlice      *P
  ,const release_resource_ptr_t  &P_release
  ,BLAS_Cpp::Transp              P_trans
  ,const MatrixOp            *A
  ,const release_resource_ptr_t  &A_release
  ,BLAS_Cpp::Transp              A_trans
  ,const Range1D                 &rng_Q
  )
{
  fully_constructed_ = false;
  TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Finish!
}

void MatrixComposite::add_matrix(
  size_type                      row_offset
  ,size_type                     col_offset
  ,value_type                    alpha
  ,const Range1D                 &rng_P_in
  ,const MatrixOp            *A
  ,const release_resource_ptr_t  &A_release
  ,BLAS_Cpp::Transp              A_trans
  ,const Range1D                 &rng_Q_in
  )
{
  namespace rcp = MemMngPack;
  using BLAS_Cpp::rows;
  using BLAS_Cpp::cols;
  using RangePack::full_range;

  TEUCHOS_TEST_FOR_EXCEPTION(
    alpha == 0.0, std::invalid_argument
    ,"MatrixComposite::add_matrix(...) : Error!" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    A == NULL, std::invalid_argument
    ,"MatrixComposite::add_matrix(...) : Error!" );

  const size_type
    A_rows   = A->rows(),
    A_cols   = A->cols(),
    opA_rows = rows(A_rows,A_cols,A_trans),
    opA_cols = cols(A_rows,A_cols,A_trans);
  const Range1D
    rng_P = full_range(rng_P_in,1,opA_rows),
    rng_Q = full_range(rng_Q_in,1,opA_cols);
  const size_type
    opPopAopQ_rows = rng_P.size(),
    opPopAopQ_cols = rng_Q.size();

  TEUCHOS_TEST_FOR_EXCEPTION(
    row_offset + opPopAopQ_rows > rows_, std::invalid_argument
    ,"MatrixComposite::add_matrix(...) : Error!" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    col_offset + opPopAopQ_cols > cols_, std::invalid_argument
    ,"MatrixComposite::add_matrix(...) : Error!" );

  fully_constructed_ = false;
  
  matrix_list_.push_back(
    SubMatrixEntry(
      row_offset+1,row_offset+opPopAopQ_rows
      ,col_offset+1,col_offset+opPopAopQ_cols
      ,alpha
      ,rng_P
      ,Teuchos::null,Teuchos::null,BLAS_Cpp::no_trans
      ,A,A_release,A_trans
      ,rng_Q
      ,Teuchos::null,Teuchos::null,BLAS_Cpp::no_trans
      )
    );

  // ToDo: We could create a sub-view of the matrix using rng_P and rng_Q
  // and then store the sub-view in the data structure?  However, this
  // would confuse an external client who wanted to iterate through
  // the matrix list using the public iterators.
}

void MatrixComposite::add_matrix(
  size_type                      row_offset
  ,size_type                     col_offset
  ,value_type                    alpha
  ,const MatrixOp            *A
  ,const release_resource_ptr_t  &A_release
  ,BLAS_Cpp::Transp              A_trans
  )
{
  namespace rcp = MemMngPack;
  using BLAS_Cpp::rows;
  using BLAS_Cpp::cols;

  TEUCHOS_TEST_FOR_EXCEPTION(
    alpha == 0.0, std::invalid_argument
    ,"MatrixComposite::add_matrix(...) : Error!" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    A == NULL, std::invalid_argument
    ,"MatrixComposite::add_matrix(...) : Error!" );

  const size_type
    A_rows   = A->rows(),
    A_cols   = A->cols(),
    opA_rows = rows(A_rows,A_cols,A_trans),
    opA_cols = cols(A_rows,A_cols,A_trans);

  TEUCHOS_TEST_FOR_EXCEPTION(
    row_offset + opA_rows > rows_, std::invalid_argument
    ,"MatrixComposite::add_matrix(...) : Error!" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    col_offset + opA_cols > cols_, std::invalid_argument
    ,"MatrixComposite::add_matrix(...) : Error!" );

  fully_constructed_ = false;

  matrix_list_.push_back(
    SubMatrixEntry(
      row_offset+1,row_offset+opA_rows
      ,col_offset+1,col_offset+opA_cols
      ,alpha
      ,Range1D()
      ,Teuchos::null,Teuchos::null,BLAS_Cpp::no_trans
      ,A,A_release,A_trans
      ,Range1D()
      ,Teuchos::null,Teuchos::null,BLAS_Cpp::no_trans
      )
    );
}

void MatrixComposite::add_matrix(
  size_type                      row_offset
  ,size_type                     col_offset
  ,value_type                    alpha
  ,const GenPermMatrixSlice      *P
  ,const release_resource_ptr_t  &P_release
  ,BLAS_Cpp::Transp              P_trans
  )
{
  namespace rcp = MemMngPack;
  using BLAS_Cpp::rows;
  using BLAS_Cpp::cols;

  TEUCHOS_TEST_FOR_EXCEPT( !(  alpha != 0.0  ) );
  TEUCHOS_TEST_FOR_EXCEPT( !(  P != NULL  ) );

  fully_constructed_ = false;

  const size_type
    P_rows   = P->rows(),
    P_cols   = P->cols(),
    opP_rows = rows(P_rows,P_cols,P_trans),
    opP_cols = cols(P_rows,P_cols,P_trans);

  TEUCHOS_TEST_FOR_EXCEPT( !(  row_offset + opP_rows <= rows_  ) );
  TEUCHOS_TEST_FOR_EXCEPT( !(  col_offset + opP_cols <= cols_  ) );

  matrix_list_.push_back(
    SubMatrixEntry(
      row_offset+1,row_offset+opP_rows,col_offset+1,col_offset+opP_cols,alpha
      ,Range1D::Invalid
      ,Teuchos::rcp(new GenPermMatrixSlice(*P)),P_release,P_trans
      ,NULL,Teuchos::null,BLAS_Cpp::no_trans
      ,Range1D()
      ,Teuchos::null,Teuchos::null,BLAS_Cpp::no_trans
      )
    );
}

void MatrixComposite::remove_matrix( matrix_list_t::iterator itr )
{
  fully_constructed_ = false;
  matrix_list_.erase(itr);
}

void MatrixComposite::finish_construction(
  const VectorSpace::space_ptr_t&  space_cols
  ,const VectorSpace::space_ptr_t& space_rows
  )
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !space_cols.get(), std::invalid_argument
    ,"MatrixComposite::finish_construction(...): Error, space_cols.get() can not be NULL" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    !space_rows.get(), std::invalid_argument
    ,"MatrixComposite::finish_construction(...): Error, space_rows.get() can not be NULL" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    space_cols->dim() != rows_, std::invalid_argument
    ,"MatrixComposite::finish_construction(...): Error, space_colss->dim() = " << space_cols->dim()
    << " != rows = " << rows_ << " where cols was passed to this->reinitialize(...)" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    space_rows->dim() != cols_, std::invalid_argument
    ,"MatrixComposite::finish_construction(...): Error, space_rows->dim() = " << space_rows->dim()
    << " != cols = " << cols_ << " where cols was passed to this->reinitialize(...)" );

  space_rows_ = space_rows;
  space_cols_ = space_cols;
  fully_constructed_ = true;
}

// Member access

int MatrixComposite::num_vectors() const
{
  return vector_list_.size();
}

MatrixComposite::vector_list_t::iterator
MatrixComposite::vectors_begin()
{
  return vector_list_.begin();
}

MatrixComposite::vector_list_t::iterator
MatrixComposite::vectors_end()
{
  return vector_list_.end();
}

MatrixComposite::vector_list_t::const_iterator
MatrixComposite::vectors_begin() const
{
  return vector_list_.begin();
}

MatrixComposite::vector_list_t::const_iterator
MatrixComposite::vectors_end() const
{
  return vector_list_.end();
}

int MatrixComposite::num_matrices() const
{
  return matrix_list_.size();
}

MatrixComposite::matrix_list_t::iterator
MatrixComposite::matrices_begin()
{
  return matrix_list_.begin();
}

MatrixComposite::matrix_list_t::iterator
MatrixComposite::matrices_end()
{
  return matrix_list_.end();
}

MatrixComposite::matrix_list_t::const_iterator
MatrixComposite::matrices_begin() const
{
  return matrix_list_.begin();
}

MatrixComposite::matrix_list_t::const_iterator
MatrixComposite::matrices_end() const
{
  return matrix_list_.end();
}

// Overridden from Matrix

size_type MatrixComposite::rows() const
{
  return fully_constructed_ ? rows_ : 0;
}

size_type MatrixComposite::cols() const
{
  return fully_constructed_ ? cols_ : 0;
}

size_type MatrixComposite::nz() const
{
  if(fully_constructed_) {
    size_type nz = 0;
    {for(matrix_list_t::const_iterator itr = matrix_list_.begin(); itr != matrix_list_.end(); ++itr) {
      if( itr->A_ != NULL )
        nz += itr->A_->nz();
      else
        nz += (*itr->P_).nz();
    }}
    {for(vector_list_t::const_iterator itr = vector_list_.begin(); itr != vector_list_.end(); ++itr) {
      nz += itr->v_->nz();
    }}
    // ToDo: Update the above code to consider DMatrixSlice and Range1D views
    return nz;
  }
  return 0;
}

// Overridden from MatrixOp

const VectorSpace& MatrixComposite::space_rows() const
{
  assert_fully_constructed();
  return *space_rows_;
}

const VectorSpace& MatrixComposite::space_cols() const
{
  assert_fully_constructed();
  return *space_cols_;
}

MatrixOp::mat_ptr_t
MatrixComposite::sub_view(const Range1D& row_rng, const Range1D& col_rng) const
{
  assert_fully_constructed();
  // For this initial implementation we will just look through the list of matrices
  // and find an exact match.  If we don't find it we will just return the result
  // from the default implementation of this method.

  // ToDo: Implement!

  // Just return the default sub-view
  return MatrixOp::sub_view(row_rng,col_rng);
}

void MatrixComposite::Vp_StMtV(
  VectorMutable* y, value_type a, BLAS_Cpp::Transp M_trans
  , const Vector& x, value_type b
  ) const
{
#ifdef PROFILE_HACK_ENABLED
  ProfileHackPack::ProfileTiming profile_timing( "MatrixComposite::Vp_StMtV(...DVectorSlice...)" );
#endif
  assert_fully_constructed();
  AbstractLinAlgPack::Vp_MtV_assert_compatibility( y, *this, M_trans, x );
  Vp_StMtV_imp(y,a,rows_,cols_,matrix_list_,vector_list_,M_trans,x,b);
}

void MatrixComposite::Vp_StMtV(
  VectorMutable* y, value_type a, BLAS_Cpp::Transp M_trans
  , const SpVectorSlice& x, value_type b
  ) const
{
#ifdef PROFILE_HACK_ENABLED
  ProfileHackPack::ProfileTiming profile_timing( "MatrixComposite::Vp_StMtV(...SpVectorSlice...)" );
#endif
  assert_fully_constructed();
  AbstractLinAlgPack::Vp_MtV_assert_compatibility( y, *this, M_trans, x );
  Vp_StMtV_imp(y,a,rows_,cols_,matrix_list_,vector_list_,M_trans,x,b);
}

void MatrixComposite::Vp_StPtMtV(
  VectorMutable* y, value_type a
  , const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
  , BLAS_Cpp::Transp M_trans
  , const Vector& x, value_type b
  ) const
{
#ifdef PROFILE_HACK_ENABLED
  ProfileHackPack::ProfileTiming profile_timing( "MatrixComposite::Vp_StPtMtV(...DVectorSlice...)" );
#endif
  assert_fully_constructed();
  MatrixOp::Vp_StPtMtV(y,a,P,P_trans,M_trans,x,b); // ToDo: Implement when needed!
}

void MatrixComposite::Vp_StPtMtV(
  VectorMutable* y, value_type a
  , const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
  , BLAS_Cpp::Transp M_trans
  , const SpVectorSlice& x, value_type b
  ) const
{	
#ifdef PROFILE_HACK_ENABLED
  ProfileHackPack::ProfileTiming profile_timing( "MatrixComposite::Vp_StPtMtV(...SpVectorSlice...)" );
#endif
  assert_fully_constructed();
  MatrixOp::Vp_StPtMtV(y,a,P,P_trans,M_trans,x,b); // ToDo: Implement when needed!
}

// private

void MatrixComposite::assert_fully_constructed() const
{
  const bool fully_constructed = fully_constructed_;
  TEUCHOS_TEST_FOR_EXCEPTION(
    !fully_constructed, std::logic_error
    ,"MatrixComposite::assert_fully_constructed() : Error, not fully constructed!");
}

} // end namespace AbstractLinAlgPack
