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

#ifndef MATRIX_SYM_POS_DEF_CHOL_FACTOR_H
#define MATRIX_SYM_POS_DEF_CHOL_FACTOR_H

#include "AbstractLinAlgPack_Types.hpp"
#include "AbstractLinAlgPack_MatrixExtractInvCholFactor.hpp"
#include "AbstractLinAlgPack_MatrixSymAddDelUpdateable.hpp"
#include "AbstractLinAlgPack_MatrixSymOpNonsingSerial.hpp"
#include "AbstractLinAlgPack_MatrixSymDenseInitialize.hpp"
#include "AbstractLinAlgPack_MatrixSymOpGetGMSSymMutable.hpp"
#include "AbstractLinAlgPack_MatrixSymSecant.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "DenseLinAlgPack_DMatrixAsTriSym.hpp"
#include "SerializationPack_Serializable.hpp"
#include "Teuchos_RCP.hpp"
#include "ReleaseResource.hpp"

namespace AbstractLinAlgPack {
/** \brief A do all class for dense symmetric positive definite matrices
 * that stores the original matrix and/or its upper cholesky factor.
 *
 * This matrix class supports a boat load of interfaces.  It is ment to be
 * a do all matrix class for all dense positive definite matrices stored as
 * the upper cholesky factor and/or the lower symmetric nonfactorized
 * portion.  It is designed to meet many different needs and apply in many different
 * contexts.  Objects from this class can represent matrices of the form:
 \verbatim
     M = scale * U' * U\endverbatim
 * where <tt>U</tt> is an upper triangular cholesky factor (i.e. the diagonal of
 * <tt>U</tt> is positive).  Also, the lower part of <tt>M</tt> may also be stored 
 * and manipulated.  The reason for maintaining the unfactorized matrix
 * also is to allow its use in contexts where the cholesky factor may not
 * be needed since linear systems are never solved for.  It is also useful
 * for (extended precision) iterative refinement.
 *
 * The purpose of <tt>scale</tt> in <tt>M = scale * U' * U</tt> is in order to allow
 * matrices of this type to represent positive definite (scale > 0)
 * or negative definite (scale < 0) matrices.
 *
 * This class allows you to do a bunch of stuff with these matrix objects:
 * \begin{itemize}
 * \item BLAS operations (\Ref{MatrixSymOp})
 * \item Solve for linear systems (\Ref{MatrixSymFactorized})
 * \item Initialize from a dense symmetric matrix, provided that the matrix
 *   is positive definite (\Ref{MatrixSymDenseInitialize}) (default implementation)
 * \item Extract a dense inverse cholesky factor (\Ref{MatrixExtractInvCholFactor})
 * \item Perform BFGS positive definite secant updates (\Ref{MatrixSymSecant})
 *       (default implementation)
 * \item Add rows/cols (provided new matrix is p.d. (scale > 0) or n.d. (scale < 0)) and
 *   remove rows/cols (\Ref{MatrixSymAddDelUpdateable}) and update the factors.
 *   Note that these operations will change the size of <tt>M</tt> and/or <tt>U</tt>.
 * \end{itemize}
 *
 * The <tt>DMatrixSliceTri</tt> view <tt>U</tt> is really a subview of another <tt>DMatrixSlice</tt>
 * <tt>MU_store</tt>.  The matrix <tt>U</tt> is defined as:
 *
 * <tt>U = tri(MU_store(U_l_r:U_l_r+M_size-1,U_l_c+1:U_l_c+M_size),upper,nonunit)</tt>.
 *
 * The <tt>DMatrixSliceSym</tt> view <tt>M</tt> is really another subview of <tt>MU_store</tt>:
 *
 * <tt>M = sym(MU_store(M_l_r+1:M_l_r+M_size,M_l_c:M_l_c+M_size-1),lower)</tt>.
 *
 * The reason for the offsets in the definition of <tt>M</tt> and <tt>U</tt> above is to keep the
 * diagonal of <tt>MU_store</tt> open for use as workspace.
 *
 * To allow for more flexible resourse management, this class maintains a
 * <tt>RCP<ReleaseResource></tt> object so whatever memory needs to be released
 * will be released when <tt>MU_store</tt> is no longer needed.  This class will also
 * allocate its own memory if none is set through the constructor or the
 * init_setup(...) method.
 *
 * In short, this class should be able to handle most uses for a dense symmetric
 * positive definite matrix stored as a cholesky factor.
 * Also, direct access is given to <tt>U</tt>, <tt>M</tt>, <tt>MU_store</tt> and the ability
 * to change the setup.  I don't know what more do you want!
 *
 * More operations will be overridden as they are needed by various applications.
 */
class MatrixSymPosDefCholFactor
  : virtual public AbstractLinAlgPack::MatrixSymOpNonsingSerial     // doxygen needs full name
  , virtual public AbstractLinAlgPack::MatrixSymDenseInitialize     // ""
  , virtual public AbstractLinAlgPack::MatrixSymOpGetGMSSymMutable  // ""
  , virtual public MatrixExtractInvCholFactor
  , virtual public MatrixSymSecant
  , virtual public MatrixSymAddDelUpdateable
  , virtual public SerializationPack::Serializable
{
public:
  
  /** @name Public types */
  //@{

  /** \brief . */
  typedef Teuchos::RCP<
    MemMngPack::ReleaseResource>  release_resource_ptr_t;

  /** \brief PostMod class to use with <tt>MemMngPack::AbstractFactorStd</tt>.
   */
  class PostMod {
  public:
    /** \brief . */
    PostMod(
      bool   maintain_original = true
      ,bool  maintain_factor   = false
      ,bool  allow_factor      = true
      )
      :maintain_original_(maintain_original)
      ,maintain_factor_(maintain_factor)
      ,allow_factor_(allow_factor)
    {}
    /** \brief . */
    void initialize(MatrixSymPosDefCholFactor* p) const
    {
      p->init_setup(
        NULL,Teuchos::null,0 // Allocate own storage
        ,maintain_original_,maintain_factor_,allow_factor_
        );
    }
  private:
    bool  maintain_original_;
    bool  maintain_factor_;
    bool  allow_factor_;
  }; // end PostMod

  //@}

    /** @name Constructors/initalizers */
  //@{

  /** \brief Initialize to uninitialized.
   *
   * By default if init_setup(...) is not called then this object
   * will allocate its own storage and maintain_original == true
   * and maintain_factor = false.
   */
  MatrixSymPosDefCholFactor();

  /** \brief Initialize with possible storage.
   *
   * This constructor just calls <tt>this->init_setup(...)</tt>.
   */
  MatrixSymPosDefCholFactor(
    DMatrixSlice                    *MU_store
    ,const release_resource_ptr_t&    release_resource_ptr = Teuchos::null
    ,size_type                        max_size             = 0
    ,bool                             maintain_original    = true
    ,bool                             maintain_factor      = false
    ,bool                             allow_factor         = true
    ,bool                             set_full_view        = true
    ,value_type                       scale                = 1.0
    );

  /** \brief Initial storage setup and possibly the view.
   *
   * @param  MU_store [in/state]
   *                  If <tt>MU_store != NULL</tt> then this matrix is 
   *                  used to store the original matrix and/or its
   *                  cholesky factor.  This matrix may or may not be
   *                  initialized and ready to go.  The maxinum size
   *                  for <tt>M</tt> that can be stored is <tt>max_size =</tt>
   *                  <tt>min(MU_store.rows(),MU_store.cols())-1</tt>.  The reason
   *                  for this is that the diagonal of <tt>MU_store</tt> is used
   *                  as workspace. If <tt>maintain_original == true</tt> then any portion
   *                  of the lower triangular region below the center
   *                  diagonal is open game to be accessed.  If
   *                  <tt>allow_factor == true</tt> then any portion of the
   *                  center diagonal or the upper triangular region
   *                  above the diagonal is fair game to be accessed.
   *                  If <tt>MU_store == NULL</tt> then any current
   *                  storage is unset and <tt>this</tt> matrix object is then
   *                  free to allocate its own storage as needed.
   * @param  release_resource_ptr
   *                  [in] Only significant if <tt>MU_store != NULL</tt>.
   *                  Points to a resource to be released when
   *                  <tt>MU_store</tt> is no longer needed.
   * @param  max_size [in] Only significant if<tt>MU_store == NULL</tt>.
   *                  If <tt>max_size > 0</tt> then this is the maximum size the matrix will be sized to.
   *                  If <tt>max_size == 0</tt> then the max size will be determined by the other initialization
   *                  functions.
   * @param  maintain_original
   *                  [in] Always significant.  If true then the original
   *                  matrix will be maintained in the lower triangle of
   *                  <tt>this->MU_store()</tt> (see intro).
   * @param  maintain_factor
   *                  [in] Always significant.  If true then the cholesky
   *                  factor will be maintained in the upper triangle of
   *                  <tt>this->MU_store()</tt> (see intro).
   * @param  allow_factor
   *                  [in] Only significant if <tt>maintain_factor == false</tt>.
   *                  If true then the factorization can be
   *                  computed in the upper triangular portion
   *                  of <tt>MU_store</tt>.  Otherwise it can not.
   * @param  set_full_view
   *                  [in] Only significant if <tt>MU_store != NULL</tt>.
   *                  If true then <tt>this</tt> will be setup as the
   *                  full view of <tt>MU_store</tt> and <tt>MU_store</tt> should
   *                  already be initialized properly.
   * @param  scale    [in] Only significant if <tt>MU_store != NULL</tt> and <tt>set_full_view == true</tt>
   *                  (see intro).
   *
   * Postconditions:<ul>
   * <li> [<tt>MU_store!=NULL && set_full_view==true</tt>]
   *      <tt>this->rows() == min( MU_store->rows(), MU_store->cols() ) - 1</tt>
   * </ul>
   */
  void init_setup(
    DMatrixSlice                    *MU_store
    ,const release_resource_ptr_t&    release_resource_ptr = Teuchos::null
    ,size_type                        max_size             = 0
    ,bool                             maintain_original    = true
    ,bool                             maintain_factor      = false
    ,bool                             allow_factor         = true
    ,bool                             set_full_view        = true
    ,value_type                       scale                = 1.0
    );
  /** \brief Resize the view <tt>M</tt> and <tt>U</tt> of <tt>MU_store</tt>.
   *
   * Preconditions:\begin{itemize}
   * \item [<tt>!allocates_storage()</tt>] <tt>MU_store().rows() > 0</tt>
   * \item <tt>M_size >= 0</tt>
   * \item [<tt>M_size > 1</tt>] <tt>scale != 0.0</tt>
   * \item <tt>maintain_original || maintain_factor</tt>
   * \item [<tt>maintain_original</tt>] <tt>1 <= M_l_r <= M_l_c</tt>
   * \item [<tt>maintain_original</tt>] <tt>M_l_r + M_size <= MU_store_.rows()</tt>
   * \item <tt>U_l_r >= U_l_c</tt>
   * \item [<tt>U_l_r > 0</tt>] <tt>U_l_c + M_size <= MU_store_.cols()</tt>
   * \end{itemize}
   *
   * @param  M_size  [in] Size of <tt>M</tt> (see intro).  If <tt>M_size == 0</tt> then the
   *                 matrix will be set to a size of zero.
   * @param  scale   [in] Only significant if <tt>M_size > 0</tt>.  If <tt>scale > 0</tt> then <tt>M</tt> must
   *                 be p.d. and if <tt>scale < 0</tt> then n.d.
   * @param  maintain_original
   *                 [in] Always significant.
   *                 If true then original <tt>M</tt> is maintained and also
   *                 elements in <tt>this->M()</tt> are expected to be initialized
   *                 already.  If false then the lower triangular portion of
   *                 <tt>MU_store()</tt> is strictly off limits and will never be
   *                 touched.
   * @param  M_l_r   [in] Only significant if <tt>M_size > 0</tt>.  Lower row index for <tt>M</tt> (see intro)
   * @param  M_l_c   [in] Only significant if <tt>M_size > 0</tt>.  Lower column index for <tt>M</tt> (see intro)
   * @param  maintain_factor
   *                 [in] Always significant.
   *                 If true then the factor <tt>U</tt> is maintained and also
   *                 the elements in <tt>this->U()</tt> are expected to be initialized
   *                 already.  If false then the factor may still be computed
   *                 and stored in the upper triangular part of <tt>MU_store()</tt>
   *                 if needed by members (such as <tt>V_InvMtV(...)</tt>) but only if
   *                 <tt>U_l_r > 0</tt>.
   * @param  U_l_r   [in] If <tt>M_size == 0</tt> then only <tt>U_l_r == 0</tt> and <tt>U_l_r > 0</tt> is significant.
   *                 Otherwise <tt>U_l_r</tt> is the lower row index for <tt>U</tt> (see intro).
   *                 If <tt>U_l_r == 0</tt> then the upper triangular portion of <tt>MU_store()</tt> is strictly
   *                 off limits and the factorization will never be computed.
   *                 Any functions that need this factorization will throw
   *                 exceptions in these cases.
   * @param  U_l_c   [in] Only significant if <tt>U_l_r > 0</tt>.
   *                  Lower column index for <tt>U</tt> (see intro).
   */
  void set_view(
    size_t            M_size
    ,value_type       scale
    ,bool             maintain_original
    ,size_t           M_l_r
    ,size_t           M_l_c
    ,bool             maintain_factor
    ,size_t           U_l_r
    ,size_t           U_l_c
    );
  /** \brief Set the default pivot tolerance.
   */
  void pivot_tols( PivotTolerances pivot_tols );
  /** \brief . */
  PivotTolerances pivot_tols() const;

  //@}

    /**  @name Access representation */
  //@{

  /** \brief Get the current setup.
   */
  void get_view_setup(
    size_t            *M_size
    ,value_type       *scale
    ,bool             *maintain_original
    ,size_t           *M_l_r
    ,size_t           *M_l_c
    ,bool             *maintain_factor
    ,size_t           *U_l_r
    ,size_t           *U_l_c
    ) const;
  /** \brief Return if this object owns and allocates storage.
   */
  bool allocates_storage() const;

  /** \brief Get access to MU_store.
   */
  DMatrixSlice& MU_store();
  /** \brief . */
  const DMatrixSlice& MU_store() const;
  /** \brief Get view of U.
   */
  DMatrixSliceTri U();
  /** \brief . */
  const DMatrixSliceTri U() const;
  /** \brief Get view of lower part of M.
   */
  DMatrixSliceSym M();
  /** \brief . */
  const DMatrixSliceSym M() const;

  //@}

  /** @name Overridden from MatrixBase */
  //@{

  /** \brief . */
  size_type rows() const;

  //@}

  /** @name Overridden from MatrixOp */
  //@{

  /** \brief . */
  void zero_out();
  /** \brief . */
  std::ostream& output(std::ostream& out) const;
  /** \brief . */
  bool Mp_StM(
    MatrixOp* m_lhs, value_type alpha
    ,BLAS_Cpp::Transp trans_rhs
    ) const;
  /** \brief . */
  bool Mp_StM(
    value_type alpha,const MatrixOp& M_rhs, BLAS_Cpp::Transp trans_rhs
    );
  /** \brief . */
  bool syrk(
    const MatrixOp      &mwo_rhs
    ,BLAS_Cpp::Transp   M_trans
    ,value_type         alpha
    ,value_type         beta
    );

  //@}

  /** @name Overridden from MatrixOpSerial */
  //@{

  /** \brief . */
  void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const DVectorSlice& vs_rhs2, value_type beta) const;
  /** \brief . */
  void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const SpVectorSlice& vs_rhs2, value_type beta) const;
  /** \brief . */
  void Vp_StPtMtV(DVectorSlice* vs_lhs, value_type alpha
    , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
    , BLAS_Cpp::Transp M_rhs2_trans
    , const DVectorSlice& vs_rhs3, value_type beta) const;
  /** \brief . */
  void Vp_StPtMtV(DVectorSlice* vs_lhs, value_type alpha
    , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
    , BLAS_Cpp::Transp M_rhs2_trans
    , const SpVectorSlice& sv_rhs3, value_type beta) const;

  //@}

  /** @name Overridden form MatrixSymOpSerial */
  //@{

  void Mp_StPtMtP( DMatrixSliceSym* sym_lhs, value_type alpha
    , EMatRhsPlaceHolder dummy_place_holder
    , const GenPermMatrixSlice& gpms_rhs, BLAS_Cpp::Transp gpms_rhs_trans
    , value_type beta ) const;

  //@}

  /** @name Overridden from MatrixNonsingSerial */
  //@{

  /// With throw exception if factorization is not allowed.
  void V_InvMtV(DVectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1
    , const DVectorSlice& vs_rhs2) const;
  /// With throw exception if factorization is not allowed.
  void V_InvMtV(DVectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1
    , const SpVectorSlice& sv_rhs2) const;

  //@}

  /** @name Overridden from MatrixSymNonsingSerial */
  //@{

  /// Will throw exception if factorization is not allowed.
  void M_StMtInvMtM(
    DMatrixSliceSym* sym_gms_lhs, value_type alpha
    , const MatrixOpSerial& mwo, BLAS_Cpp::Transp mwo_trans, EMatrixDummyArg
    ) const;

  //@}

  /** @name Overridden from MatrixSymDenseInitialize */
  //@{

  /// Will resize view of matrices and reset scale
  void initialize( const DMatrixSliceSym& M );

  //@}

  /** @name Overridden from MatrixSymOpGetGMSSym */
  //@{

  /** \brief . */
  const DenseLinAlgPack::DMatrixSliceSym get_sym_gms_view() const;
  /** \brief . */
  void free_sym_gms_view(const DenseLinAlgPack::DMatrixSliceSym* sym_gms_view) const;

  //@}

  /** @name Overridden from MatrixSymOpGetGMSSymMutable */
  //@{

  /** \brief . */
  DenseLinAlgPack::DMatrixSliceSym get_sym_gms_view();
  /** \brief . */
  void commit_sym_gms_view(DenseLinAlgPack::DMatrixSliceSym* sym_gms_view);

  //@}

  /** @name Overridden from MatrixExtractInvCholFactor */
  //@{

  /** \brief . */
  void extract_inv_chol( DMatrixSliceTriEle* InvChol ) const;
  
  //@}

  /** @name Overridden from MatrixSymSecantUpdateble */
  //@{

  /// Will reset view and set scale
  void init_identity( const VectorSpace& space_diag, value_type alpha );
  /// Will reset view and set scale
  void init_diagonal( const Vector& diag );
  /// Must agree with current scale
  void secant_update(VectorMutable* s, VectorMutable* y, VectorMutable* Bs);

  //@}

  /** @name Overridden from MatrixSymAddDelUpdateble */
  //@{

  /// Will reset view of U and M and reset scale
  void initialize(
    value_type         alpha
    ,size_type         max_size
    );
  /// Will reset view of U and M and reset scale
  void initialize(
    const DMatrixSliceSym      &A
    ,size_type         max_size
    ,bool              force_factorization
    ,Inertia           inertia
    ,PivotTolerances   pivot_tols
    );
  /** \brief . */
  size_type max_size() const;
  /// Will be (rows(),0,0) if scale < 0 or (0,0,rows()) if scale > 0.
  Inertia inertia() const;
  /// Will set rows() == 0
  void set_uninitialized();
  /// Will throw exceptions if not p.d. (scale > 0) or n.d. (scale < 0).
  void augment_update(
    const DVectorSlice  *t
    ,value_type        alpha
    ,bool              force_refactorization
    ,EEigenValType     add_eigen_val
    ,PivotTolerances   pivot_tols
    );
  /// Should always succeed unless user gives wrong value for drop_eigen_val
  void delete_update(
    size_type          jd
    ,bool              force_refactorization
    ,EEigenValType     drop_eigen_val
    ,PivotTolerances   pivot_tols
    );

  //@}


  /** @name Overridden from Serializable */
  //@{

  /** \brief . */
  void serialize( std::ostream &out ) const;
  /** \brief . */
  void unserialize( std::istream &in );

  //@}

private:
  
  // /////////////////////////////
  // Private data members

  bool                    maintain_original_;
  bool                    maintain_factor_;
  bool                    factor_is_updated_;    // Is set to true if maintain_factor_=true
  bool                    allocates_storage_;    // If true then this object allocates the storage
  release_resource_ptr_t  release_resource_ptr_;
  DMatrixSlice            MU_store_;
  size_t                  max_size_;
  size_t                  M_size_,               // M_size == 0 is flag that we are completely uninitialized
                          M_l_r_,
                          M_l_c_,
                          U_l_r_,
                          U_l_c_;
  value_type              scale_;
  bool                    is_diagonal_;
  PivotTolerances         pivot_tols_;
  DVector                 work_; // workspace.

  // /////////////////////////////
  // Private member functions

  void assert_storage() const;
  void allocate_storage(size_type max_size) const;
  void assert_initialized() const;
  void resize_and_zero_off_diagonal(size_type n, value_type scale);
  void update_factorization() const;
  std::string build_serialization_string() const;
  static void write_matrix( const DMatrixSlice &Q, BLAS_Cpp::Uplo Q_uplo, std::ostream &out );
  static void read_matrix( std::istream &in, BLAS_Cpp::Uplo Q_uplo, DMatrixSlice *Q );

}; // end class MatrixSymPosDefCholFactor

// ///////////////////////////////////////////////////////
// Inline members for MatrixSymPosDefCholFactor

inline
bool MatrixSymPosDefCholFactor::allocates_storage() const
{
  return allocates_storage_;
}

inline
DMatrixSlice& MatrixSymPosDefCholFactor::MU_store()
{
  return MU_store_;
}

inline
const DMatrixSlice& MatrixSymPosDefCholFactor::MU_store() const
{
  return MU_store_;
}

inline
void MatrixSymPosDefCholFactor::get_view_setup(
  size_t            *M_size
  ,value_type       *scale
  ,bool             *maintain_original
  ,size_t           *M_l_r
  ,size_t           *M_l_c
  ,bool             *maintain_factor
  ,size_t           *U_l_r
  ,size_t           *U_l_c
  ) const
{
  *M_size               = M_size_;
  *scale                = scale_;
  *maintain_original    = maintain_original_;
  *M_l_r                = maintain_original_ ? M_l_r_ : 0;
  *M_l_c                = maintain_original_ ? M_l_c_ : 0;
  *maintain_factor      = maintain_factor_;
  *U_l_r                = maintain_factor_ ? U_l_r_ : 0;
  *U_l_c                = maintain_factor_ ? U_l_c_ : 0;
}

inline
DMatrixSliceTri MatrixSymPosDefCholFactor::U()
{
  return DenseLinAlgPack::nonconst_tri(
    MU_store_(U_l_r_,U_l_r_+M_size_-1,U_l_c_+1,U_l_c_+M_size_)
    , BLAS_Cpp::upper, BLAS_Cpp::nonunit
    );
}

inline
const DMatrixSliceTri MatrixSymPosDefCholFactor::U() const
{
  return DenseLinAlgPack::tri(
    MU_store_(U_l_r_,U_l_r_+M_size_-1,U_l_c_+1,U_l_c_+M_size_)
    , BLAS_Cpp::upper, BLAS_Cpp::nonunit
    );
}

inline
DMatrixSliceSym MatrixSymPosDefCholFactor::M()
{
  return DenseLinAlgPack::nonconst_sym(
    MU_store_(M_l_r_+1,M_l_r_+M_size_,M_l_c_,M_l_c_+M_size_-1)
    , BLAS_Cpp::lower
    );
}

inline
const DMatrixSliceSym MatrixSymPosDefCholFactor::M() const
{
  return DenseLinAlgPack::sym(
    MU_store_(M_l_r_+1,M_l_r_+M_size_,M_l_c_,M_l_c_+M_size_-1)
    , BLAS_Cpp::lower
    );
}

} // end namespace AbstractLinAlgPack

#endif // MATRIX_SYM_POS_DEF_CHOL_FACTOR_H
