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

#ifndef	ALAP_DIRECT_SPARSE_SOLVER_H
#define ALAP_DIRECT_SPARSE_SOLVER_H

#include "AbstractLinAlgPack_Types.hpp"
#include "AbstractLinAlgPack_MatrixConvertToSparse.hpp"
#include "AbstractLinAlgPack_MatrixNonsing.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_AbstractFactory.hpp"

namespace AbstractLinAlgPack {

/** \brief Abstract interface to serial direct sparse linear solvers.
 *
 * This interface is designed to accommodate both unsymmetic and
 * symmetric direct solvers.  The motivation for using one common
 * interface is that we would like to be able to easily use an
 * unsymmetic solver on a symmetric system.  Of couse this would also
 * allow a client to attempt to solve an unsymmetric system with a
 * symmetric solver which of course will not work.  A symmetric
 * implementation can of course check that <tt>A.rows()==A.cols()</tt>
 * but it would be much more difficult to determine if the matrix was
 * indeed symmetric.  It is likely that the symmetric solver would
 * only extract the upper or lower triangular region of an unsymmetric
 * matrix and would never check the other triangular region for
 * compatibility.
 *
 * The interface is designed to allow for maximum memory sharing
 * and recycling.  For example, we would like for the client to be able
 * to maintain more than one set of nonzero factors using the same sparsity
 * and fill-in structure determined by the \c analyze_and_factor()
 * method.
 *
 * The basic idea is that given an <tt>m x n</tt> rectangular matrix \c A,
 * where <tt>m <= n</tt>, objects of this type will find
 * row and column permutations \c P and \c Q and a rank \c r such that
 * <tt>C = (P'*A*Q)(1:r,1:r)</tt> is the largest square nonsingular
 * (well conditioned) matrix that can be found.  This basis matrix
 * \c C is represented as a \c BasisMatrix object that supports the
 * \c MatrixNonsingSerial interface.
 *
 * All the information needed to solve for a linear
 * system and to factor other matrices with the same structure
 * is contained in the \c BasisMatrix object.  Therefore,
 * a \c DirectSparseSolver (\c DSS) object can
 * be considered to be stateless if needed.  The first point to make
 * clear is that once the client has a reference to a \c BasisMatrix object
 * in a smart pointer like \c C1_ptr, that \c BasisMatrix object will not be altered
 * by any future operations on the \c DSS object that created it unless the
 * client specifically requests an alteration.  This behavior allows \c BasisMatrix
 * objects to be completely decoupled from the \c DSS object that created them.
 * This behavior is explained in more detail later.
 *
 * It is through the \c MatrixNonsing interface of the \c BasisMatrix object
 * that clients can solve for linear systems.
 *
 * The usage of this interface will be illustrated by considering
 * a few different scenarios:
 *
 * 1) Maintain one factorization structure of an unsymmetric matrix.
 *
 \code

  // Storage for the permutations and rank
  IVector row_perm, col_perm;
  size_type rank;
  
  // Analyze and factor A1
  typedef DirectSparseSolver::basis_matrix_ptr_t C_ptr_t;
  C_ptr_t  C1_ptr = direct_solver.basis_matrix_factory()->create();
  direct_solver.analyze_and_factor( A1, &row_perm, &col_perm, &rank, C1_ptr.get() );
  // After the above method finishes direct_solver.get_fact_struc().get() points
  // to the same factorization structure as C1_ptr->get_fact_struc().get().

  // Solve for x1 = inv(C1)*b1
  VectorMutableDense x1(C1_ptr->cols());
  V_InvMtV( &x1, *C1_ptr, BLAS_Cpp::no_trans, b1 );

  // Factor another matrix with the same structure.
  // The current factorization structure in direct_solver.get_fact_struc() is used.
  C_ptr_t  C2_ptr = direct_solver.basis_matrix_factory()->create();
  direct_solver.factor( A2, C2_ptr.get() );

  // Solve for x2 = inv(C2')*b2
  VectorMutableDense x2(C2_ptr->rows());
  V_InvMtV( &x2, *C2_ptr, BLAS_Cpp::trans, b2 );

 \endcode
 *
 * 2) Maintain the factorization of three unsymmetic matrices, two with the
 * same structure and one with a different structure.
 *
 \code
  
  // Storage for the permutations and rank
  IVector row_perm1, col_perm1;
  size_type rank1;
  IVector row_perm2, col_perm2;
  size_type rank2;

  // Analyze and factor A1
  typedef DirectSparseSolver::basis_matrix_ptr_t C_ptr_t;
  C_ptr_t  C1_ptr = direct_solver.basis_matrix_factory()->create();
  direct_solver.analize_and_factor( A1, &row_perm1, &col_perm1, &rank1, C1_ptr.get() );

  // Solve for x1 = inv(C1)*b1
  VectorMutableDense x1(C1_ptr->cols());
  V_InvMtV( &x1, *C1_ptr, BLAS_Cpp::no_trans, b1 );

  // Factor another matrix A2 with a different structure from A1.
  C_ptr_t  C2_ptr = direct_solver.basis_matrix_factory()->create();
  direct_solver.analize_and_factor( A2, &row_perm2, &col_perm2, &rank2, C2_ptr.get() );
  // In the above method, a new current factorization structure will be computed and
  // used for C2.  This is because the current factorization structure before the
  // method is called is being used by the basis matrix object *C1_ptr and must
  // be preserved.

  // Solve for x2 = inv(C2)*b2
  VectorMutableDense x2(C2_ptr->cols());
  V_InvMtV( &x2, *C2_ptr, BLAS_Cpp::no_trans, b2 );

  // Factor another matrix A3 with the same structure of as A1 but preserve
  // the factorization nonzeros of A1.
  C_ptr_t  C3_ptr = direct_solver.basis_matrix_factory()->create();
  direct_solver.factor( A3, C3_ptr.get(), C1_ptr->get_fact_struc() );
  // Above, the current factorization structure direct_solver.get_fact_struc()
  // will be preserved before and after this method call.

  // Solve for x3 = inv(C3)*b3
  VectorMutableDense x3(C3_ptr->cols());
  V_InvMtV( &x3, *C3_ptr, BLAS_Cpp::no_trans, b3 );

 \endcode
 *
 * 3) Factor of two unsymmetic matrices with different structure and then
 * recycle the storage of the first matrix to factor a third matrix with
 * an all new structure.
 *
 \code
  
  // Storage for the permutations and rank
  IVector row_perm1, col_perm1;
  size_type rank1;
  IVector row_perm2, col_perm2;
  size_type rank2;

  typedef DirectSparseSolver::basis_matrix_ptr_t  C_ptr_t;

  // Analyze and factor A1
  C_ptr_t  C1_ptr = direct_solver.basis_matrix_factory()->create();
  direct_solver.analize_and_factor( A1, &row_perm1, &col_perm1, &rank1, C1_ptr.get() );

  // Solve for x1 = inv(C1)*b1
  VectorMutableDense x1(C1_ptr->cols());
  V_InvMtV( &x1, *C1_ptr, BLAS_Cpp::no_trans, b1 );

  // Factor another matrix A2 with a different structure from A1.
  C_ptr_t  C2_ptr = direct_solver.basis_matrix_factory()->create();
  direct_solver.analize_and_factor( A2, &row_perm2, &col_perm2, &rank2, C2_ptr.get() );

  // Solve for x2 = inv(C2)*b2
  VectorMutableDense x2(C2_ptr->cols());
  V_InvMtV( &x2, *C2_ptr, BLAS_Cpp::no_trans, b2 );

  // Analyze and factor another matrix A3 with an all new structure and recycle storage of C1.
  C_ptr_t  C3_ptr = C1_ptr;
  C1_ptr = Teuchos::null;
  direct_solver.analyze_and_factor( A3, &row_perm1, &col_perm1, &rank1, C3_ptr.get() );

  // Solve for x3 = inv(C3)*b3
  VectorMutableDense x3(C1_ptr->cols());
  V_InvMtV( &x3, *C3_ptr, BLAS_Cpp::no_trans, b3 );

  // Factor another matrix A4 with the same structure as A2 but preserve
  // its factorization in the basis matrix C2.
  C_ptr_t  C4_ptr = direct_solver.basis_matrix_factory()->create();
  direct_solver.factor( A4, C4_ptr.get(), C2_ptr->get_fact_struc() );

  // Solve for x4 = inv(C4)*b4
  VectorMutableDense x4(C4_ptr->cols());
  V_InvMtV( &x4, *C4_ptr, BLAS_Cpp::no_trans, b4 );

  // Analyze and factor another matrix A5 with an all new structure
  // and recycle nonzero storage of C4.  Note that the behavoir of C2_ptr must
  // not change by this operation.  Therefore, a new factorization structure
  // must be allocated in this case.
  C_ptr_t  C5_ptr = C4_ptr;
  C4_ptr = Teuchos::null;
  direct_solver.analyze_and_factor( A5, &row_perm5, &col_perm5, &rank5, C4_ptr.get() );

 \endcode
 *
 * If by some major mistake the client were to pass in \c basis_matrix
 * or \c fact_struc objects of the incorrect type from non-compatible
 * DSS objects then the \c InvalidObjectType exception would be thrown.
 * If the client follows the protocal defined in this interface, this will
 * never happen.
 */
class DirectSparseSolver {
public:

  /** @name Public Types */
  //@{

  /** \brief Abstract class for objects that represent the factorization
    * structure of a particular matrix.
    *
    * This structure can be reused over and over again for factorizing matrices
    * with the same structure but different nonzero elements.
    */
  class FactorizationStructure {
  public:
    /** \brief . */
    virtual ~FactorizationStructure() {}
  };

  /** \brief Abstract class for objects that represent the factorized
    * matrix and can be used to solve for different right-hand-sides.
    *
    * This object encapsulates the factorzation structure and the nonzero
    * values of the factorized basis matrix.
    */
  class BasisMatrix : public AbstractLinAlgPack::MatrixNonsing {
  public:
    /** \brief . */
    typedef Teuchos::RCP<FactorizationStructure>  fact_struc_ptr_t;
    /** \brief Return a reference to a smart pointer to the object that represents
     * the factorization structure.
     *
     * Returning a reference to a \c RCP<> object verses returning
     * a \c RCP<> object itself is critical so that we can rely on
     * \c RCP<>::count() to tell us how many clients have a reference
     * to this object.
     */
    virtual const fact_struc_ptr_t&  get_fact_struc() const = 0;
  };

  /** \brief . */
  typedef Teuchos::RCP<
    const Teuchos::AbstractFactory<BasisMatrix> >   basis_matrix_factory_ptr_t;

  /** \brief . */
  class UnsymmetricRankDeficientException : public std::logic_error
  {public: UnsymmetricRankDeficientException (const std::string& what_arg)
  : std::logic_error(what_arg) {}};

  /** \brief . */
  class IncompatibleMatrixStructureException : public std::logic_error
  {public: IncompatibleMatrixStructureException (const std::string& what_arg)
  : std::logic_error(what_arg) {}};

  /** \brief . */
  class InvalidObjectType : public std::logic_error
  {public: InvalidObjectType (const std::string& what_arg)
  : std::logic_error(what_arg) {}};

  /** \brief . */
  class NoCurrentBasisException : public std::logic_error
  {public: NoCurrentBasisException (const std::string& what_arg)
  : std::logic_error(what_arg) {}};

  /** \brief . */
  class FactorizationFailure : public std::logic_error
  {public: FactorizationFailure (const std::string& what_arg)
  : std::logic_error(what_arg) {}};

  //@}
  
  /** \brief . */
  virtual ~DirectSparseSolver() {}

  /** \brief Return a factory object that can create the basis matrix.
   */
  virtual const basis_matrix_factory_ptr_t basis_matrix_factory() const = 0;

  /** \brief Set the estimate of the fill-in ratio of the factorization.
   *
   * This value is used to set the initial guess at the amount
   * of storage needed for the factorization of the matrix
   * plus some "elbow" room for the factorization algorithm.
   *
   * The estimated amount of nonzeros for the factorization is:
   *
   * num_factor_nonzeros = estimated_fillin_ratio * A.nz()
   */
  virtual void estimated_fillin_ratio( value_type estimated_fillin_ratio ) = 0;

  /** \brief Analyze and factor a matrix.
   *
   * @param  A         [in] Matrix who's elements are extracted in order
   *                   to form the factorization.  Here <tt>A.rows() <= A.cols()</tt>
   *                   must be \c true.
   * @param  row_perm  [out] On output holds an array (size \c A.rows()) of row
   *                   permutations \c P used to select the basis matrix \c C.
   *                   The ith row of \c P'*A*Q is row \c (*row_perm)(i) of \c A.
   * @param  col_perm  [out] On output holds an array (size \c A.cols()) of column
   *                   permutations \c Q used to select the basis matrix \c C.
   *                   The jth column of \c P'*A*Q is column \c (*col_perm)(j) of \c A.
   *                   If the client is solving a symmetric system then
   *                   \c col_perm may be set to \c NULL on input which forces
   *                   the implementation to use symmetric permutations
   *                   (see above discussion).
   * @param  rank      [out] On output holds the rank \c r (<tt>r <= A.rows()</tt>)
   *                   of the basis matrix \c C.
   * @param  basis_matrix
   *                   [in/out] On input, this can be a previously initialized basis matrix
   *                   object in which case its storage for the structure of the factorization
   *                   and its nonzeros can be recycled if possible.
   * @param  out       [in/out] If <tt>out!=NULL</tt>, then some information about the operations performed
   *                   internally may be printed to \c *out.  The amount of this output should be
   *                   very minimal and should not significantly increase with the size of the problem
   *                   being solved.
   *
   * Preconditions:<ul>
   * <li> <tt>A.rows() <= A.cols()</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>basis_matrix != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditions: <ul>
   * <li> \c *row_perm, \c *col_perm and \c *rank give the basis selection determined by the
   *      direct sparse solver implementation.
   * <li> <tt>this->get_fact_struc().get() == basis_matrix->get_fact_struc().get()</tt> points
   *       to the factorizztion structure for this matrix.
   * </ul>
   *
   * This function extracts the nonzero elements from the matrix \c A (using the
   * \c AbstractLinAlgPack::MatrixConvertToSparse interface) and then finds a nonsingular
   * basis matrix and factors it (see intro).
   *
   * Given the <tt>m x n</tt> matrix \c A (<tt>m <= n</tt>), this method finds the
   * row and column permutations \c P and \c Q respectively and the rank \c r such that
   * the basis matrix <tt>C = (P'*A*Q)(1:r,1:r)</tt> is nonsingular.
   *
   * If a symmetric system is being solved and it is feared that \c A is not full rank
   * and the client wants to use symmetric pivoting (<tt>P = Q</tt>) to find the basis
   * \c C, then the client can pass in \c NULL for \c col_perm and force the solver to
   * return a symmetric permutation.  Note that the row and column permutations of
   * \c row_perm(k) and \c col_perm(k) for <tt>k = 1...r</tt> does not neccesarly indicate
   * the exact permutations used within the sparse solver.  Instead, what is significant is
   * the partitioning of rows and columns between <tt>k = 1..r</tt> and <tt>k = r+1..m</tt>
   * for row permutations and <tt>k = r+1..n</tt> for column permutations.  Therefore,
   * if <tt>col_perm==NULL</tt> on input and an unsymmetric solver is used and the matrix
   * is not full rank, then since in reality unsymmetric pivoting is being used internally,
   * the implementation must throw an \c UnsymmetricRankDeficientException exception.
   * This convention is needed to allow an unsymmetric solver to be used for a symmetric
   * system but still meet the needs of users for symmetric systems.  In summary, as long
   * as the symmetric system is full rank, then an unsymmetic solver can always be used and
   * the client need not care how the system is solved.  On the other hand, if the client does
   * not care if unsymmetric pivoting is used on a symmetric rank deficient system, then
   * the client can set \c col_perm to point to a valid \c IVector object and deal
   * with the unsymmetic permutations that are returned.  Who knows what the use of this
   * would be but there is no reason not to allow it.
   *
   * This method allows for the careful sharing and recycling of the factorization structure
   * and nonzeros of a factorized matrix while also maintaining a "current" factorization
   * structure.  To explain how this works let \c C1_ptr be set to the smart pointer to a basis
   * matrix object created by \c this->basis_matrix_factory()->create().
   * 
   * If the client whats to recycle the factorization structure and nonzeros that are
   * contained in \c *C1_ptr to analyze and factor \c A3, then the client would perform:
   *
   \code
   direct_solver.analyze_and_factor(A3,&row_perm,&col_perm,&rank,C1_ptr.get());
   \endcode
   *
   * If the factorization can not be performed for some reason then the exception
   * \c FactorizationFailure will be thrown.
   *
   * Scenarios / behavior:
   * <ul>
   * <li> \c basis_matrix->get_fact_struc().get()==NULL
   *      <ul>
   *      <li> \c this->get_fact_struc().get()==NULL
   *           <ul><li> New storage for the factorization structure will have to be
   *                    allocated. </ul>
   *      <li> \c this->get_fact_struc().get()!=NULL
   *           <ul>
   *           <li> \c this->get_fact_struc().count()==1
   *                <ul><li> No other basis matrix object is using this factorization structure
   *                 so the storage in this object can be recycled and used here.</ul>
   *           <li> \c this->get_fact_struc().count()>1
   *                <ul><li> Some other basis matrix object or some other client is using
   *                this factorization structure so new storage for the factorization
   *                structure must be allocated.</ul>
   *           </ul>
   *      <li> New storage for factorization <em>nonzeros</em> will have to be allocated.
   *      </ul>
   * <li> \c basis_matrix->get_fact_struc().get()!=NULL
   *      <ul>
   *      <li> \c this->get_fact_struc().get()==NULL
   *           <ul>
   *           <li> \c basis_matrix->get_fact_struc().count()==1
   *                <ul><li> No other basis matrix ojbect is using this factorization structure
   *                 so the storage in \c *basis_matrix->get_fact_struc() can be recycled and used here.</ul>
   *           <li> \c basis_matrix->get_fact_struc().count()>1
   *                <ul><li> Some other basis matrix object or some other client is using
   *                this factorization structure so new storage for the factorization
   *                structure must be allocated.</ul>
   *           </ul>
   *      <li> \c this->get_fact_struc().get()!=NULL
   *           <ul>
   *           <li> <tt>basis_matrix->%get_fact_struc().get() == this->%get_fact_struc().get()</tt>
   *                <ul>
   *                <li> \c basis_matrix->get_fact_struc().count()==2
   *                     <ul><li> These are the same factorization structure objects and
   *                     no other basis matrix object is using this factorization structure
   *                     so the storage in this object can be recycled and used here.</ul>
   *                <li> \c basis_matrix->get_fact_struc().count()>2
   *                     <ul><li> Some other basis matrix object or some other client is using
   *                     this factorization structure so new storage for the factorization
   *                     structure must be allocated.</ul>
   *                </ul>
   *           <li> <tt>basis_matrix->%get_fact_struc().get() != this->%get_fact_struc().get()</tt>
   *                <ul>
   *                <li> \c this->get_fact_struc().count()==1
   *                     <ul>
   *                     <li> \c basis_matrix->get_fact_struc().count()==1
   *                          <ul><li> One of the factorization structure objects can be
   *                          discarded and the other can be recycled here since no other
   *                          basis system object or other client has a reference to these
   *                          objects.</ul>
   *                     <li> \c basis_matrix->get_fact_struc().count()>1
   *                          <ul><li> The factorization structure associated with the basis
   *                          matrix object \c basis_matrix is being shared by some other
   *                          basis matrix object or other client but the current factorization
   *                          structure in \c this is not.  Therefore, we can recycle the storage
   *                          in \c *this->get_fact_struc() here.</ul>
   *                     </ul>
   *                <li> \c this->get_fact_struc().count()>1
   *                     <ul>
   *                     <li> \c basis_matrix->get_fact_struc().count()==1
   *                          <ul><li> The factorization structure associated with the basis
   *                          matrix object \c basis_matrix is not being shared by some other
   *                          basis matrix object or other client but the current factorization
   *                          structure in \c this is.  Therefore, we can recycle the storage
   *                          in \c *basis_matrix->get_fact_struc() here.</ul>
   *                     <li> \c basis_matrix->get_fact_struc().count()>1
   *                          <ul><li> The factorization structure objects associated with the basis
   *                          matrix object \c basis_matrix and \c this are being shared by some other
   *                          basis matrix objects or other clients.  Therefore, we will have to
   *                          allocate new factorization storage here.</ul>
   *                     </ul>
   *                </ul>
   *           </ul>
   *      <li> Current storage for the factorization nonzeros in \c *basis_matrix will be recycled here.
   *      </ul>
   * </ul>
   */
  virtual void analyze_and_factor(
    const AbstractLinAlgPack::MatrixConvertToSparse   &A
    ,DenseLinAlgPack::IVector                            *row_perm
    ,DenseLinAlgPack::IVector                            *col_perm
    ,size_type                                      *rank
    ,BasisMatrix                                    *basis_matrix
    ,std::ostream                                   *out            = NULL
    ) = 0;
  
  /** \brief Factor a matrix given its factorization structure.
   *
   * @param  A            [in] Matrix to be factored given a previously determined
   *                      factorization structure (see \c fact_struc).
   * @param  basis_matrix
   *                      [in/out] On input, this can be a previously initialized basis matrix
   *                      object in which case its storage for the factorization nonzeros can be
   *                      recycled if possible.
   * @param  fact_struc   [in] Smart pointer to the factorization structure to use.  If
   *                      \c fact_struc.get()==NULL then \c this->get_fact_struc() is used in its place.
   * @param  out          [in/out] If <tt>out != NULL</tt>, then some information about the
   *                      operations performed internally may be printed to \c *out.  The amount
   *                      of this output should be very minimal and should not significantly increase
   *                      with the size of the problem being solved.
   *
   * Preconditions:<ul>
   * <li> <tt>A</tt> has the same structure as used to form \c fact_struc (or \c this->get_fact_struc()
   *      if \c fact_struc.get()==NULL )
   * <li> <tt>this->get_fact_struc().get() != NULL || fact_struc.get() != NULL</tt>
   *      (throw <tt>NoCurrentBasisException</tt>)
   * <li> <tt>basis_matrix != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> [<tt>fact_struc.get() != NULL</tt>] <tt>basis_matrix->get_fact_struc().get() == fact_struc.get()</tt>
   *      and <tt>fact_struc.count()</tt> is incremented by one on output.
   * <li> [<tt>fact_struc.get() == NULL</tt>] <tt>basis_matrix->get_fact_struc().get() == this->get_fact_struc().get()</tt>
   *      and <tt>this->get_fact_struc().count()</tt> is incremented by one on output.
   * <li> <tt>*this->get_fact_struc()</tt> is unchanged in any way (however,
   *      <tt>this->get_fact_struc().cout()</tt> will be incremented by one
   *      if <tt>fact_struct.get() == NULL</tt>).
   * <li> If <tt>basis_matrix->get_fact_struct().get() != fact_struc.get()</tt>
   *      (or \c this->get_fact_struc().get() if \c fact_struc.get()==NULL) on
   *      input then if <tt>basis_matrix->get_fact_struct().count() == 1</tt>
   *      on input then this factorization structure will be discarded on output.
   * <li> If <tt>basis_matrix->get_fact_struc().get() == NULL</tt> on input, then
   *      new factorization <em>nonzeros</em> are allocated on output to hold the
   *      matrix factors.  If <tt>basis_matrix->get_fact_struc().get() != NULL</tt>
   *      on input, then the storage for nonzeros in \c *basis_matrix is recycled.
   * </ul>
   *
   * This method allows clients to factor a matrix given it predetermined
   * factorization structure.  The client can use the factorization
   * structure stored internally in \c this->fact_struc()
   * (\c fact_struc.get()==NULL) or use the factorization structure from another
   * precomputed basis matrix (\c fact_struc.get()!=NULL).  If the passed in matrix
   * \c A is obviously not compatible with the precomputed factorization structure
   * being used then an \c IncompatibleMatrixStructureException exception will be
   * thrown.  If \c A1 was the matrix originally passed to \c this->analyze_and_factor(A1,...)
   * and \c A2 is the matrix being passed to \c this->factor(A2,...) then if
   * \c A2.rows()!=A1.rows() or \c A2.cols()!=A1.cols() or \c A2.nz()!=A1.nz()
   * then the matrices are obviously not compatible and the exception will
   * be thrown.
   *
   * The argument \c basis_matrix can be used to recycle
   * the nonzero values of the factorization of an existing basis matrix.
   * For example, suppose that \c C1_ptr and \c C2_ptr point to two different
   * basis matrices with different structure (i.e. <tt>C1_ptr->get_fact_struc().get()
   * != C2_ptr->get_fact_struc().get()</tt>).  Now suppose we would like to factor
   * another matrix \c A3 that has the same structure of \c A2 by reusing the nonzero
   * values referenced in \c C1_ptr.  To do this the client would call:
   *
   \code
   direct_solver.factor( A3, C1_ptr.get(), C2_ptr->get_fact_struc() );
   \endcode
   *
   * Now consider the case where we would like to recycle the storage
   * for the nonzeros in \c C1_ptr and also use the same structure as C1.
   * To do this the client would call:
   *
   \code
   direct_solver.factor(A3,C1_ptr.get(),C1_ptr->get_fact_struc());
   \endcode
   *
   * Above, the client must explicitly pass C1_ptr->get_fact_struc() to use
   * this structure or the current interally stored structure given by
   * \c this->get_fact_struc() would be used instead.
   *
   * If the factorization can not be performed for some reason then the exception
   * \c FactorizationFailure will be thrown.
   *
   */
  virtual void factor(
    const AbstractLinAlgPack::MatrixConvertToSparse   &A
    ,BasisMatrix                                    *basis_matrix
    ,const BasisMatrix::fact_struc_ptr_t            &fact_struc    = Teuchos::null
    ,std::ostream                                   *out           = NULL
    ) = 0;

  /** \brief Get a smart pointer object to the current factorization structure.
   *
   * This returns a smart reference counted pointer to the factorization structure
   * computed from the last call to \c this->analyze_and_factor().  If no successful
   * call to \c this->analyze_and_factor() has been made yet then this function
   * will return <tt>return.get() == NULL</tt>.
   */
  virtual const BasisMatrix::fact_struc_ptr_t& get_fact_struc() const = 0;

  /** \brief Release all allocated memory.
   *
   * Postconditions:<ul>
   * <li> \c this->get_fact_struc().get()==NULL.
   * </ul>
   */
  virtual void set_uninitialized() = 0;

};	// end class DirectSparseSolver 

}	// end namespace AbstractLinAlgPack 

#endif	// ALAP_DIRECT_SPARSE_SOLVER_H
