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

#ifndef ABSTRACT_LIN_ALG_PACK_BASIS_SYSTEM_H
#define ABSTRACT_LIN_ALG_PACK_BASIS_SYSTEM_H

#include "AbstractLinAlgPack_Types.hpp"
#include "Teuchos_AbstractFactory.hpp"
#include "Teuchos_RCP.hpp"

namespace AbstractLinAlgPack {

/** \brief Interface for the creation and maintainance of a basis matrix for a decomposition of
 * linearlized constriants.
 *
 * <b>Overview:</b>
 *
 * This interface is designed to take the Jacobian for a sub-set of equality constraints
 * \f$ \nabla c \f$ and to create a basis matrix.
 * Assume we have the folloing linealrized equality constraints:
 *
 \f[ \nabla c^T d + c = 0 \f]
 *
 * The C++ identifier given to \f$ \nabla c \f$ is <tt>Gc</tt>.
 *
 * In this basis interface we will assume that <tt>d</tt>, <tt>c</tt> and <tt>h</tt> are
 * sorted such that we define the following sets (given the partitioning matrices
 * \f$ Q_{x} = \left[\begin{array}{c} Q_{xD} \\ Q_{xI} \end{array}\right] \f$,
 * \f$ Q_{c} = \left[\begin{array}{c} Q_{cD} \\ Q_{cU} \end{array}\right] \f$,
 * ):
 * <ul>
 * <li> <tt>d(var_dep)</tt> (\f$ d_D = Q_{xD} d\f$) : Dependent (i.e. basis) variables.
 * <li> <tt>d(var_indep)</tt> (\f$ d_I = Q_{xI} d\f$) : Independent (i.e. nonbasic) variables.
 * <li> <tt>c(equ_decomp)</tt> (\f$ c_D = Q_{cD} c\f$) : Decomposed equality constriants.
 * <li> <tt>c(equ_undecomp)</tt> (\f$ c_U = Q_{cU} c\f$): Undecomposed equality constriants.
 * </ul>
 * Given these partitionings we can define a basis matrix \a C for the
 * following Jacobian sub-matrices (in mathematical and Matlab-like notation):
 \f[
    C = Q_{cD} \nabla c_T (Q_{xD})^T
 \f]
 \verbatim

 C = Gc(var_dep,equ_decomp)'
 \endverbatim
 * We can also define a nonbasis matrix \a N for the decomposed constraints as:
 \f[
    N = Q_{cD} \nabla c^T (Q_{xI})^T
 \f]
 \verbatim

 N = Gc(var_indep,equ_decomp)'
 \endverbatim
 * Given the definitions of \a C and \a N above, we can define the following
 * direct-sensitivity matrix <tt>D</tt>:
 \f[
   D = - C^{-1} N
 \f]
 \verbatim

  D = -inv(C)*N
 \endverbatim
 * Given this matrix \a D, we can define another projected sensistivity matrix:
 * <ul>
 * <li> <tt>GcUP = Gc(var_indep,equ_undecomp)'   + Gc(var_dep,equ_undecomp)'   * D</tt>
 * </ul>
 *
 * This interface allows a client to create the basis matrix <tt>C</tt> and optionally
 * the direct sensitivity matrix <tt>D = -inv(C)*N</tt> and the auxiliary projected
 * sensistivity matrix <tt>GcUP</tt> (shown above).  These matrix
 * objects are independent from \c this \c BasisSystem object or from other \a C, \a D,
 * or \c GcUP objects.  Therefore, a <tt>%BasisSystem</tt> object can be thought of
 * as an "Abstract Factory" for basis matrices and auxillary matrices.  Note that
 * a <tt>%BasisSystem</tt> object will not compute the matrices \c D, \c GcUP
 * unless specifically asked.
 *
 * Note that the purpose of this interface is to abstract client code away from the
 * details of how the matrix \c Gc is represented and implemented and how
 * the basis matrix \a C is formed and implemented.  The complexity of these matrices could
 * vary from simple dense serial matrices all the way up massively parallel matrices using
 * iterative solvers for \a C running on computers with thousands of nodes.
 *
 * This interface also allows clients to compute the matrices <tt>J = D'*D</tt> and
 * <tt>S = I + D'*D</tt>.
 *
 * <b>Client usage:</b>
 *
 * The matrix objects for <tt>C</tt>, <tt>D</tt>, <tt>GcUP</tt>,
 * <tt>D'*D</tt> and <tt>S=I+D'*D</tt> are created by the client
 * using the \c AbstractFactory<> objects returned from
 * <tt>factory_C()</tt>, <tt>factory_D()</tt>,
 * <tt>factory_GcUP()</tt>, <tt>factory_transDtD()</tt> and
 * <tt>factory_S()</tt> respectively.  These methods return smart
 * pointers to these matrix factory objects and these objects are ment
 * to have a lifetime that extends up to and beyond the lifetime of
 * the <tt>%BasisSystem</tt> object that created them.  Note that the
 * matrix objects returned by these matrix factory objects are not to
 * be considered usable until they have passed through
 * <tt>update_basis()</tt> or receive some other appropriate
 * initialization.
 *
 * The ranges of the dependent and independent variables, and decomposed and undecomposed
 * equality constriants are returned by the methods \c var_dep(), \c var_indep(),
 * \c equ_decomp() and\c equ_undecomp() respectively.  There are a few obvious assertions
 * for the values that these ranges can take on.  Assuming that \c Gc is non-null
 * when passed to \c update_basis(), the following assertions apply:
 *
 * <A NAME="ranges_assertions"></A>
 * Assertions:<ul>
 * <li> <tt>var_dep().size() == equ_decomp().size() + inequ_decomp().size()</tt>
 * <li> <tt>var_dep().size() + var_indep().size() == Gc.rows()</tt>
 * <li> <tt>equ_decomp().size() + equ_undecomp().size() == Gc.cols()</tt>
 * </ul>
 *
 * Note that the client should not rely on \c var_dep(), \c
 * var_indep(), \c equ_decomp(), or \c equ_undecomp() until after the
 * first call to \c update_basis().  This allows a
 * <tt>%BasisSystem</tt> object to adjust itself to accommodate the
 * input matrix \c Gc.
 *
 * A fully initialized <tt>%BasisSystem</tt> object will be setup to
 * work with specific types and sizes of input matrices \c Gc and \c
 * Gh.  Therefore, the client should be able to get accrate values
 * from \c var_dep(), \c var_indep(), \c equ_decomp(), or \c
 * equ_undecomp() even before the first call to \c update_basis().
 * The <tt>%BasisSystem</tt> object must therefore be initialized in
 * some way to accommodate input matrices \c Gc and \c Gh of a
 * specific dimension.
 *
 * Note that This interface is completely worthless unless
 * \c var_dep() returns some valid range (i.e. a basis matrix exists).
 * If <tt>var_dep().size() == 0</tt> then this is an indication that
 * \c this is uninitialzed and therefore only the factory methods can
 * be called!
 *
 * The method \c update_basis() is used by the client to update the
 * basis matrix \a C and perhaps the direct sensitivity matrix \a D
 * and it's auxillary projected sensistivity matrix \c GcUP
 * Strictly speaking, it would be possible to form the matrix
 * \a D externally through the <tt>MatrixNonsing</tt> interface
 * using the returned \a C and an \a N matrix object, but this may not
 * take advantage of any special application specific tricks that can
 * be used to form \a D.  Note that this interface does not return a
 * nonbasis matrix object for \a N.  However, this matrix object will
 * be needed for an implicit \a D matrix object that the client will
 * undoubtably want to create.  Creating such a matrix object is
 * simple given the method <tt>MatrixOp::sub_view()</tt>.  The
 * following code example shows how to create a matrix object for \a N
 * (given the matrix \c Gc input to <tt>bs.update_basis(Gc,...)</tt> and \c bs):
 \code
Teuchos::RCP<const MatrixOp>
create_N(
    const AbstractLinAlgPack::MatrixOp       &Gc
    ,const AbstractLinAlgPack::BasisSystem   &bs
    )
{
    namespace mmp = MemMngPack;
  return Teuchos::rcp(
        new MatrixOpSubView(
            Gc.sub_view(bs.var_indep(),bs.equ_decomp()), BLAS_Cpp::trans
            )
        );
}
 \endcode
 * Given the nonbasis matrix object for \a N returned by the above function, this matrix object could be used
 * to form an explicit \a D matrix object (but perhaps not very efficiently) or be used to implicitly implement
 * matrix vector products with \a D as:
 \verbatim

 op(D)*v = op(-inv(C)*N)*v = -inv(C)*(N*v) or -N'*(inv(C')*v)
 \endverbatim
 *
 * The client can also form matrices of the form <tt>S = I + D'*D</tt> as follows:
 \code
 Teuchos::RCP<MatrixSymOpNonsing>
     S = basis_sys.factory_S()->create();
 Teuchos::dyn_cast<MatrixSymInitDiag>(*S).init_identity(D.space_rows());
 syrk(D,BLAS_Cpp::trans,1.0,1.0,S.get();
 \endcode
 * The matrix <tt>S</tt> must then be fully initialized and ready to go.
 *
 * <b>Subclass developer's notes:</b>
 *
 * The default implementation (of the methods that have default implementations) assume
 * that there are no undecomposed equality constriants.
 * 
 * ToDo: Finish documentation!
 *
 */
class BasisSystem {
public:

  /** @name Public types */
  //@{

  /** \brief . */
  typedef Teuchos::RCP<
    const Teuchos::AbstractFactory<MatrixOpNonsing> >    mat_nonsing_fcty_ptr_t;
  /** \brief . */
  typedef Teuchos::RCP<
    const Teuchos::AbstractFactory<MatrixOp> >           mat_fcty_ptr_t;
  /** \brief . */
  typedef Teuchos::RCP<
    const Teuchos::AbstractFactory<MatrixSymOp> >        mat_sym_fcty_ptr_t;
  /** \brief . */
  typedef Teuchos::RCP<
    const Teuchos::AbstractFactory<MatrixSymOpNonsing> > mat_sym_nonsing_fcty_ptr_t;
  /** \brief . */
  class SingularBasis : public std::runtime_error
  {public: SingularBasis(const std::string& what_arg) : std::runtime_error(what_arg) {}};
  /** \brief . */
  enum EMatRelations { MATRICES_INDEP_IMPS, MATRICES_ALLOW_DEP_IMPS };

  //@}

  /** \brief Required constructor (calls <tt>initialize()</tt>).
   */
  BasisSystem(
    const mat_sym_fcty_ptr_t             &factory_transDtD
    ,const mat_sym_nonsing_fcty_ptr_t    &factory_S
    );

  /** \brief Initialize the factory objects for the special matrices for <tt>D'*D</tt> and <tt>S = I + D'*D</tt>.
   *
   * Postconditions:<ul>
   * <li>this->factory_transDtD().get() == factory_transDtD.get()</tt>
   * <li>this->factory_S().get() == factory_S.get()</tt>
   * </ul>
   */
  virtual void initialize(
    const mat_sym_fcty_ptr_t             &factory_transDtD
    ,const mat_sym_nonsing_fcty_ptr_t    &factory_S
    );

  /** \brief . */
  virtual ~BasisSystem() {}

  /** @name Matrix factories */
  //@{

  /** \brief Return a matrix factory object for basis <tt>C = [ Gc(var_dep,equ_decomp)';  Gh(var_dep,inequ_decomp)' ]</tt>.
   */
  virtual const mat_nonsing_fcty_ptr_t factory_C() const = 0;
  
  /** \brief Return a matrix factory object for sensitivity matrix <tt>D = -inv(C)*N</tt>.
   *
   * It is allowed for this to return \c NULL in which case \c update_basis() will not
   * accept a \c D matrix to be computed.
   */
  virtual const mat_fcty_ptr_t factory_D() const = 0;

  /** \brief Return a matrix factory object for auxiliary sensitivity matrix <tt>GcUP = Gc(var_indep,equ_undecomp)' + Gc(var_dep,equ_undecomp)'*D</tt>.
   *
   * It is allowed for this to return \c NULL in which case \c update_basis() will not
   * accept a \c GcUP matrix to be computed.
   */
  virtual const mat_fcty_ptr_t factory_GcUP() const;

  /** \brief Returns a matrix factory for the result of <tt>J = D'*D</tt>
   * 
   * The resulting matrix is symmetric but is assumed to be singular.
   *
   * Postconditions:<ul>
   * <li> The function <tt>AbstractLinAlgPack::syrk(D,trans,alpha,beta,return->create().get())</tt>
   *      must not throw an exception once <tt>D</tt> has been initialized by <tt>this</tt>.
   * </ul>
   */
  virtual const mat_sym_fcty_ptr_t factory_transDtD() const;
  
  /** \brief Returns a matrix factory for the result of <tt>S = I + D'*D</tt>
   * 
   * The resulting matrix is symmetric and is guarrenteed to be nonsingular.
   *
   * Postconditions:<ul>
   * <li><tt>dynamic_cast<MatrixSymInitDiag*>(return->create().get()) != NULL</tt>
   * </ul>
   */
  virtual const mat_sym_nonsing_fcty_ptr_t factory_S() const;
  
  //@}

  /** @name Return the ranges for variable and constraint partitioning */
  //@{

  /** \brief Range of dependent (basic) variables.
   *
   * If there are no dependent variables then <tt>return.size() == 0</tt>.
   * This would be a strange case where there was no basis matrix in which
   * case this whole interface would be worthless.  Therefore, to be useful
   * <tt>return.size() > 0</tt> must be true.
   *
   * If \c var_dep().size() returns 0, then this is an indication that
   * \c *this is uninitialized and only the factory methods can be
   * called.
   */
  virtual Range1D var_dep() const = 0;
  /** \brief Range of independnet (nonbasic) variables.
   *
   * It is possible that the basis matrix may take up all of the degrees of
   * freedom with <tt>var_dep().size() == Gc->rows()</tt>.  In this case, there
   * is no nonbasis matrix \a N and no direct sensitivity matrix \a D.
   * In this case <tt>return.size() == 0</tt>.  In the more general case
   * however, <tt>return.size() > 0</tt>.
   */
  virtual Range1D var_indep() const = 0;
  /** \brief Range of decomposed general equality constraints.
   *
   * If there are no decomposed general equality constriants then
   * <tt>return.size() == 0</tt>.  Otherwise, <tt>return.size() > 0</tt>.
   *
   * The default implementation return <tt>Range1D(1,this->var_dep().size())</tt>
   */
  virtual Range1D equ_decomp() const;
  /** \brief Range of undecomposed general equality constriants.
   *
   * If there are no undecomposed equality constriants then
   * <tt>return.size() == 0</tt>.  Otherwise, <tt>return.size() > 0</tt>.
   *
   * The default implementation return <tt>Range1D::Invalid</tt>
   */
  virtual Range1D equ_undecomp() const;

  //@}

  /** @name Update matrices */
  //@{

  /** \brief Update a basis and posssibly the direct sensitivity matrix for a 
   * set of Jacobian matrices.
   *
   * @param  Gc    [in] Jacobian of the equality constriants.
   * @param  C     [out] Basis matrix.  If <tt>C == NULL</tt> on input, then this
   *               quantity is not updated.  If <tt>C != NULL</tt> then this must
   *               have been created by <tt>this->factory_C()->create()</tt>.
   *               This basis matrix object must be independent of the input
   *               matrices \c Gc and/or \c Gh.  Therefore, it must be legal to
   *               destroy \c Gc and/or \c Gh without affecting the behavior of
   *               the basis matrix object \c C.
   * @param  D     [out] Direct sensitivity matrix <tt>D = -inv(C)*N</tt>.  If
   *               <tt>D == NULL</tt> on input then this quantity is not updated.
   *               If <tt>D != NULL</tt> then this must have been created by
   *               <tt>this->factory_D()->create()</tt>.  This matrix object
   *               is meaningless if <tt>this->var_indep() == Range1D::Invalid</tt>
   *               on return.
   *               This matrix object must be independent matrices \c Gc and/or \c Gh
   *               Therefore, it must be legal to
   *               destroy \c Gc and/or \c Gh without affecting the behavior of
   *               the direct sensitivity matrix object \c D.
   * @param  GcUP  [out] Auxiliary sensistivity matrix
   *               <tt>GcUP = Gc(var_indep,equ_undecomp)' + Gc(var_dep,equ_undecomp)'*D</tt>.
   *               If <tt>GcUP == NULL</tt> on input then this quantity is not updated.
   *               If <tt>GcUP != NULL</tt> then this must have been created by
   *               <tt>this->factory_GcUP()->create()</tt>.  This matrix object
   *               is meaningless if <tt>this->var_indep() == Range1D::Invalid</tt>
   *               on return.
   *               This matrix object must be independent of the matrices \c Gc and/or \c Gh
   *               and/or \c D.  Therefore, it must be legal to destroy \c Gc and/or \c Gh
   *               and/or \c D without affecting the behavior of the matrix object \c GcUP.
   * @param mat_rel
   *               [in] Determines if the matrix objects must be completely independent or not.
   *               <ul>
   *               <li> MATRICES_INDEP_IMPS: The matrix objects must have independent implementations (default).
   *               <li> MATRICES_ALLOW_DEP_IMPS: The matrix objects can have implementation dependencies.
   *               </ul>
   * @param  out   [in/out] If <tt>out!=NULL</tt>, then some information about the operations performed
   *               internally may be printed to \c *out.  The amount of this output should be
   *               very minimal and should not significantly increase with the size of the problem
   *               being solved.
   *
   * Preconditions:<ul>
   * <li> <tt>Gc != NULL || Gh != NULL</tt>
   * <li> [<tt>Gc != NULL && Gh != NULL</tt>]
   *      <tt>Gc->space_cols().is_compatible(Gh->space_cols()) == true</tt>
   * <li> [<tt>Gc != NULL</tt>] <tt>Gc->space_cols().sub_space(var_dep()).get() != NULL</tt>
   * <li> [<tt>Gc != NULL</tt>] <tt>Gc->space_cols().sub_space(var_indep()).get() != NULL</tt>
   * <li> [<tt>Gc != NULL</tt>] <tt>Gc->space_rows().sub_space(equ_decomp()).get() != NULL</tt>
   * <li> [<tt>Gc != NULL && equ_decomp().size() > 0 </tt>]
   *      <tt>Gc.space_rows().sub_space(equ_decomp()).get() != NULL</tt>
   * <li> [<tt>Gc != NULL && equ_undecomp().size() > 0 </tt>]
   *      <tt>Gc.space_rows().sub_space(equ_undecomp()).get() != NULL</tt>
   * <li> <tt>C != NULL || D != NULL || GcUP != NULL</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->var_dep() != Range1D::Invalid && !this->var_dep().full_range()</tt>
   * <li> [<tt>C != NULL && Gc != NULL && equ_decomp().size() > 0</tt>]
   *      <tt>C->space_cols().sub_space(equ_decomp())->is_compatible(Gc->space_rows().sub_space(equ_decomp()))
   *      && C->space_rows().is_compatible(Gc->space_cols().sub_space(var_dep()))</tt>
   * <li> [<tt>C != NULL && Gh != NULL && inequ_decomp().size() > 0</tt>]
   *      <tt>C->space_cols().sub_space(equ_decomp().size()+inequ_decomp())->is_compatible(Gh->space_rows().sub_space(inequ_decomp()))
   *      && C->space_rows().is_compatible(Gh->space_cols().sub_space(var_dep()))</tt>
   * <li> [<tt>D != NULL && Gc != NULL && var_indep().size() > 0 && equ_decomp().size() > 0</tt>]
   *      <tt>D->space_cols().sub_space(equ_decomp())->is_compatible(Gc->space_rows().sub_space(equ_decomp()))
   *      && D->space_rows().is_compatible(Gc->space_cols().sub_space(var_indep()))</tt>
   * <li> [<tt>D != NULL && Gh != NULL && var_indep().size() > 0 && inequ_decomp().size() > 0</tt>]
   *      <tt>D->space_cols().sub_space(equ_decomp().size()+inequ_decomp())->is_compatible(Gh->space_rows().sub_space(inequ_decomp()))
   *      && D->space_rows().is_compatible(Gh->space_cols().sub_space(var_indep()))</tt>
   * <li> [<tt>GcUP != NULL && var_indep().size() > 0 && equ_undecomp().size() > 0</tt>]
   *      <tt>GcUP->space_rows()->is_compatible(Gc->space_cols().sub_space(var_indep()))
   *      && GcUP->space_cols()->is_compatible(Gc->space_rows().sub_space(equ_undecomp()))</tt>
   * </ul>
   *
   * This method with throw a \c SingularBasis exception if the updated basis matrix \a C is too close
   * (as defined by the underlying implementation by some means) to being numerically singular.
   */
  virtual void update_basis(
    const MatrixOp          &Gc
    ,MatrixOpNonsing        *C
    ,MatrixOp               *D
    ,MatrixOp               *GcUP
    ,EMatRelations          mat_rel = MATRICES_INDEP_IMPS
    ,std::ostream           *out    = NULL
    ) const = 0;

  //@}

private:
  mat_sym_fcty_ptr_t             factory_transDtD_;
  mat_sym_nonsing_fcty_ptr_t     factory_S_;

  // not defined and not to be called
  BasisSystem();

}; // end class BasisSystem

}  // end namespace AbstractLinAlgPack

#endif // ABSTRACT_LIN_ALG_PACK_BASIS_SYSTEM_H
