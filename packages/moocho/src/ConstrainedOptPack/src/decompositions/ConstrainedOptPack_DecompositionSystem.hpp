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

#ifndef DECOMPOSITION_SYSTEM_H
#define DECOMPOSITION_SYSTEM_H

#include <stdexcept>

#include "ConstrainedOptPack_Types.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"

namespace ConstrainedOptPack {

/** \brief This class abstracts a decomposition choice for the quasi-range
 * space \a Y and null space \a Z matrices for a linearly independent
 * set of columns of \a Gc.
 *
 * <tt>Gc = [ Gc(:,equ_decomp),  Gc(:,equ_undecomp) ]</tt>
 *
 * where \c Gc is <tt>n x m</tt>, \c Gc(:,equ_decomp) is <tt>n x r</tt> and
 * \c Gc(:,equ_undecomp) is <tt>n x (m - r)</tt>.
 *
 * Note that the columns in <tt>Gc(:,equ_undecomp)</tt> may be
 * linearly dependent with the columns in <tt>Gc(:,equ_undecomp)</tt>
 * or they may just be undecomposed linearly independent equality
 * constraints.
 *
 * The decomposition formed by subclasses must have the properties:
 \verbatim
   Z s.t. Gc(:,equ_decomp)' * Z = 0
   Y s.t. [Z  Y] is nonsingular
   R = Gc(:,equ_decomp)' * Y is nonsingular
   Uz = Gc(:,equ_undecomp)' * Z
   Uy = Gc(:,equ_undecomp)' * Y
 \endverbatim
 *
 * The matrix factory objects returned by ??? are ment to have a lifetime that is
 * independent of \c this.
 *
 * The decomposition matrices \c Z, \c Y, \c R, \c Uz and \c Uy which
 * are updated in <tt>this->update_decomp()</tt> must be completely independent from
 * \c this and from each other and \c Gc that they based on.  For example,
 * Once \c update_decomp() is called, \c this, \c Gc can be destroyed and
 * the behaviors of the decomposition matrices must not be altered.  In this respect
 * the <tt>%DecompositionSystem</tt> interface is really nothing more than a "Strategy"
 * interface (with some state data of course) for computing range/null decompositions.
 * This gives the client great flexibility in how the decomposition matrices are used. 
 *
 * ToDo: Finish documentation!
 */
class DecompositionSystem {
public:

  /** @name Public types */
  //@{

  /** \brief . */
  typedef Teuchos::RCP<
    const Teuchos::AbstractFactory<MatrixOpNonsing> >    mat_nonsing_fcty_ptr_t;
  /** \brief . */
  typedef Teuchos::RCP<
    const Teuchos::AbstractFactory<MatrixOp> >               mat_fcty_ptr_t;
  /** \brief . */
  class SingularDecomposition : public std::logic_error
  {public: SingularDecomposition(const std::string& what_arg) : std::logic_error(what_arg) {}};
  /** \brief . */
  class InvalidMatrixType : public std::logic_error
  {public: InvalidMatrixType(const std::string& what_arg) : std::logic_error(what_arg) {}};
  /** \brief . */
  class TestFailed : public std::runtime_error
  {public: TestFailed(const std::string& what_arg) : std::runtime_error(what_arg) {}};
  /// Enumeration for the amount of output to create from <tt>update_decomp()</tt>.
  enum EOutputLevel {
    PRINT_NONE          = 0,
    PRINT_BASIC_INFO    = 1,
    PRINT_MORE_INFO     = 2,
    PRINT_VECTORS       = 3,
    PRINT_EVERY_THING   = 4
    };
  /// Enumeration for if to run internal tests or not.
  enum ERunTests { RUN_TESTS, NO_TESTS };
  /** \brief . */
  enum EMatRelations { MATRICES_INDEP_IMPS, MATRICES_ALLOW_DEP_IMPS };

  //@}

  /** \brief . */
  virtual ~DecompositionSystem() {}

  /** @name Dimensionality of the decomposition */
  //@{

  /** \brief Return the number of rows in \c Gc.
   *
   * Postconditions:<ul>
   * <li> <tt>n > m</tt>
   * </ul>
   *
   * The default implementation returns
   * <tt>this->space_range()->dim() + this->space_null()->dim()</tt>.
   */
  virtual size_type n() const;

  /** \brief Return the number of columns in \c Gc.
   *
   * Postconditions:<ul>
   * <li> <tt>m > 0</tt>
   * </ul>
   */
  virtual size_type m() const = 0;

  /** \brief Returns the rank of \c Gc(:,equ_decomp()).
   *
   * Postconditions:<ul>
   * <li> <tt>r =< m</tt>
   * </ul>
   *
   * The default implementation returns
   * <tt>this->space_range()->dim()</tt>.
   */
  virtual size_type r() const;

  /** \brief Returns the range of the decomposed equalities.
   *
   * The default implementation returns <tt>Range1D(1,this->r())</tt>.
   */
  virtual Range1D equ_decomp() const;

  /** \brief Returns the range of the undecomposed equalities.
   *
   * The default implementation returns <tt>Range1D(this->r()+1,this->m())</tt>
   * or <tt>Range1D::Invalid</tt> if <tt>this->r() == this->m()<tt>
   */
  virtual Range1D equ_undecomp() const;

  //@}

  /** @name Range and null vector spaces */
  //@{

  /** \brief Return a \c VectorSpace object for the range space.
   *
   * Postconditions:<ul>
   * <li> <tt>return.get() != NULL</tt>
   * <li> <tt>return->dim() == this->r()</tt>
   * </ul>
   */
  virtual const VectorSpace::space_ptr_t space_range() const = 0;

  /** \brief Return a \c VectorSpace object for the range space.
   *
   * Postconditions:<ul>
   * <li> <tt>return.get() != NULL</tt>
   * <li> <tt>return->dim() == this->n() - this->r()</tt>
   * </ul>
   */
  virtual const VectorSpace::space_ptr_t space_null() const = 0;

  //@}

  /** @name Matrix factories */
  //@{

  /** \brief Return a matrix factory object for <tt>Z</tt>
   */
  virtual const mat_fcty_ptr_t factory_Z() const = 0;

  /** \brief Return a matrix factory object for <tt>Y</tt>
   */
  virtual const mat_fcty_ptr_t factory_Y() const = 0;

  /** \brief Return a matrix factory object for <tt>R</tt>.
   */
  virtual const mat_nonsing_fcty_ptr_t factory_R() const = 0;
  
  /** \brief Return a matrix factory object for <tt>Uz</tt>
   */
  virtual const mat_fcty_ptr_t factory_Uz() const = 0;

  /** \brief Return a matrix factory object for <tt>Uy</tt>
   */
  virtual const mat_fcty_ptr_t factory_Uy() const = 0;

  //@}

  /** @name Update range/null decomposition */
  //@{

  /** \brief Creates the range/null decomposition for <tt>Gc(:,equ_decomp)'</tt>.
   *
   * The decomposition is based on the linearly independent columns \c Gc(:,equ_decomp)
   * of \c Gc
   *
   * <tt>Gc = [ Gc(:,equ_decomp),  Gc(:,equ_undecomp) ]</tt>
   *
   * Specifically this operation finds the matrices:
   \verbatim
   Z s.t. Gc(:,equ_deomp)' * Z = 0
   Y s.t. [Z  Y] is nonsingular
   R = Gc(:,equ_decomp)' * Y is nonsingular
   Uz = Gc(:,equ_undecomp)' * Z
   Uy = Gc(:,equ_undecomp)' * Y
   \endverbatim
   * If there is some problem creating the decomposition then exceptions
   * with the base class \c std::exception may be thrown.  The meaning
   * of these exceptions are more associated with the subclasses
   * that implement this operation.
   *
   * The concrete types for <tt>Gc</tt>, <tt>Z</tt>, <tt>Y</tt>,
   * <tt>Uz</tt> and <tt>Uy</tt> must be compatable with
   * the concrete implementation of \c this or an <tt>InvalidMatrixType</tt> exeption
   * will be thrown..
   *
   * @param  out [out] If <tt>out!=NULL</tt> then output is printed to this stream
   *             depending on the value of \c olevel.
   * @param olevel
   *             [in] Determines the amount of output to print to \c *out.
   *             The exact type of output is determined by the implementing
   *             subclass but here is the sugguested behavior:<ul>
   *             <li> \c PRINT_NONE : Don't print anything (same as <tt>out==NULL</tt>).
   *             <li> \c PRINT_BASIC_INFO : Only print basic information about
   *                  how the decomposition is formed.
   *                  Amount of output = \c O(1).
   *             <li> \c PRINT_VECTORS : Prints out important vectors computed
   *                  durring the computations (usually only durring testing).
   *                  This level is only useful for debugging.
   *                  Amount of output = \c O(n).
   *             <li> \c PRINT_EVERY_THING : Print out nearly every important
   *                  quantity that is computed (except for the output matrices
   *                  themselves, clients can do that) while the matrices are being
   *                  formed or tests are being conducted.  This level is only
   *                  useful for debugging.
   *                  Amount of output = <tt>O(m*n)</tt>
   *             </ul>
   * @param test_what
   *             [in] Determines if internal validation tests are performed.
   *             The post conditions for the output matrices are not checked
   *             internally.  This is something that client can (and should)
   *             do independently (see \c DecompositionSystemTester).  Values:<ul>
   *             <li> \c RUN_TESTS : As many validation/consistency tests
   *                  are performed internally as possible.  If a test
   *                  fails then a \c TestFailed execption will be thrown.
   *                  The subclasses determine what the tests are and
   *                  what failing a test means.
   *             <li> \c NO_TEST : No tests are performed internally.  This is
   *                  to allow the fastest possible execution.
   *             </ul>
   *             If a test fails, then a \c TestFailed exception will be thrown with
   *             a helpful error message.
   * @param  Gc  [in] The matrix for which the range/null decomposition is defined.
   * @param  Z   [out] On output represents the <tt>n x (n-r)</tt> null space	matrix such that
   *             <tt>Gc(:,equ_decomp) * Z == 0</tt>.  This matrix object must have been created
   *             by <tt>this->factory_Z()->create()</tt>.
   * @param  Y   [out] On output represents the <tt>n x r</tt> range space matrix	such that
   *             <tt>[ Y  Z ]</tt> is nonsingular.  This matrix object must have been created
   *             by <tt>this->factory_Y()->create()</tt>.
   * @param  R   [out] On output represents the nonsingular <tt>r x r</tt> matrix <tt>Gc(:,equ_decomp) * Y</tt>.
   *             This matrix object must have been created by <tt>this->factory_R()->create()</tt>.
   * @param  Uz  [in/out] If <tt>Uz != NULL</tt> (<tt>this->m() > this->r()</tt> only) then on output
   *             <tt>*Uz</tt> represents the <tt>(m-r) x (n-r)</tt> matrix <tt>Gc(:,equ_undecomp) * Z</tt>.
   *             If <tt>this->m() == this->r()</tt> then <tt>Uz == NULL</tt> must be true.
   *             If <tt>Uz!=NULL</tt>, then this matrix object must have been created by
   *             <tt>this->factory_Uz()->create()</tt>.
   * @param  Uy  [in/out] If <tt>Uy != NULL</tt> (<tt>this->m() > this->r()</tt> only) then on output
   *             <tt>*Uy</tt> represents the <tt>(m-r) x r</tt> matrix <tt>Gc(:,equ_undecomp) * Y</tt>.
   *             If <tt>this->m() == this->r()</tt> then <tt>Uy == NULL</tt> must be true.
   *             If <tt>Uy!=NULL</tt>, then this matrix object must have been created by
   *             <tt>this->factory_Uy()->create()</tt>.
   * @param mat_rel
   *             [in] Determines if the ouput matrix objects must be completely independent or not.
   *             <ul>
   *             <li> MATRICES_INDEP_IMPS: The matrix objects must have independent implementations (default).
   *             <li> MATRICES_ALLOW_DEP_IMPS: The matrix objects can have implementation dependencies.
   *             </ul>
   *
   * Preconditions:<ul>
   * <li> <tt>Gc.rows() == this->n()</tt> (throw \c std::invalid_argument)
   * <li> <tt>Gc.cols() == this->m()</tt> (throw \c std::invalid_argument)
   * <li> [<tt>this->m() == this->r()</tt>] <tt>Uz == NULL</tt> (throw \c std::invalid_argument)
   * <li> [<tt>this->m() == this->r()</tt>] <tt>Uy == NULL</tt> (throw \c std::invalid_argument)
   * <li> <tt>Z!=NULL || Y!=NULL || R!=NULL || Uz!=NULL || Uy!=NULL</tt>
   *      (throw \c std::invalid_argument)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> [<tt>Z != NULL</tt>] <tt>Z.space_cols().is_compatible(Gc.space_cols()) == true)</tt>
   * <li> [<tt>Z != NULL</tt>] <tt>Z.space_rows().is_compatible(*space_null()) == true</tt>
   * <li> [<tt>Z != NULL</tt>] <tt>Gc(:,equ_decomp())' * Z == 0</tt>
   * <li> [<tt>Y != NULL</tt>] <tt>Y.space_cols().is_compatible(Gc.space_cols()) == true)</tt>
   * <li> [<tt>Y != NULL</tt>] <tt>Y.cols() == this->r()</tt>
   * <li> [<tt>Y != NULL</tt>] <tt>[ Y  Z ]</tt> is nonsingular
   * <li> [<tt>R != NULL</tt>] <tt>R->space_cols().is_compatible(*Gc.space_cols()->sub_space(equ_decomp())) == true</tt>
   * <li> [<tt>R != NULL</tt>] <tt>R->space_rows().is_compatible(*space_range()) == true</tt>
   * <li> [<tt>R != NULL</tt>] <tt>R == Gc(:,equ_decomp())'*Y</tt>
   * <li> [<tt>Uz != NULL</tt>] <tt>Uz.space_cols().is_compatible(*Gc.space_rows()->sub_space(equ_undecomp())) == true</tt>
   * <li> [<tt>Uz != NULL</tt>] <tt>Uz.space_rows().is_compatible(*space_null()) == true</tt>
   * <li> [<tt>Uz != NULL</tt>] <tt>Uz == Gc(:,equ_undecomp())'*Z</tt>
   * <li> [<tt>Uy != NULL</tt>] <tt>Uy.space_cols().is_compatible(*Gc.space_rows()->sub_space(equ_undecomp())) == true</tt>
   * <li> [<tt>Uy != NULL</tt>] <tt>Uy.space_rows().is_compatible(*space_range()) == true</tt>
   * <li> [<tt>Uy != NULL</tt>] <tt>Uy == Gc(:,equ_undecomp())'*Y</tt>
   * <li> [<tt>mat_rel == MATRICES_INDEP_IMPS</tt>] The behaviors of all of the participating output matrix
   *      objects must not be altered by changes to other matrix objects.
   * <li> [<tt>mat_rel == MATRICES_ALLOW_DEP_IMPS</tt>] The behaviors of all of the participating output matrix
   *      objects may change when other matrix objects or \c this is altered.
   * </ul>
   *
   * Note that this method requires that all of the output matrix objects \c Z, \c Y, \c R,
   * \c Uz, \c Uy, \c Vz and \c Vy must be independent of \c this and of each other.
   * For example, the behavior of \c R must be be altered if \c this is destroyed or if
   * \c Z is modified of destroyed.  This requirment constrains the implementations somewhat
   * but makes things much easier for the client and gives the client much more power.
   */
  virtual void update_decomp(
    std::ostream          *out
    ,EOutputLevel         olevel
    ,ERunTests            test_what
    ,const MatrixOp       &Gc
    ,MatrixOp             *Z
    ,MatrixOp             *Y
    ,MatrixOpNonsing      *R
    ,MatrixOp             *Uz
    ,MatrixOp             *Vy
    ,EMatRelations        mat_rel = MATRICES_INDEP_IMPS
    ) const = 0;
  
  /** \brief Print the sub-algorithm by which the decomposition is formed
   *
   */
  virtual void print_update_decomp(
    std::ostream& out, const std::string& leading_str ) const = 0;

  //@}
  
};	// end class DecompositionSystem

}	// end namespace ConstrainedOptPack

#endif // DECOMPOSITION_SYSTEM_H
