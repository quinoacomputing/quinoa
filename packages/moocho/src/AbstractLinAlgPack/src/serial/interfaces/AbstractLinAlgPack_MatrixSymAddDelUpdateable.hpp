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

#ifndef MATRIX_SYM_ADD_DEL_UPDATEABLE_H
#define MATRIX_SYM_ADD_DEL_UPDATEABLE_H

#include "AbstractLinAlgPack_Types.hpp"
#include "AbstractLinAlgPack_MatrixBase.hpp"

namespace AbstractLinAlgPack {

/** \brief Mix-in Interface for updating a serial symmetric matrix by adding and deleting
 * rows and columns.
 *
 * This interface is designed to allow objects of this type to be used
 * in several different situations.  Generally, the matrix being updated
 * would be nonsingular but it does not have to be.  An important property
 * of symmetric matrix is its inertia.  The inertia of a matrix is the
 * number of negative, zero, and positive eigenvalues respectively.
 * While the eigenvalues themselves can be very difficult to compute,
 * in many cases the inertia is very easy to determine.  For instance,
 * if a <tt>A = L*D*L'</tt> factorization is being used, the inertia of the diagonal
 * matrix <tt>D</tt> is the same as for <tt>A</tt>.  Likewise, if a cholesky factorization
 * <tt>A = (+-)C*C'</tt> is used then it is easy to prove if the matrix is positive
 * definite (all positive eigen values) or negative definite (all negative
 * eigen values).  With other factorizations (such as QR for instance)
 * it is more difficult to determine the inertia and therefore it may
 * not be available.
 *
 * In inexact floating point arithmetic, it can be difficult to distingish
 * if a matrix is singular or nonsingular or if a matrix really has the
 * wrong inertia or is just singular.  In order to do this, tolerances
 * have to be identified.  In many contexts it is important for the client
 * to be able to specify these tolerances.  In order to establish a frame
 * of reference that is independent of the actual implementation of the
 * factorizations for the subclasses of this interface, we will use the
 * LU factorization:
 *
 * <tt>A = L*U</tt>
 *
 * where <tt>L</tt> is lower unit triangular and <tt>U</tt>
 * is upper nonunit triangular.  We will define the quantity:
 *
 * <tt>gamma = min{|U(i,i)|,i=1..n}/max{|U(i,i)|,i=1..n}</tt>
 *
 * as a measure of singularity.  Of course subclasses may not actually
 * use a LU factorization but by establishing this simple frame of
 * reference the tolerances can be properly interpreted by the subclasses.
 *
 * The classification of an factorized matrix will be as follows:
 \verbatim
  if (correct inertia) and (gamma > warning_tol) then
      The matrix is nonsingular and has the correct inertia, the
    initialization or update will succeed and all is good :-)
  elseif (correct inertia) and (singular_tol < gamma <= warning_tol) then
      The matrix will be considered nonsingular and the initialization
    or the update will succeed but a WarnNearSingularUpdateException
    will be thrown containing gamma and a warning message.
  elseif (correct inertia) and (0.0 < gamma <= singular_tol) then
      The matrix is considered singular, the initialization or update
    will not succeeed and a SingularUpdateException will be thrown
    containing gamma and an error message.
  elseif (gamma == 0.0) then
      The matrix is exactly singular, the initialization or update
    will not succeed and a SingularUpdateException will be thrown
    containing gamma and an error message.
  elseif (incorrect inertia) and (0.0 < gamma < wrong_inertia_tol) then
      The matrix will be considered singular, the initialization or update
    will not succeed and a SingularUpdateException will be thrown
    containing gamma and an error message.
  elseif (incorrect inertia) and (gamma >= wrong_inertia_tol) then
      The matrix is considered to be nonsingular but to have the wrong inertia,
    the initialization or update will not succeed and a WrongInertiaException
    will be thrown containing gamma and an error message.
  endif
 \endverbatim
 * The tolerances <tt>warning_tol</tt>, <tt>singular_tol</tt> and <tt>wrong_inertia_tol</tt> are
 * passed in as part of the struct <tt>PivotTolerances</tt> to each of the
 * methods that need them.  The default initialization for these tolerances
 * is <tt>PivotTolerances::UNKNOWN</tt> which means that the matrix object can do what
 * ever it wants to do.  It may use its own tolerances or none at all.  In
 * other words, the behavior is completely up to the subclass.
 * The idea of <tt>warning_tol</tt> and throwing an exception <tt>except</tt> of type 
 * <tt>WarnNearSingularUpdateException</tt> is to allow the updates to succeed but
 * return a warning message and the value of <tt>gamma</tt> as information to
 * the user.
 */
class MatrixSymAddDelUpdateable
  : public virtual AbstractLinAlgPack::MatrixBase // doxygen needs full name
{
public:

  /** @name Public types */
  //@{

  /** \brief . */
  enum EEigenValType { EIGEN_VAL_POS, EIGEN_VAL_NEG, EIGEN_VAL_ZERO, EIGEN_VAL_UNKNOWN };
  /** \brief Struct for the inertia of the matrix.
   *
   * Any or all of the values <tt>neg_eigens</tt>, <tt>zero_eigens</tt> or <tt>pos_eigens</tt>
   * may be <tt>UNKNOWN</tt>.
   */
  struct Inertia {
        enum { UNKNOWN = -1 };
    Inertia(
      int neg_eigen_vals    = UNKNOWN
      ,int zero_eigen_vals  = UNKNOWN
      ,int pos_eigen_vals  = UNKNOWN
      )
      : neg_eigens(neg_eigen_vals)
      ,zero_eigens(zero_eigen_vals)
      ,pos_eigens(pos_eigen_vals)
      {}
    /** \brief . */
    int  neg_eigens;
    /** \brief . */
    int  zero_eigens;
    /** \brief . */
    int  pos_eigens;
  }; 
  /** \brief Struct for pivot tolerances to be used when initializing, and augmenting
   * and deleting rows and columns.
   */
  struct PivotTolerances {
    enum { UNKNOWN = -1 };
    PivotTolerances()              // 2001/03/08: g++ 2.95.2 requries separate
      :warning_tol(UNKNOWN)      // constructor for use in default argument
      ,singular_tol(UNKNOWN)     // or you get internalcomplier error later?
      ,wrong_inertia_tol(UNKNOWN)
      {}
    PivotTolerances(
      value_type  _warning_tol
      ,value_type _singular_tol
      ,value_type _wrong_inertia_tol
      )
      :warning_tol(_warning_tol)
      ,singular_tol(_singular_tol)
      ,wrong_inertia_tol(_wrong_inertia_tol)
      {}
    /** \brief . */
    value_type warning_tol;
    /** \brief . */
    value_type singular_tol;
    /** \brief . */
    value_type wrong_inertia_tol;
  };
  /// Thrown if the matrix is near singular as a warning.
  class WarnNearSingularUpdateException : public std::logic_error	{
  public:
    WarnNearSingularUpdateException(const std::string& what_arg,value_type _gamma)
      : std::logic_error(what_arg), gamma(_gamma) {}
    value_type gamma;
  };
  /// Thrown if the matrix is singular and should not have been.
  class SingularUpdateException : public std::logic_error	{
  public:
    SingularUpdateException(const std::string& what_arg,value_type _gamma)
      : std::logic_error(what_arg), gamma(_gamma) {}
    value_type gamma;
  };
  /// Thrown if matrix has the wrong inertia from what was expected.
  class WrongInertiaUpdateException : public std::logic_error	{
  public:
    WrongInertiaUpdateException(const std::string& what_arg,value_type _gamma)
      : std::logic_error(what_arg), gamma(_gamma) {}
    value_type gamma;
  };
  /// Thrown if the maximum size is exceeded in augment_update(...).
  class MaxSizeExceededException : public std::logic_error
  {public: MaxSizeExceededException(const std::string& what_arg) : std::logic_error(what_arg) {}};

  //@}

  /** @name Public members to be overridden */
  //@{

  /** \brief . */
  virtual ~MatrixSymAddDelUpdateable()
  {}

  /** \brief Initialize to a 1x1 matrix.
   *
   * Since this is a 1x1 matrix the inetia is given by the sign
   * of alpha.
   *
   * @param  alpha    [in] The single entry in the 1x1 matrix to initialize.
   * @param  max_size [in] The maximum size for <tt>rows()</tt> and <tt>cols()</tt> the
   *                  maxtix is allowed to become.
   */
  virtual void initialize(
    value_type    alpha
    , size_type   max_size
    ) = 0;

  /** \brief Initialize given a symmetric matrix.
   *
   * The behavior of this function will vary based on the subclass that implements it.
   * Some subclasses may require that <tt>A</tt> be nonsingular and therefore <tt>inertia.zero_eigens</tt>
   * should be zero.  
   *
   * @param  A        [in] Symetric matrix that <tt>this</tt> is initialized with.
   * @param  max_size [in] The maximum size <tt>rows()</tt> and <tt>cols()</tt> can become.
   * @param  force_factorization
   *                  [in] If true, the factorization of the matrix will be forced and
   *                  any possible exceptions will be thrown.  If false then the factorization
   *                  may not be forced, in which case the client may not know immediatly
   *                  that the matrix is singular or has the wrong inertia.
   * @param  inertia  [in] The estimated inertia of the matrix.  If the user knows any
   *                  of the members of inertia then they should be set.  Some subclasses
   *                  may rely on this estimate of the inertia to determine what should be
   *                  done.
   * @param  pivot_tols
   *                  [in] Tolerances to use to determine singularity, nonsingularity etc.
   *                  See the intro.  Default is no tolerances.
   */
  virtual void initialize(
    const DMatrixSliceSym      &A
    ,size_type         max_size
    ,bool              force_factorization
    ,Inertia           inertia
    ,PivotTolerances   pivot_tols            = PivotTolerances()
    ) = 0;

  /** \brief Return the maximum size the matrix is allowed to become.
   */
  virtual size_type max_size() const = 0;

  /** \brief Return the inertia of the matrix (if it is known).
   * If any of the members of the inertia is not known then
   * they may be set to <tt>Inertia::UNKNOWN</tt>.  If the matrix is
   * nonsingular then <tt>return.zero_eigens == 0</tt> will be true.
   */
  virtual Inertia inertia() const = 0;

  /** \brief Set the matrix to uninitialized.
   */
  virtual void set_uninitialized() = 0;

  /** \brief Update by adding a symmetric row and column.
   *
   * The update performed is:
   \verbatim

   [ A     t   ]       
   [ t'  alpha ] ==>  A_new

   \endverbatim
   * Preconditions:<br>
   * \begin{itemize}
   * \item <tt>[t != NULL] t->size() == this->rows()</tt> (throw <tt>std::length_error</tt>)
   * \item <tt>this->rows() < this->max_size()</tt> (throw <tt>MaxSizeExceededException</tt>)
   * \end{itemize}
   *
   * Postcondiditons:<br>
   * The update gives a legal update depending on the
   * context of the subclass (nonsigular, positive definite etc.).
   * If the subclass requires the matrix to be nonsingular but 
   * <tt>inertia.zero_eigens == 0</tt> or the matrix is determined to be singular
   * then the exception <tt>SingularUpdateException</tt> will be thrown.
   * If the matrix is found to not have the propper inertia then the
   * exception <tt>WrongInertiaUpdateException</tt> will be thrown.  This subclass
   * may not be able to determine the inertia in which case this exception
   * will never be thrown.
   * If no exceptions are thrown then <tt>this->rows()</tt> and <tt>this->cols()</tt>
   * will increase by one and <tt>this->inertia()</tt> will return the new inertia
   * if it is known.
   *
   * @param  t       [in] DVectorSlice (size == <tt>rows()</tt>) where <tt>t</tt> may be <tt>NULL</tt> in which
   *                 case t is considered zero.
   * @param  alpha   [in] Scalar added.
   * @param  force_refactorization
   *                 [in] If true, then the factorization of the matrix will
   *                 be performed before the function returns.  If something
   *                 goes wrong then an exeception will be thrown here.
   * @param  add_eigen_val
   *                 [in] Gives the estimate of the new eigen value added
   *                 to the matrix.  If the matrix does not agree with this
   *                 then an exception will be thrown.
   * @param  pivot_tols
   *                  [in] Tolerances to use to determine singularity, nonsingularity etc.
   *                  See the intro.  Default is no tolerances.
   */
  virtual void augment_update(
    const DVectorSlice  *t
    ,value_type        alpha
    ,bool              force_refactorization = true
    ,EEigenValType     add_eigen_val         = EIGEN_VAL_UNKNOWN
    ,PivotTolerances   pivot_tols            = PivotTolerances()
    ) = 0;

  /** \brief Update by deleteing a symmetric row and column.
   *
   \verbatim
   
                 jd
       [ A11    a12    A13  ]
   A = [ a12'   a22    a23' ] jd  ==>  A_new = [ A11  A13  ]
       [ A13'   a23    A33  ]                  [ A13'  A33 ]
   
   \endverbatim
   *
   * Preconditions:<br>
   * \begin{itemize}
   * \item <tt>1 <= jd && jd <= this->rows()</tt> (throw <tt>std::out_of_range</tt>)
   * \end{itemize}
   *
   * Postcondiditons:<br>
   * The update give a legal update depending on the
   * context of the subclass (nonsigular, positive definite etc.).
   * Also <tt>rows()</tt> and <tt>cols()</tt> will decrease by one so this_after<tt>->rows()</tt> == this_before<tt>->rows()</tt> - 1.
   *
   * @param  jd      [in] The jth row and column to be removed from the matrix.
   * @param  force_refactorization
   *                 [in] If true, then the factorization of the matrix will
   *                 be performed before the function returns.  If something
   *                 goes wrong then an exeception will be thrown here.
   * @param  drop_eigen_val
   *                 [in] Gives the estimate of the eigen value dropped
   *                 from the matrix.  If the matrix does not agree with this
   *                 then an exception will be thrown.
   * @param  pivot_tols
   *                  [in] Tolerances to use to determine singularity, nonsingularity etc.
   *                  See the intro.  Default is no tolerances.
   */
  virtual void delete_update(
    size_type          jd
    ,bool              force_refactorization = true
    ,EEigenValType     drop_eigen_val        = EIGEN_VAL_UNKNOWN
    ,PivotTolerances   pivot_tols            = PivotTolerances()
    ) = 0;

  //@}

}; // end class MatrixSymAddDelUpdateable

}  // namespace AbstractLinAlgPack 

#endif // MATRIX_SYM_ADD_DEL_UPDATEABLE_H
