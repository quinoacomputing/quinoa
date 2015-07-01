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

#ifndef MATRIX_KKT_FULL_SPACE_RELAXED_H
#define MATRIX_KKT_FULL_SPACE_RELAXED_H

#include "ConstrainedOptPack_Types.hpp"
#include "AbstractLinAlgPack/src/MatrixWithOpFactorized.hpp"
#include "AbstractLinAlgPack/src/MatrixConvertToSparseFortranCompatible.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace ConstrainedOptPack {

/** \brief Implementation of a KKT matrix factorized in the full space.
  *
  * This class is used to represent the KKT matrix of the following
  * relaxed QP:
  *
  \begin{verbatim}
  min     [ g'  M ] * [  d  ] + 1/2 * [ d'  eta ] * [ G      ] * [  d  ]
                      [ eta ]                       [     M  ]   [ eta ]

  s.t.    [ A'  -c ] * [  d  ] + c = 0
                       [ eta ]
  \end{verbatim}
  *
  * The only matrix actually factorized is:
  *
  \begin{verbatim}
  K_bar = [ G  A ]
          [ A'   ]
  \end{verbatim}
  *
  * The class has two modes.
  *
  * First mode is to not include the relaxation
  * term and therefore the KKT matrix is:
  *
  \begin{verbatim}
  K = [ G  A ]
      [ A'   ]
  \end{verbatim}
  *
  * The second mode is the use the relaxation and he represented matrix is:
  *
  \begin{verbatim}
      [ G       A   ]
  K = [     M   -c' ]
      [ A'  -c      ]
  \end{verbatim}
  *
  * This class uses an aggregate DirectSparseFortranCompatibleSolver (DSFCS) 
  * object to factorize K above and then to solve for the linear systems
  * involving K.
  */
class MatrixKKTFullSpaceRelaxed
  : public MatrixWithOpFactorized
  , public MatrixConvertToSparseFortranCompatible
{
public:
  
  /** \brief . */
  typedef AbstractLinAlgPack::DirectSparseFortranCompatibleSolver
    DirectSparseFortranCompatibleSolver;

  /** \brief . */
  class NotInitializedException : public std::logic_error
  {public: NotInitializedException (const std::string& what_arg) : std::logic_error(what_arg) {}};

  /** \brief . */
  class SingularMatrixException : public std::logic_error
  {public: SingularMatrixException (const std::string& what_arg) : std::logic_error(what_arg) {}};

  /** \brief . */
  class InvalidMatrixType : public std::logic_error
  {public: InvalidMatrixType (const std::string& what_arg) : std::logic_error(what_arg) {}};

  /** \brief . */
  enum ERunTests { RUN_TESTS, NO_TESTS };

  /** \brief . */
  enum EPrintMoreOrLess { PRINT_MORE, PRINT_LESS };

  /// <<std comp>> members for the direct sparse linear solver
  STANDARD_COMPOSITION_MEMBERS( DirectSparseFortranCompatibleSolver, direct_solver );

  /** \brief . */
  MatrixKKTFullSpaceRelaxed( const direct_solver_ptr_t& direct_solver = 0 );

  /** @name Initialize the relaxed or unrelaxed KKT matrix.
    *
    * These operations will factorize the matrix K.  If the matrix K
    * is not full rank then a SingularMatrixException exception will
    * be thrown.  The objects G and A must support the
    * MatrixConvertToSparseFortranCompatible (MCTSFC) interface or
    * the exception InvalidMatrixType will be thrown.
    *
    * Some of the common arguments that these initialization methods share
    * are:
    *
    * \begin{itemize}
    *	\item G [I] Hessian matrix ( must support MESFCE interface ).
    *	\item A [I] Gradient of constraints matrix ( must support MESFCE interface ).
    *	\item out [O] Output stream to print to.  This stream may be used for
    *		output after initialization also so make sure that it remains valid
    *		as long as this matrix object is is use.  For no output set out=NULL.
    *	\item run_test [I] If set the true then many (expensive) tests will be
    *		preformed to ensure that everything is working properly. 
    *	\item print_more [I] If set the true then a lot more output may be produced
    *		expecially if some error occurs. 
    * \end{itemize}
    *
    * Important: It is vital that the definitions of G and A do not change
    * externally while this object is being used.  To do so may invalidate
    * the behavior of this object (especially the MatrixOp functions).
    *
    * This class will try to reuse the factorization structure from the
    * last call to initialze(...) or initialize_relaxed(...) when possible.
    * Namely if G and A have the same dimensions and same number of nonzeros
    * of the matrices previously factorized, it will be assumed that the
    * structure will be the same.  If this is not the case then the
    * client should call release_memory(...) to wipe the slate clean and
    * start over before calling initialize...(...) again.
    */
  //@{

  /** \brief Initialize the nonrelaxed matrix.
    *
    */
  void initialize( const MatrixOp& G, const MatrixOp& A
    , std::ostream* out = 0, EPrintMoreOrLess print_what = PRINT_LESS
    , ERunTests test_what = NO_TESTS );

  /** \brief Initialize the relaxed matrix.
    *
    * If the unrelaxed QP is well scaled (near 1.0) then a reasonable
    * value for bigM = M might be 1e+10 however this is problem specific.
    */
  void initialize_relaxed( const MatrixOp& G, const MatrixOp& A
    , const DVectorSlice& c, value_type bigM = 1e+10
    , std::ostream* out = 0, EPrintMoreOrLess print_what = PRINT_LESS
    , ERunTests test_what = NO_TESTS );

  /** \brief Set the matrix to uninitialized.
    *
    * The purpose of this method is for the client to specifically state that
    * it is done using this object for now.  This is to avoid problems where
    * the definitions of G and A might change and then another client unknowingly
    * trys to use this object.
    *
    * Note that this does not erase storage of the factorization structure
    * for example.
    */
  void set_uninitialized();

  /** \brief Clear all allocated storage.
    *
    * The client should call this routine if he wants the new KKT matrix
    * to be reanalyze and factorized the next time initialize...(...) is
    * called.
    */
  void release_memory();

  //@}

  // /////////////////////////////////////////////////////
  // Overridden from Matrix

  /** \brief . */
  size_type rows() const;

  /** \brief . */
  size_type cols() const;

  // /////////////////////////////////////////////////////////
  // Overridden from MatrixOp

  /** \brief . */
  std::ostream& output(std::ostream& out) const;

  /** \brief . */
  MatrixOp& operator=(const MatrixOp& m);

  /// (2) vs_lhs = alpha * op(M_rhs1) * vs_rhs2 + beta * vs_lhs (BLAS xGEMV)
  void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const DVectorSlice& vs_rhs2, value_type beta) const;

  // ////////////////////////////////////////////////////////////
  // Overridden from MatrixFactorized

  /// (1) v_lhs	= inv(op(M_rhs1)) * vs_rhs2
  void V_InvMtV( DVectorSlice* v_lhs, BLAS_Cpp::Transp trans_rhs1
    , const DVectorSlice& vs_rhs2) const;

  // ////////////////////////////////////////////////////////////
  // Overridden from MatrixConvertToSparseFortranCompatible

  /** \brief . */
  FortranTypes::f_int num_nonzeros( EExtractRegion extract_region ) const;

  /** \brief . */
  void coor_extract_nonzeros(
      EExtractRegion extract_region
    , const FortranTypes::f_int len_Aval
      , FortranTypes::f_dbl_prec Aval[]
    , const FortranTypes::f_int len_Aij
      , FortranTypes::f_int Arow[]
      , FortranTypes::f_int Acol[]
      , const FortranTypes::f_int row_offset
      , const FortranTypes::f_int col_offset
     ) const;

private:

  // //////////////////////////////
  // Private data members

  bool				initialized_;
  size_type			n_;		// Number of rows and columns in G and number of rows in A.
  size_type			m_;		// Number of columns in A
  bool				use_relaxation_;
  value_type			bigM_;
  EPrintMoreOrLess	print_what_;
  ERunTests			test_what_;
  std::ostream		*out_;
  const MatrixOp	*G_;
  const MatrixConvertToSparseFortranCompatible
            *convG_;
  size_type			G_nz_;	// Remember the number of nonzeros of G
  const MatrixOp	*A_;
  const MatrixConvertToSparseFortranCompatible
            *convA_;
  size_type			A_nz_;	// Remember the number of nonzeros of A

  // //////////////////////////////
  // Private member functions

  /** \brief . */
  void assert_matrices_set() const;

  /** \brief . */
  void assert_initialized() const;

  /** \brief Validate the types and sizes of G and A, set the member pointers G_ and A_
    * and return the conversion interfaces convG and convA.
    */
  void validate_and_set_matrices( const MatrixOp& G, const MatrixOp& A );

};	// end class MatrixKKTFullSpaceRelaxed

}	// end namespace ConstrainedOptPack

#endif	// MATRIX_KKT_FULL_SPACE_RELAXED_H
