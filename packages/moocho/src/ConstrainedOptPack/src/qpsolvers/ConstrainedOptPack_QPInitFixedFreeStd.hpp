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

#ifndef QP_INIT_FIXED_FREE_STD_H
#define QP_INIT_FIXED_FREE_STD_H

#include "ConstrainedOptPack_QPSchur.hpp"

namespace ConstrainedOptPack {
namespace QPSchurPack {

/** \brief General (and flexible) implementation class for a QPSchur QP
 * problem.
 *
 * The basic idea of this class is to just build the QP from its
 * various components in a way that is easy and flexible for the
 * client.  The class will also do consistency testing if asked
 * to.
 */
class QPInitFixedFreeStd : public QP {
public:

  /// Construct uninitialized
  QPInitFixedFreeStd();

  /** \brief Initialize.
   *
   * The pointers and references to the objects pointed to by the
   * arguments to this function must not be modified by the caller.
   * Copies of these objects are not made internally so these
   * objects must remain valid while this object is in use.
   *
   * If the sizes of the arguments do not match up or some consistency
   * test fails then exceptions may be thrown with (hopefully) helpful
   * messages.
   *
   *	@param	g 	[in] vector (size <tt>n</tt>): objective gradient
   *	@param	G	[in] matrix (size <tt>n x n</tt>): objective Hessian
   *	@param	A	[in] matrix (size <tt>n x m</tt>): full rank equality constraints
   *					in <tt>Ko</tt>.  If <tt>A==NULL</tt> then there are no equality constraints
   *					in <tt>Ko</tt> and m will be zero.
   *	@param	n_R	[in] number of initially free variables
   *	@param	i_x_free
   *				[in] array (size <tt>n_R</tt>): <tt>i_x_free[l-1], l = 1...n_R</tt> defines
   *					the matrix <tt>Q_R</tt> as:<br>
   *					<tt>Q_R(:,l) = e(i_x_free[l-1]), l = 1...n_R</tt><br>
   *					The ordering of these indices is significant.  It is allowed
   *                 for <tt>i_x_free == NULL</tt> in which case it will be
   *                 considered to be identity.
   *	@param	i_x_fixed
   *				[in] array (size <tt>n_X = n - n_R</tt>):
   *					<tt>i_x_fixed[l-1], l = 1...n_X</tt> defines the matrix <tt>Q_X</tt> as:<br>
   *					<tt>Q_X(:,l) = e(i_x_fixed[l-1]), l = 1...n_X</tt><br>
   *					The ordering of these indices is significant.
   *	@param	bnd_fixed
   *				[in] array (size <tt>n_X = n - n_R</tt>):
    *					<tt>bnd_fixed[l-1], l = 1...n_X</tt> defines the initial active set as:<br>
   \begin{verbatim}
                      / LOWER : b_X(l) = xL(i_x_fixed[l-1])
    bnd_fixed[l-1] = |  UPPER : b_X(l) = xU(i_x_fixed[l-1])
                      \ EQUALITY : b_X(l) = xL(i) = xU(i) (i = i_x_fixed[l-1])
   \end{verbatim}
   *	@param	b_X	[in] vector (size <tt>n_X = n - n_R</tt>):
   *				Initial varaible bounds (see <tt>bnd_fixed</tt>)
   *	@param	Ko	[in] matrix (size <tt>(n_R+m) x (n_R+m)</tt>):  Initial KKT matrix
   *	@param	fo 	[in] vector (size <tt>n_R + m</tt>): Initial KKT system rhs vector
   *	@param	constraints
   *				[in] Constraints object for the extra constraints
   *					<tt>cL_bar <= A_bar'*x <= cU_bar</tt>
   *	@param	out	[out] If <tt>out!=NULL</tt>, then any warning or error messages will
   *					be printed here.
   *	@param	test_setup
   *				[in] If set to true, then consistency checks will be
   *					made on all the input arguments.  The cost of the
   *					tests will not be too excessive in runtime or
   *					storge costs and do not completly validate everything
   *	@param	waring_tol
   *				[in] Warning tolerance for tests.
   *	@param	error_tol
   *				[in] Error tolerance for tests.  If the relative error
   *					of any test exceeds this limit, then an error
   *					message will be printed to out (if <tt>out!=NULL</tt>) and then
   *					a runtime exception will be thrown.
   *	@param	print_all_warnings
   *				[in] If set to <tt>true</tt>, then any relative errors for tests
   *					that are above <tt>warning_tol</tt> will be printed to
   *					<tt>out</tt> (if <tt>out!= NULL</tt>) (O(<tt>n</tt>) output).
   *					Otherwise, if <tt>false</tt>, then
   *					only the number of violations and the maximum
   *					violation will be printed (O(1) output).
   */
  void initialize(
    const DVectorSlice						&g
    ,const MatrixSymOp					&G
    ,const MatrixOp						*A
    ,size_type								n_R
    ,const size_type						i_x_free[]
    ,const size_type						i_x_fixed[]
    ,const EBounds							bnd_fixed[]
    ,const DVectorSlice						&b_X
    ,const MatrixSymOpNonsing		&Ko
    ,const DVectorSlice						&fo
    ,Constraints							*constraints
    ,std::ostream							*out				= NULL
    ,bool									test_setup			= false
    ,value_type								warning_tol			= 1e-10
    ,value_type								error_tol			= 1e-5
    ,bool									print_all_warnings	= false
    );

  /** @name Overridden from QP */
  //@{ 

  /** \brief . */
  size_type n() const;
  /** \brief . */
  size_type m() const;
  /** \brief . */
  const DVectorSlice g() const;
  /** \brief . */
  const MatrixSymOp& G() const;
  /** \brief . */
  const MatrixOp& A() const;
  /** \brief . */
  size_type n_R() const;
  /** \brief . */
  const x_init_t& x_init() const;
  /** \brief . */
  const l_x_X_map_t& l_x_X_map() const;
  /** \brief . */
  const i_x_X_map_t& i_x_X_map() const;
  /** \brief . */
  const DVectorSlice b_X() const;
  /** \brief . */
  const GenPermMatrixSlice& Q_R() const;
  /** \brief . */
  const GenPermMatrixSlice& Q_X() const;
  /** \brief . */
  const MatrixSymOpNonsing& Ko() const;
  /** \brief . */
  const DVectorSlice fo() const;
  /** \brief . */
  Constraints& constraints();
  /** \brief . */
  const Constraints& constraints() const;

  //@}

private:

  // ///////////////////////////////////
  // Private types

  typedef std::vector<size_type>		row_i_t;
  typedef std::vector<size_type>		col_j_t;
  
  // ///////////////////////////////////
  // Private data members


  size_type				n_;
  size_type				n_R_;
  size_type				m_;
  DVectorSlice				g_;	// will not be modified!
  const MatrixSymOp	*G_;
  const MatrixOp		*A_;	// If NULL not no equalities in Ko
  x_init_t				x_init_;
  l_x_X_map_t				l_x_X_map_;
  i_x_X_map_t				i_x_X_map_;
  DVectorSlice				b_X_;	// will not be modified!
  GenPermMatrixSlice		Q_R_;
  row_i_t					Q_R_row_i_;
  col_j_t					Q_R_col_j_;
  GenPermMatrixSlice		Q_X_;
  row_i_t					Q_X_row_i_;
  col_j_t					Q_X_col_j_;
  const MatrixSymOpNonsing
              *Ko_;
  DVectorSlice				fo_;	// will not be modified
  Constraints				*constraints_;

  // Private member function
  void assert_initialized() const;

};	// end class QPInitFixedFreeStd

}	// end namespace QPSchurPack
}	// end namespace ConstrainedOptPack 

#endif	// QP_INIT_FIXED_FREE_STD_H
