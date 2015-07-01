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
//
// C declarations for MA28 functions.  These declarations should not have to change
// for different platforms.  As long as the fortran object code uses capitalized
// names for its identifers then the declarations in Teuchos_F77_wrappers.h should be
// sufficent for portability.

#ifndef MA28_CPPDECL_H
#define MA28_CPPDECL_H

#include "Teuchos_F77_wrappers.h"

namespace MA28_CppDecl {

// Declarations that will link to the fortran object file.
// These may change for different platforms

using FortranTypes::f_int;			// INTEGER
using FortranTypes::f_real;			// REAL
using FortranTypes::f_dbl_prec;		// DOUBLE PRECISION
using FortranTypes::f_logical;		// LOGICAL

namespace Fortran {
extern "C" {

// analyze and factorize a matrix
FORTRAN_FUNC_DECL_UL(void,MA28AD,ma28ad) (const f_int& n, const f_int& nz, f_dbl_prec a[], const f_int& licn
  , f_int irn[], const f_int& lirn, f_int icn[], const f_dbl_prec& u, f_int ikeep[], f_int iw[]
  , f_dbl_prec w[], f_int* iflag);
  
// factor using previous analyze
FORTRAN_FUNC_DECL_UL(void,MA28BD,ma28bd) (const f_int& n, const f_int& nz, f_dbl_prec a[], const f_int& licn
  , const f_int ivect[], const f_int jvect[], const f_int icn[], const f_int ikeep[], f_int iw[]
  , f_dbl_prec w[], f_int* iflag);

// solve for rhs using internally stored factorized matrix
FORTRAN_FUNC_DECL_UL(void,MA28CD,ma28cd) (const f_int& n, const f_dbl_prec a[], const f_int& licn, const f_int icn[]
  , const f_int ikeep[], f_dbl_prec rhs[], f_dbl_prec w[], const f_int& mtype);

// /////////////////////////////////////////////////////////////////////////////////////////
// Declare structs that represent the MA28 common blocks.  
// These are the common block variables that are ment to be accessed by the user
// Some are used to set the options of MA28 and others return information
// about the attempts to solve the system.
// I want to provide the access functions that allow all of those common block
// variables that are ment to be accessed by the user to be accessable.
// For each of the common data items there will be a get operation that 
// returns the variable value.  For those items that are ment to be
// set by the user there will also be set operations.

//  COMMON /MA28ED/ LP, MP, LBLOCK, GROW
//  INTEGER LP, MP
//  LOGICAL LBLOCK, GROW
struct MA28ED_struct {
  f_int		lp;
  f_int		mp;
  f_logical	lblock;
  f_logical	grow;
};
extern MA28ED_struct FORTRAN_NAME_UL(MA28ED,ma28ed); // link to fortan common block

//  COMMON /MA28FD/ EPS, RMIN, RESID, IRNCP, ICNCP, MINIRN, MINICN,
// * IRANK, ABORT1, ABORT2
//  INTEGER IRNCP, ICNCP, MINIRN, MINICN, IRANK
//  LOGICAL ABORT1, ABORT2
//  REAL EPS, RMIN, RESID
struct MA28FD_struct {
  f_dbl_prec	eps;
  f_dbl_prec	rmin;
  f_dbl_prec	resid;
  f_int		irncp;
  f_int		icncp;
  f_int		minirn;
  f_int		minicn;
  f_int		irank;
  f_logical	abort1;
  f_logical	abort2;
};
extern MA28FD_struct FORTRAN_NAME_UL(MA28FD,ma28fd); // link to fortan common block


//  COMMON /MA28GD/ IDISP
//  INTEGER IDISP
struct MA28GD_struct {
  f_int		idisp[2];
};
extern MA28GD_struct FORTRAN_NAME_UL(MA28GD,ma28gd); // link to fortan common block

//  COMMON /MA28HD/ TOL, THEMAX, BIG, DXMAX, ERRMAX, DRES, CGCE,
// * NDROP, MAXIT, NOITER, NSRCH, ISTART, LBIG
//  INTEGER NDROP, MAXIT, NOITER, NSRCH, ISTART
//  LOGICAL LBIG
//  REAL TOL, THEMAX, BIG, DXMAX, ERRMAX, DRES, CGCE
struct MA28HD_struct {
  f_dbl_prec	tol;
  f_dbl_prec	themax;
  f_dbl_prec	big;
  f_dbl_prec	dxmax;
  f_dbl_prec	errmax;
  f_dbl_prec	dres;
  f_dbl_prec	cgce;
  f_int		ndrop;
  f_int		maxit;
  f_int		noiter;
  f_int		nsrch;
  f_int		istart;
  f_logical	lbig;
};
extern MA28HD_struct FORTRAN_NAME_UL(MA28HD,ma28hd); // link to fortan common block

//  COMMON /MA30ED/ LP, ABORT1, ABORT2, ABORT3
//  INTEGER LP
//  LOGICAL ABORT1, ABORT2, ABORT3
struct MA30ED_struct {
  f_int		lp;
  f_logical	abort1;
  f_logical	abort2;
  f_logical	abort3;
};
extern MA30ED_struct FORTRAN_NAME_UL(MA30ED,ma30ed); // link to fortan common block

//  COMMON /MA30FD/ IRNCP, ICNCP, IRANK, IRN, ICN
//  INTEGER IRNCP, ICNCP, IRANK, IRN, ICN
struct MA30FD_struct {
  f_int		irncp;
  f_int		icncp;
  f_int		irank;
  f_int		minirn;
  f_int		minicn;
};
extern MA30FD_struct FORTRAN_NAME_UL(MA30FD,ma30fd); // link to fortan common block

//  COMMON /MA30GD/ EPS, RMIN
//  DOUBLE PRECISION EPS, RMIN
struct MA30GD_struct {
  f_dbl_prec	eps;
  f_dbl_prec	rmin;
};
extern MA30GD_struct FORTRAN_NAME_UL(MA30GD,ma30gd); // link to fortan common block

//  COMMON /MA30HD/ RESID
//  DOUBLE PRECISION RESID
struct MA30HD_struct {
  f_dbl_prec	resid;
};
extern MA30HD_struct FORTRAN_NAME_UL(MA30HD,ma30hd); // link to fortan common block

//  COMMON /MA30ID/ TOL, BIG, NDROP, NSRCH, LBIG
//  INTEGER NDROP, NSRCH
//  LOGICAL LBIG
//  DOUBLE PRECISION TOL, BIG
struct MA30ID_struct {
  f_dbl_prec	tol;
  f_dbl_prec	big;
  f_int		ndrop;
  f_int		nsrch;
  f_logical	lbig;
};
extern MA30ID_struct FORTRAN_NAME_UL(MA30ID,ma30id); // link to fortan common block

//  COMMON /MC23BD/ LP,NUMNZ,NUM,LARGE,ABORT
//  INTEGER LP, NUMNZ, NUM, LARGE
//  LOGICAL ABORT
struct MC23BD_struct {
  f_int		lp;
  f_int		numnz;
  f_int		num;
  f_int		large;
  f_logical	abort;
};
extern MC23BD_struct FORTRAN_NAME_UL(MC23BD,mc23bd); // link to fortan common block


} // end extern "C"
} // end namespace Fortran


/* * @name {\bf MA28 C++ Declarations}.
  *
  * These the C++ declarations for MA28 functions and common block data.
  * These declarations will not change for different platforms.
  * All of these functions are in the C++ namespace #MA28_CppDecl#.
  *
  * These functions perform the three phases that are normally associated
  * with solving sparse systems of linear equations; analyze (and factorize)
  * , factorize, and solve.  The MA28 interface uses a coordinate format
  * (aij, i, j) for the sparse matrix.
  * 
  * There are three interface routienes that perform these steps:
  * \begin{description}
  * \item[#ma28ad#] Analyzes and factorizes a sparse matrix stored in #a#, #irn#
  *		, #icn# and returns the factorized matrix data structures in #a#, #icn#
  *		, and #ikeep#.
  * \item[#ma28bd#] Factorizes a matrix with the same sparsity structure that was
  *		previously factorized by #ma28ad#.  Information about the row and column
  *		permutations, fill-in elements etc. from the previous analyze and factorization
  *		is passed in the arguments #icn#, and #ikeep#.  The matrix to be
  *		factorized is passed in #a#, #ivect#, and #jvect# and the non-zero
  *		elements of the factorization are returned in #a#.
  * \item[#ma28cd#] Solves for a dense right hand side (rhs) given a matrix factorized
  *		by #ma28ad# or #ma28bd#.  The rhs is passed in #rhs# and the solution
  *		is returned in #rhs#.  The factorized matrix is passed in by #a#, #icn#
  *		and #ikeep#.  The transposed or the non-transposed system can be solved
  *		for by passing in #mtype != 1# and #mtype == 1# respectively.
  */

// @{
//		begin MA28 C++ Declarations

/* * @name {\bf MA28 / MA30 Common Block Access}.
  *
  * These are references to structures that allow C++ users to set and retrive
  * values of the MA28xD and MA30xD common blocks.   Some of the common block
  * items listed below for MA28xD are also present in MA30xD.  The control
  * parameters (abort1, eps, etc.) for MA28xD are transfered to the equivalent
  * common block variables in the #ma28ad# function but not in any of the other
  * functions.
  *
  * The internal states for MA28, MA30, and MC23 are determined by the 
  * values in these common block variables as there are no #SAVE# variables
  * in any of the functions.  So to use MA28 with more than
  * one sparse matrix at a time you just have to keep copies of these
  * common block variable for each system and then set them when every
  * you want to work with that system agian.  This is very good news.
  *
  * These common block variables are:
  * \begin{description}
  * \item[lp, mp]
  *   Integer: Used by the subroutine as the unit numbers for its warning
  *   and diagnostic messages. Default value for both is 6 (for line
  *   printer output). the user can either reset them to a different
  *   stream number or suppress the output by setting them to zero.
  *   While #lp# directs the output of error diagnostics from the
  *   principal subroutines and internally called subroutines, #mp#
  *   controls only the output of a message which warns the user that he
  *   has input two or more non-zeros a(i), . . ,a(k) with the same row
  *   and column indices.  The action taken in this case is to proceed
  *   using a numerical value of a(i)+...+a(k). in the absence of other
  *   errors, #iflag# will equal -14 on exit.
  * \item[lblock]
  *   Logical: Controls an option of first
  *   preordering the matrix to block lower triangular form (using
  *   harwell subroutine mc23a). The preordering is performed if #lblock#
  *   is equal to its default value of #true# if #lblock# is set to
  *   #false# , the option is not invoked and the space allocated to
  *   #ikeep# can be reduced to 4*n+1.
  * \item[grow]
  *    Logical: If it is left at its default value of
  *   #true# , then on return from ma28a/ad or ma28b/bd, w(1) will give
  *   an estimate (an upper bound) of the increase in size of elements
  *   encountered during the decomposition. If the matrix is well
  *   scaled, then a high value for w(1), relative to the largest entry
  *   in the input matrix, indicates that the LU decomposition may be
  *   inaccurate and the user should be wary of his results and perhaps
  *   increase u for subsequent runs.  We would like to emphasise that
  *   this value only relates to the accuracy of our LU decomposition
  *   and gives no indication as to the singularity of the matrix or the
  *   accuracy of the solution.  This upper bound can be a significant
  *   overestimate particularly if the matrix is badly scaled. If an
  *   accurate value for the growth is required, #lbig# (q.v.) should be
  *   set to #true#
  * \item[eps, rmin]
  *   Double Precision:  If on entry to ma28b/bd, #eps# is less
  *   than one, then #rmin# will give the smallest ratio of the pivot to
  *   the largest element in the corresponding row of the upper
  *   triangular factor thus monitoring the stability of successive
  *   factorizations. if rmin becomes very large and w(1) from
  *   ma28b/bd is also very large, it may be advisable to perform a
  *    new decomposition using ma28a/ad.
  * \item[resid]
  *   Double Precision:  On exit from ma28c/cd gives the value
  *   of the maximum residual over all the equations unsatisfied because
  *   of dependency (zero pivots).
  * \item[irncp,icncp]
  *   Integer:  Monitors the adequacy of "elbow
  *   room" in #irn# and #a#/#icn# respectively. If either is quite large (say
  *   greater than n/10), it will probably pay to increase the size of
  *   the corresponding array for subsequent runs. if either is very low
  *   or zero then one can perhaps save storage by reducing the size of
  *   the corresponding array.
  * \item[minirn, minicn]
  *   Integer: In the event of a
  *   successful return (#iflag# >= 0 or #iflag# = -14) give the minimum size
  *   of #irn# and #a#/#icn# respectively which would enable a successful run
  *   on an identical matrix. On an exit with #iflag# equal to -5, #minicn#
  *   gives the minimum value of #icn# for success on subsequent runs on
  *   an identical matrix. in the event of failure with #iflag# = -6, -4,
  *   -3, -2, or -1, then #minicn# and #minirn# give the minimum value of
  *   #licn# and #lirn# respectively which would be required for a
  *   successful decomposition up to the point at which the failure
  *   occurred.
  * \item[irank]
  *   Integer:  Gives an upper bound on the rank of the matrix.
  * \item[abort1]
  *   Logical:  Default value #true#.  If #abort1# is
  *   set to #false# then ma28a/ad will decompose structurally singular
  *   matrices (including rectangular ones).
  * \item[abort2]
  *   Logical:   Default value #true#.  If #abort2# is
  *   set to #false# then ma28a/ad will decompose numerically singular
  *   matrices.
  * \item[idisp]
  *   Integer[2]:  On output from ma28a/ad, the
  *   indices of the diagonal blocks of the factors lie in positions
  *   idisp(1) to idisp(2) of #a#/#icn#. This array must be preserved
  *   between a call to ma28a/ad and subsequent calls to ma28b/bd,
  *   ma28c/cd or ma28i/id.
  * \item[tol]
  *   Double Precision:  If it is set to a positive value, then any
  *   non-zero whose modulus is less than #tol# will be dropped from the
  *   factorization.  The factorization will then require less storage
  *   but will be inaccurate.  After a run of ma28a/ad with #tol# positive
  *   it is not possible to use ma28b/bd and the user is recommended to
  *   use ma28i/id to obtain the solution.  The default value for #tol# is
  *   0.0.
  * \item[themax]
  *   Double Precision:  On exit from ma28a/ad, it will hold the
  *   largest entry of the original matrix.
  * \item[big]
  *   Double Precision:  If #lbig# has been set to #true#, #big# will hold
  *   the largest entry encountered during the factorization by ma28a/ad
  *   or ma28b/bd.
  * \item[dxmax]
  *   Double Precision:  On exit from ma28i/id, #dxmax# will be set to
  *   the largest component of the solution.
  * \item[errmax]
  *   Double Precision:  On exit from ma28i/id, If #maxit# is
  *   positive, #errmax# will be set to the largest component in the
  *   estimate of the error.
  * \item[dres]
  *   Double Precision:  On exit from ma28i/id, if #maxit# is positive,
  *   #dres# will be set to the largest component of the residual.
  * \item[cgce]
  *   Double Precision:  It is used by ma28i/id to check the
  *   convergence rate.  if the ratio of successive corrections is
  *   not less than #cgce# then we terminate since the convergence
  *   rate is adjudged too slow.
  * \item[ndrop]
  *   Integer:  If #tol# has been set positive, on exit
  *   from ma28a/ad, #ndrop# will hold the number of entries dropped from
  *   the data structure.
  * \item[maxit]
  *   Integer:  It is the maximum number of iterations
  *   performed by ma28i/id. It has a default value of 16.
  * \item[noiter]
  *   Integer:  It is set by ma28i/id to the number of
  *   iterative refinement iterations actually used.
  * \item[nsrch]
  *   Integer:  If #nsrch# is set to a value less than #n#,
  *   then a different pivot option will be employed by ma28a/ad.  This
  *   may result in different fill-in and execution time for ma28a/ad.
  *   If #nsrch# is less than or equal to #n#, the workspace array #iw# can be
  *   reduced in length.  The default value for nsrch is 32768.
  * \item[istart]
  *   Integer:  If #istart# is set to a value other than
  *   zero, then the user must supply an estimate of the solution to
  *   ma28i/id.  The default value for istart is zero.
  * \item[lbig]
  *   Logical:  If #lbig# is set to #true#, the value of the
  *   largest element encountered in the factorization by ma28a/ad or
  *   ma28b/bd is returned in #big#.  setting #lbig# to #true#  will
  *   increase the time for ma28a/ad marginally and that for ma28b/bd
  *   by about 20%.  The default value for #lbig# is #false#.
  */

// @{
//		begin MA28 Common Block Access

// / Common block with members: #lp#, #mp#, #lblock#, #grow#
static MA28ED_struct &ma28ed_cb = FORTRAN_COMMMON_BLOCK_NAME_UL(MA28ED,ma28ed);
// / Common block with members: #eps#, #rmin#, #resid#, #irncp#, #icncp#, #minirc#, #minicn#, #irank#, #abort1#, #abort2#
static MA28FD_struct &ma28fd_cb = FORTRAN_COMMMON_BLOCK_NAME_UL(MA28FD,ma28fd);
// / Common block with members: #idisp#
static MA28GD_struct &ma28gd_cb = FORTRAN_COMMMON_BLOCK_NAME_UL(MA28GD,ma28gd);
// / Common block with members: #tol#, #themax#, #big#, #bxmax#, #errmax#, #dres#, #cgce#, #ndrop#, #maxit#, #noiter#, #nsrch#, #istart#, #lbig#
static MA28HD_struct &ma28hd_cb = FORTRAN_COMMMON_BLOCK_NAME_UL(MA28HD,ma28hd);
// / Common block with members: #lp#, #abort1#, #abort2#, #abort3#
static MA30ED_struct &ma30ed_cb = FORTRAN_COMMMON_BLOCK_NAME_UL(MA30ED,ma30ed);
// / Common block with members: #irncp#, #icncp#, #irank#, #irn#, #icn#
static MA30FD_struct &ma30fd_cb = FORTRAN_COMMMON_BLOCK_NAME_UL(MA30FD,ma30fd);
// / Common block with members: #eps#, #rmin#
static MA30GD_struct &ma30gd_cb = FORTRAN_COMMMON_BLOCK_NAME_UL(MA30GD,ma30gd);
// / Common block with members: #resid#
static MA30HD_struct &ma30hd_cb = FORTRAN_COMMMON_BLOCK_NAME_UL(MA30HD,ma30hd);
// / Common block with members: #tol#, #big#, #ndrop#, #nsrch#, #lbig#
static MA30ID_struct &ma30id_cb = FORTRAN_COMMMON_BLOCK_NAME_UL(MA30ID,ma30id);
// / Common block with members: #lp#, #numnz#, #num#, #large#, #abort#
static MC23BD_struct &mc23bd_cb = FORTRAN_COMMMON_BLOCK_NAME_UL(MC23BD,mc23bd);

  // The reason that these are declared static is because I need to
  // make sure that these references are initialized before they are
  // used in other global defintions.  This means that every translation
  // unit will have their own copy of this data.  To reduce this code
  // blot you could declair them as pointers then set then using the 
  // trick of an initialization class (Myers).

//		end MA28 Common Block Access
// @}

// /
/* * Analyze and factor a sparse matrix.
  *
  * This function analyzes (determines row and column pivots to minimize
  * fill-in and result in a better conditioned factorization) then factors
  * (calculates upper and lower triangular factors for the determined
  * row and column pivots) the permuted system.  The function takes a sparse
  * matrix stored in coordinate form ([Aij, i, j] => [#a#, #irn#, #icn#]) and
  * factorizes it.  On entry, the first #nz# elements of #a#, #irn#, and
  * #icn# hold the matrix elements.  The remaining entires of #a#, #irn#
  * , and #icn# hold the fill-in entries after the matrix factorization is
  * complete.
  *
  * The amount of fill-in is influenced by #u#.  A value of #u = 1.0# gives partial
  * pivoting where fill-in is sacrificed for the sake of numerical stability and
  * #u = 0.0# gives pivoting that strictly minimizes fill-in.
  *
  * The parameters #ikeep#, #iw#, and #w# are used as workspace but
  * #ikeep# contains important information about the factoriation
  * on return.
  *
  * The parameter #iflag# is used return error information about the attempt
  * to factorized the system.
  *
  * @param	n	[input] Order of the system being solved.
  * @param	nz  [input] Number of non-zeros of the input matrix (#nz >= n#).
  *				The ratio #nz/(n*n)# equals the sparsity fraction for the matrix.
  * @param	a	[input/output] Length = #licn#.  The first #nz# entries hold the
  *				non-zero entries of the input matrix on input and
  *				the non-zero entries of the factorized matrix on exit.
  * @param	licn	[input] length of arrays #a# and #icn#.  This
  *					is the total amount of storage advalable for the
  *					non-zero entries of the factorization of #a#.
  *					therefore #licn# must be greater than #nz#. How
  *					much greater depends on the amount of fill-in.
  * @param	irn	[input/modifed] Length = #ircn#.  The first #nz# entries hold
  *				the row indices of the matrix #a# on input.
  * @param	ircn	[input] Lenght of irn.
  * @param	icn	[input/output] Length = #licn#.  Holds column indices of #a# on 
  *				entry and the column indices of the reordered
  *				#a# on exit.
  * @param	u	[input] Controls partial pivoting.
  *				\begin{description}
  *				\item[#* u >= 1.0#]
  *					Uses partial pivoting for maximum numerical stability
  *					at the expense of some extra fill-in.
  *				\item[#* 1.0 < u < 0.0#]
  *					Balances numerical stability and fill-in with #u# near 1.0
  *					favoring stability and #u# near 0.0 favoring less fill-in.
  *				\item[#* u <= 0.0#]
  *					Determines row and column pivots to minimize fill-in
  *					irrespective of the numerical stability of the 
  *					resulting factorization.
  *				\end{description}
  * @param	ikeep	[output] Length = 5 * #n#.  On exist contains information about 
  *					the factorization.
  *					\begin{description}
  *					\item[#* #ikeep(:,1)]
  *						Holds the total length of the part of row i in
  *						the diagonal block.
  *					\item[#* #ikeep(:,2)]
  *						Holds the row pivots.  Row #ikeep(i,2)# of the
  *						input matrix is the ith row of the pivoted matrix
  *						which is factorized.
  *					\item[#* #ikeep(:,3)]
  *						Holds the column pivots.  Column #ikeep(i,3)# of the
  *						input matrix is the ith column of the pivoted matrix
  *						which is factorized.
  *					\item[#* #ikeep(:,4)]
  *						Holds the length of the part of row i in the L
  *						part of the L/U decomposition.
  *					\item[#* #ikeep(:,5)]
  *						Holds the length of the part of row i in the
  *						off-diagonal blocks.  If there is only one
  *						diagonal block, #ikeep(i,5)# is set to -1.
  *					\end{description}
  * @param	iw	[] Length = #8*n#.  Integer workspace.
  * @param	w	[] Length = #n#.  Real workspace.
  * @param	iflag	[output] Used to return error condtions.
  *					\begin{description}
  *					\item[#* >= 0#] Success.
  *					\item[#* < 0#] Some error has occured.
  *					\end{description}
  */
inline void ma28ad(const f_int& n, const f_int& nz, f_dbl_prec a[], const f_int& licn
  , f_int irn[], const f_int& lirn, f_int icn[], const f_dbl_prec& u, f_int ikeep[], f_int iw[]
  , f_dbl_prec w[], f_int* iflag)
{	Fortran::FORTRAN_FUNC_CALL_UL(MA28AD,ma28ad) (n,nz,a,licn,irn,lirn,icn,u,ikeep,iw,w,iflag);	}

// /
/* * Factor a sparse matrix using previous analyze pivots.
  *
  * This function uses the pivots determined from a previous factorization
  * to factorize the matrix #a# again.  The assumption is that the 
  * sparsity structure of #a# has not changed but only its numerical 
  * values.  It is therefore possible that the refactorization may 
  * become unstable.
  *
  * The matrix to be refactorized on order #n# with #nz# non-zero elements
  * is input in coordinate format in #a#, #ivect# and #jvect#. 
  *
  * Information about the factorization is contained in the #icn# and
  * #ikeep# arrays returned from \Ref{ma28ad}.
  *
  * @param	n	[input] Order of the system being solved.
  * @param	nz  [input] Number of non-zeros of the input matrix
  * @param	a	[input/output] Length = #licn#.  The first #nz# entries hold the
  *				non-zero entries of the input matrix on input and
  *				the non-zero entries of the factorized matrix on exit.
  * @param	licn	[input] length of arrays #a# and #icn#.  This
  *					is the total amount of storage avalable for the
  *					non-zero entries of the factorization of #a#.
  *					therefore #licn# must be greater than #nz#. How
  *					much greater depends on the amount of fill-in.
  * @param	icn	[input] Length = #licn#.  Same array output from #ma28ad#.
  *				It contains information about the analyze step.
  * @param	ikeep	[input] Length = 5 * #n#.  Same array output form #ma28ad#.
  *					It contains information about the analyze step.
  * @param	iw	[] Length = #8*n#.  Integer workspace.
  * @param	w	[] Length = #n#.  Real workspace.
  * @param	iflag	[output] Used to return error condtions.
  *					\begin{description}
  *					\item[#* >= 0#] Success
  *					\item[#* < 0#] Some error has occured.
  *					\end{description}
  */
inline void ma28bd(const f_int& n, const f_int& nz, f_dbl_prec a[], const f_int& licn
  , const f_int ivect[], const f_int jvect[], const f_int icn[], const f_int ikeep[], f_int iw[]
  , f_dbl_prec w[], f_int* iflag)
{	Fortran::FORTRAN_FUNC_CALL_UL(MA28BD,ma28bd) (n,nz,a,licn,ivect,jvect,icn,ikeep,iw,w,iflag);	}

// /
/* * Solve for a rhs using a factorized matrix.
  *
  * This function solves for a rhs given a matrix factorized by
  * #ma28ad# or #ma28bd#.  The right hand side (rhs) is passed
  * in in #rhs# and the solution is return in #rhs#.  The 
  * factorized matrix is passed in in #a#, #icn#, and #ikeep#
  * which were set by #ma28ad# and/or #ma28bd#. The 
  *
  * The matrix or its transpose can be solved for by selecting
  * #mtype == 1# or #mtype != 1# respectively.
  *
  * @param	n	[input] Order of the system being solved.
  * @param	a	[input] Length = #licn#.  Contains the non-zero
  *				elements of the factorized matrix.
  * @param	licn	[input] length of arrays #a# and #icn#.  This
  *					is the total amount of storage avalable for the
  *					non-zero entries of the factorization of #a#.
  *					therefore #licn# must be greater than #nz#. How
  *					much greater depends on the amount of fill-in.
  * @param	icn	[input] Length = #licn#.  Same array output from #ma28ad#.
  *				It contains information about the analyze step.
  * @param	ikeep	[input] Length = 5 * #n#.  Same array output form #ma28ad#.
  *					It contains information about the analyze step.
  * @param	w	[] Length = #n#.  Real workspace.
  * @param	mtype	[input] Instructs to solve using the matrix or its transpoze.
  *					\begin{description}
  *					\item[#* mtype == 1#] Solve using the non-transposed matrix.
  *					\item[#* mtype != 1#] Solve using the transposed matrix.
  *					\end{description} 
  */
inline void ma28cd(const f_int& n, const f_dbl_prec a[], const f_int& licn, const f_int icn[]
  , const f_int ikeep[], f_dbl_prec rhs[], f_dbl_prec w[], const f_int& mtype)
{	Fortran::FORTRAN_FUNC_CALL_UL(MA28CD,ma28cd) (n,a,licn,icn,ikeep,rhs,w,mtype);	}

//		end MA28 C++ Declarations
// @}

} // end namespace MA28_CDecl

#endif // MA28_CPPDECL_H
