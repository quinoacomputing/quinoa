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
// This is a class for encapsulsting the MA28 package as
// and object so that multiple sparse matrix systems can
// be solved at the same time.

#ifndef SSP_MA28_SOLVER_H
#define SSP_MA28_SOLVER_H

#include "AbstractLinAlgPack_MA28CommonBlockEncap.hpp"

namespace MA28_Cpp {

// Adobt some of the declarations from MA29_CppDecl
using MA28_CppDecl::f_int;
using MA28_CppDecl::f_logical;
using MA28_CppDecl::f_real;
using MA28_CppDecl::f_dbl_prec;

/** \brief MA28 Basic Encapsulation Class.
  *
  * Each object of this class represents a specific MA28 package.
  * Each object encapsulates the common block data for MA28.
  * This class has the same interface functions as MA28 (#ma28ad#
  * , #ma28bd#, and #ma28cd#).  It also has functions for
  * getting and retrieving the common block data from the 
  * MA28 common blocks.
  *
  */
class MA28Solver {
public:

  /// Construct a solver object that is initialized with the default common block data variables.
  MA28Solver();

  /// Construct a solver object that is initialized with the common block data of another solver.
  MA28Solver(const MA28Solver& s);

  // MA28 interface functions

  /** \brief . */
  void ma28ad(const f_int& n, const f_int& nz, f_dbl_prec a[], const f_int& licn
    , f_int irn[], const f_int& lirn, f_int icn[], const f_dbl_prec& u
    , f_int ikeep[], f_int iw[], f_dbl_prec w[], f_int* iflag)
  {	
    set_common_block_data();
    MA28_CppDecl::ma28ad(n,nz,a,licn,irn,lirn,icn,u,ikeep,iw,w,iflag);
    get_common_block_data();
  }

  /** \brief . */
  void ma28bd(const f_int& n, const f_int& nz, f_dbl_prec a[], const f_int& licn
    , const f_int ivect[], const f_int jvect[], const f_int icn[]
    , const f_int ikeep[], f_int iw[], f_dbl_prec w[], f_int* iflag)
  {	
    set_common_block_data();
    MA28_CppDecl::ma28bd(n,nz,a,licn,ivect,jvect,icn,ikeep,iw,w,iflag);
    get_common_block_data();
  }

  /** \brief . */
  void ma28cd(const f_int& n, const f_dbl_prec a[], const f_int& licn, const f_int icn[]
    , const f_int ikeep[], f_dbl_prec rhs[], f_dbl_prec w[], const f_int& mtype)
  {	
    set_common_block_data();
    MA28_CppDecl::ma28cd(n,a,licn,icn,ikeep,rhs,w,mtype);
    get_common_block_data();
  }

  // Common block data setting and retrieval functions

  /** \brief . */
  void		lp(f_int lp)
  {	changed_=true; common_blocks_.ma28ed_.lp=lp; }
  /** \brief . */
  f_int		lp()
  {	return common_blocks_.ma28ed_.lp; }
  /** \brief . */
  void		mp(f_int mp)
  {	changed_=true; common_blocks_.ma28ed_.mp=mp; }
  /** \brief . */
  f_int		mp()
  {	return common_blocks_.ma28ed_.mp; }
  /** \brief . */
  void		lblock(f_logical lblock)
  {	changed_=true; common_blocks_.ma28ed_.lblock=lblock; }
  /** \brief . */
  f_logical	lblock()
  {	return common_blocks_.ma28ed_.lblock; }
  /** \brief . */
  void		grow(f_logical grow)
  {	changed_=true; common_blocks_.ma28ed_.grow=grow; }
  /** \brief . */
  f_logical	grow()
  {	return common_blocks_.ma28ed_.grow; }
  /** \brief . */
  void		eps(f_dbl_prec eps)
  {	changed_=true; common_blocks_.ma28fd_.eps=eps; }
  /** \brief . */
  f_dbl_prec	eps()
  {	return common_blocks_.ma28fd_.eps; }
  /** \brief . */
  void		rmin(f_dbl_prec rmin)
  {	changed_=true; common_blocks_.ma28fd_.rmin=rmin; }
  /** \brief . */
  f_dbl_prec	rmin()
  {	return common_blocks_.ma28fd_.rmin; }
  /** \brief . */
  void		resid(f_dbl_prec resid)
  {	changed_=true; common_blocks_.ma28fd_.resid=resid; }
  /** \brief . */
  f_dbl_prec	resid()
  {	return common_blocks_.ma28fd_.resid; }
  /** \brief . */
  void		irncp(f_int irncp)
  {	changed_=true; common_blocks_.ma28fd_.irncp=irncp; }
  /** \brief . */
  f_int		irncp()
  {	return common_blocks_.ma28fd_.irncp; }
  /** \brief . */
  void		icncp(f_int icncp)
  {	changed_=true; common_blocks_.ma28fd_.icncp=icncp; }
  /** \brief . */
  f_int		icncp()
  {	return common_blocks_.ma28fd_.icncp; }
  /** \brief . */
  void		minirn(f_int minirn)
  {	changed_=true; common_blocks_.ma28fd_.minirn=minirn; }
  /** \brief . */
  f_int		minirn()
  {	return common_blocks_.ma28fd_.minirn; }
  /** \brief . */
  void		minicn(f_int minicn)
  {	changed_=true; common_blocks_.ma28fd_.minicn=minicn; }
  /** \brief . */
  f_int		minicn()
  {	return common_blocks_.ma28fd_.minicn; }
  /** \brief . */
  void		irank(f_int irank)
  {	changed_=true; common_blocks_.ma28fd_.irank=irank; }
  /** \brief . */
  f_int		irank()
  {	return common_blocks_.ma28fd_.irank; }
  /** \brief . */
  void		abort1(f_logical abort1)
  {	changed_=true; common_blocks_.ma28fd_.abort1=abort1; }
  /** \brief . */
  f_logical	abort1()
  {	return common_blocks_.ma28fd_.abort1; }
  /** \brief . */
  void		abort2(f_logical abort2)
  {	changed_=true; common_blocks_.ma28fd_.abort2=abort2; }
  /** \brief . */
  f_logical	abort2()
  {	return common_blocks_.ma28fd_.abort2; }
  /** \brief . */
  void		idisp(f_int val, f_int i)
  {	changed_=true; common_blocks_.ma28gd_.idisp[i]=val; }
  /** \brief . */
  f_int		idisp(f_int i)
  {	return common_blocks_.ma28gd_.idisp[i]; }
  /** \brief . */
  void		tol(f_dbl_prec tol)
  {	changed_=true; common_blocks_.ma28hd_.tol=tol; }
  /** \brief . */
  f_dbl_prec	tol()
  {	return common_blocks_.ma28hd_.tol; }
  /** \brief . */
  void		themax(f_dbl_prec themax)
  {	changed_=true; common_blocks_.ma28hd_.themax=themax; }
  /** \brief . */
  f_dbl_prec	themax()
  {	return common_blocks_.ma28hd_.themax; }
  /** \brief . */
  void		big(f_dbl_prec big)
  {	changed_=true; common_blocks_.ma28hd_.big=big; }
  /** \brief . */
  f_dbl_prec	big()
  {	return common_blocks_.ma28hd_.big; }
  /** \brief . */
  void		dxmax(f_dbl_prec dxmax)
  {	changed_=true; common_blocks_.ma28hd_.dxmax=dxmax; }
  /** \brief . */
  f_dbl_prec	dxmax()
  {	return common_blocks_.ma28hd_.dxmax; }
  /** \brief . */
  void		errmax(f_dbl_prec errmax)
  {	changed_=true; common_blocks_.ma28hd_.errmax=errmax; }
  /** \brief . */
  f_dbl_prec	errmax()
  {	return common_blocks_.ma28hd_.errmax; }
  /** \brief . */
  void		dres(f_dbl_prec dres)
  {	changed_=true; common_blocks_.ma28hd_.dres=dres; }
  /** \brief . */
  f_dbl_prec	dres()
  {	return common_blocks_.ma28hd_.dres; }
  /** \brief . */
  void		cgce(f_dbl_prec cgce)
  {	changed_=true; common_blocks_.ma28hd_.cgce=cgce; }
  /** \brief . */
  f_dbl_prec	cgce()
  {	return common_blocks_.ma28hd_.cgce; }
  /** \brief . */
  void		ndrop(f_int ndrop)
  {	changed_=true; common_blocks_.ma28hd_.ndrop=ndrop; }
  /** \brief . */
  f_int		ndrop()
  {	return common_blocks_.ma28hd_.ndrop; }
  /** \brief . */
  void		maxit(f_int maxit)
  {	changed_=true; common_blocks_.ma28hd_.maxit=maxit; }
  /** \brief . */
  f_int		maxit()
  {	return common_blocks_.ma28hd_.maxit; }
  /** \brief . */
  void		noiter(f_int noiter)
  {	changed_=true; common_blocks_.ma28hd_.noiter=noiter; }
  /** \brief . */
  f_int		noiter()
  {	return common_blocks_.ma28hd_.noiter; }
  /** \brief . */
  void		nsrch(f_int nsrch)
  {	changed_=true; common_blocks_.ma28hd_.nsrch=nsrch; }
  /** \brief . */
  f_int		nsrch()
  {	return common_blocks_.ma28hd_.nsrch; }
  /** \brief . */
  void		istart(f_int istart)
  {	changed_=true; common_blocks_.ma28hd_.istart=istart; }
  /** \brief . */
  f_int		istart()
  {	return common_blocks_.ma28hd_.istart; }
  /** \brief . */
  void		lbig(f_logical lbig)
  {	changed_=true; common_blocks_.ma28hd_.lbig=lbig; }
  /** \brief . */
  f_logical	lbig()
  {	return common_blocks_.ma28hd_.lbig; }

  /// Dump the common block infomation for this solver object.
  void dump_common_blocks(std::ostream& o)
  {	common_blocks_.dump_values(o); }

  /// Copy the state of one solver to another
  MA28Solver& operator=(const MA28Solver& solver)
  {	changed_ = true; common_blocks_ = solver.common_blocks_; return *this; }

  // ///////////////////////////////////
  // Static member functions

  /// Dump the common block infomation for ma28 common blocks
  static void dump_ma28_common_blocks(std::ostream& o)
  {	ma28_common_blocks_.dump_values(o); }

private:

  // ////////////////////////////////////
  // Private member functions

  // Copy the local copy the common block data to MA28 before a MA28 call.
  void set_common_block_data();

  // Retrieve the common block data after a ma28 call.
  void get_common_block_data();
  
  // ///////////////////////////////////
  // Private member data

  // Common block data for this solver object
   MA28CommonBlockStorage common_blocks_;

  // Flag for if the common bock data has changed
  bool changed_;

  // ///////////////////////////////////
  // Static member data

  // Copies of the default values for the 
  // common block data.
  static MA28CommonBlockStorage default_common_blocks_;

  // References to the MA28 common blocks
  static MA28CommonBlockReferences ma28_common_blocks_;

  // Pointer variable who's purpose it to identify
  // what solver object is the current one.
  static MA28Solver* curr_solver_;
  
};

}	// end namespace MA28_Cpp

#endif // SSP_MA28_SOLVER_H
