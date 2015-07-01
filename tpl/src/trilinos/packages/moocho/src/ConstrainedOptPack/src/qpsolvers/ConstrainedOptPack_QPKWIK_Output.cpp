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


#include "Moocho_ConfigDefs.hpp"


#ifdef CONSTRAINED_OPTIMIZATION_PACK_USE_QPKWIK


#include "ConstrainedOptPack_QPKWIK_Output.hpp"
#include "Teuchos_F77_wrappers.h"

namespace QPKWIK_Output {
  std::ostream* out = 0;

  set_output::set_output(std::ostream* _out)
  {	out = _out; }
  set_output::~set_output()
  {	out = 0; }

}	// end namespace QPKWIK_Output

// implementations.
namespace {

// scalar
template<class T>
inline
void output(const char name[], const T& val) {
  *QPKWIK_Output::out << name << " = " << val << "\n";
}

// array
template<class T>
inline
void output(const char name[], const int n, const T array[]) {
  *QPKWIK_Output::out << name << " =\n";
  for(const T* itr = array; itr != array + n; )
    *QPKWIK_Output::out << "\t" << *itr++;
  *QPKWIK_Output::out << "\n";
}

// matrix
template<class T>
inline
void output(const char name[], const int m, const int n, const T matrix[]) {
  *QPKWIK_Output::out << name << " =\n";
  for(int i = 0; i < m; ++i) {
    for(int j = 0; j < n; ++j) {
      *QPKWIK_Output::out << "\t" << matrix[ i + j * m ];
    }
    *QPKWIK_Output::out << "\n";
  }
}

} // end namespace

namespace QPKWIK_Print_Decl {

using FortranTypes::f_int;
using FortranTypes::f_dbl_prec;

// Functions that are called by Fortran QPKWIK

extern "C" {

/// Print input data from QPKWIK
FORTRAN_FUNC_DECL_UL_(void,QPKWIK_PRINT_INPUT,qpkwik_print_input) ( const f_int& N, const f_int& M1
  , const f_int& M2, const f_int& M3, const f_int& M1D, const f_int& M2D
  , const f_int& M3D, const f_dbl_prec GRAD[], const f_dbl_prec Z[]
  , const f_int& LDZ, const f_int IBND[]
  , const f_dbl_prec BL[], const f_dbl_prec BU[], const f_dbl_prec A[], const f_int& LDA
  , const f_dbl_prec YPY[], const f_int& INF, const f_dbl_prec& SMALL
  , const f_dbl_prec& VSMALL, const f_dbl_prec& VLARGE, const f_int& N1
  , const f_int& M12, const f_int& M23, const f_int& M123 )
{
  using QPKWIK_Output::out;
  if(!out) return;

  *out	<< "\n*** Printing QPKWIK input\n";
  output("N",N);
  output("M1",M1);
  output("M2",M2);
  output("M3",M3);
  output("M1D",M1D);
  output("M2D",M2D);
  output("M3D",M3D);
  output("GRAD",N,GRAD);
  output("Z",LDZ,N1,Z);
  output("IBND",M1D,IBND);
  output("BL",M1D,BL);
  output("BU",M3D,BU);
  output("A",LDA,N,A);
//	output("YPY",M1D,YPY);
  output("INF",INF);
  output("SMALL",SMALL);
  output("VSMALL",VSMALL);
  output("VLARGE",VLARGE);
  output("N1",N1);
  output("M12",M12);
  output("M23",M23);
  output("M123",M123);
}

/// Print sparsity info for
FORTRAN_FUNC_DECL_UL_(void,QPKWIK_PRINT_SPARSITY,qpkwik_print_sparsity) ( const f_int& N, const f_int& M2D
  , const f_int& ISPARSE, const f_int ISTART[], const f_int IPOINT[] )
{
  using QPKWIK_Output::out;
  if(!out) return;

  *out	<< "\n*** Printing QPKWIK sparsity data\n";
//	output("ISPARSE",ISPARSE);
  output("ISTART",M2D+1,ISTART);
  output("IPOINT",M2D*N,IPOINT);
}

/// Print iteration data from QPKWIK
FORTRAN_FUNC_DECL_UL_(void,QPKWIK_PRINT_ITERATION_INFO,qpkwik_print_iteration_info) (
    const f_int& CALLING_LABLE,  const f_int& N, const f_int& M1
  , const f_int& M2, const f_int& M3, const f_int& M1D, const f_int& M2D
  , const f_int& M3D, const f_dbl_prec X[]
  , const f_int& NACT, const f_int IACT[], const f_dbl_prec UR[]
  , const f_int IACTSTORE[], const f_dbl_prec Z[], const f_int& LDZ, const f_dbl_prec AINV[]
  , const f_dbl_prec T1[], const f_dbl_prec T2[], const f_dbl_prec R[]
  , const f_dbl_prec XX[], const f_int& IYPY, const f_dbl_prec& EXTRA
  , const f_int& WARM, const f_int& NACTSTORE, const f_dbl_prec& SUMY
  , const f_int& ICHECK, const f_int& I, const f_int& J, const f_int& II
  , const f_dbl_prec& SUM, const f_int& KDROP, const f_int& IFLAG
  , const f_int& KSTART, const f_dbl_prec& SUMNORM, const f_dbl_prec& CVMAX
  , const f_dbl_prec& RES, const f_int& KNEXT, const f_int& IFINISH
  , const f_int& IBEGIN, const f_dbl_prec& TEMP, const f_int& INDEX
  , const f_dbl_prec& PARNEW, const f_int& LFLAG, const f_dbl_prec& SUMA
  , const f_dbl_prec& SUMB, const f_dbl_prec& SUMC, const f_dbl_prec& TEMPA
  , const f_dbl_prec& TEMPB, const f_int& IKNEXT, const f_int& JJ, const f_int& JN
  , const f_dbl_prec& PARINC, const f_dbl_prec& STEP
  , const f_dbl_prec& RATIO, const f_int& ICOUNT, const f_dbl_prec& XMIN
  , const f_dbl_prec& BOTTOM, const f_int& IWARM, const f_int& ITEMP
  , const f_int& ITEMPP )
{
  using QPKWIK_Output::out;
  if(!out) return;

  *out	<< "\n*** QPKWIK Iteration info, CALLING_LABLE = " << CALLING_LABLE << "\n";
  output("X",N,X);
  output("NACT",NACT);
  output("IACT",NACT,IACT);
  output("UR",NACT,UR);
  output("IACTSTORE",NACTSTORE,IACTSTORE);
  output("Z",LDZ,N+1,Z);
  output("AINV",M3D+1,AINV);
  output("T1",N+1,T1);
  output("T2",N+1,T2);
  output("R",(3*(N+1)+(N+1)*(N+1))/2,R);
  output("XX",N,XX);
  output("IYPY",IYPY);
  output("EXTRA",EXTRA);
  output("WARM",WARM);
  output("NACTSTORE",NACTSTORE);
  output("SUMY",SUMY);
  output("ICHECK",ICHECK);
  output("I",I);
  output("J",J);
  output("II",II);
  output("SUM",SUM);
  output("KDROP",KDROP);
  output("IFLAG",IFLAG);
  output("KSTART",KSTART);
  output("SUMNORM",SUMNORM);
  output("CVMAX",CVMAX);
  output("RES",RES);
  output("KNEXT",KNEXT);
  output("IFINISH",IFINISH);
  output("IBEGIN",IBEGIN);
  output("TEMP",TEMP);
  output("INDEX",INDEX);
  output("PARNEW",PARNEW);
  output("LFLAG",LFLAG);
  output("SUMA",SUMA);
  output("SUMB",SUMB);
  output("SUMC",SUMC);
  output("TEMPA",TEMPA);
  output("TEMPB",TEMPB);
  output("IKNEXT",IKNEXT);
  output("JJ",JJ);
  output("JN",JN);
  output("PARINC",PARINC);
  output("STEP",STEP);
  output("RATIO",RATIO);
  output("ICOUNT",ICOUNT);
  output("XMIN",XMIN);
  output("BOTTOM",BOTTOM);
  output("IWARM",IWARM);
  output("ITEMP",ITEMP);
  output("ITEMPP",ITEMPP);
}


}	// end extern "C"


}	// end namespace QPKWIK_Print_Decl


#endif // CONSTRAINED_OPTIMIZATION_PACK_USE_QPKWIK
