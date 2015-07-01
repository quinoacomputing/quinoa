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

#include <math.h>

#include <limits>

#include "ConstrainedOptPack_QPSolverRelaxedTester.hpp"
#include "AbstractLinAlgPack_MatrixSymOp.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "AbstractLinAlgPack_VectorAuxiliaryOps.hpp"

namespace {

//
const char* solution_type_str( ConstrainedOptPack::QPSolverStats::ESolutionType solution_type )
{
  typedef ConstrainedOptPack::QPSolverStats qpst;
  switch( solution_type ) {
  
  case qpst::OPTIMAL_SOLUTION:
    return "OPTIMAL_SOLUTION";
  case qpst::PRIMAL_FEASIBLE_POINT:
    return "PRIMAL_FEASIBLE_POINT";
  case qpst::DUAL_FEASIBLE_POINT:
    return "DUAL_FEASIBLE_POINT";
  case qpst::SUBOPTIMAL_POINT:
    return "SUBOPTIMAL_POINT";
  default:
    TEUCHOS_TEST_FOR_EXCEPT(true);
  }
  return "";	// will never be executed.
}

/* ToDo: Update this code!

// Compute the scaled complementarity conditions.
// 
// If uplo == upper then:
// 
//                / gamma(i) * constr_resid(i) / ( 1 + |constr(i)| + opt_scale ), for gamma(i) > 0
// comp_err(i) = |
//                \ 0 otherwise
// 
// If uplo == lower then:
// 
//                / gamma(i) * constr_resid(i) /  ( 1 + |constr(i)| + opt_scale ), for gamma(i) < 0
// comp_err(i) = |
//                \ 0 otherwise
// 
// 
void set_complementarity(
  const AbstractLinAlgPack::SpVector	&gamma
  ,const DenseLinAlgPack::DVectorSlice		&constr_resid
  ,const DenseLinAlgPack::DVectorSlice     &constr
  ,const DenseLinAlgPack::value_type      opt_scale
  ,BLAS_Cpp::Uplo					uplo
  ,DenseLinAlgPack::DVector			 	*comp_err
  )
{
  TEUCHOS_TEST_FOR_EXCEPT( !(  gamma.size() == constr_resid.size() && gamma.size() == constr.size()  ) );
  comp_err->resize( gamma.size() );
  *comp_err = 0.0;
  const AbstractLinAlgPack::SpVector::difference_type o = gamma.offset();
  if( gamma.nz() ) {
    for( AbstractLinAlgPack::SpVector::const_iterator itr = gamma.begin(); itr != gamma.end(); ++itr ) {
      const DenseLinAlgPack::size_type i = itr->indice() + o;
      if( itr->value() > 0 && uplo == BLAS_Cpp::upper )
        (*comp_err)(i) = itr->value() * constr_resid(i) / ( 1.0 + ::fabs(constr(i)) + opt_scale );
      else if( itr->value() < 0 && uplo == BLAS_Cpp::lower )
        (*comp_err)(i) = itr->value() * constr_resid(i) / ( 1.0 + ::fabs(constr(i)) + opt_scale );
    }
  }
}

*/

// Handle the error reporting
void handle_error(
  std::ostream                            *out
  ,const AbstractLinAlgPack::value_type   err
  ,const char                             err_name[]
  ,const AbstractLinAlgPack::value_type   error_tol
  ,const char                             error_tol_name[]
  ,const AbstractLinAlgPack::value_type   warning_tol
  ,const char                             warning_tol_name[]
  ,bool                                   *test_failed
  )
{
  if( err >= error_tol ) {
    if(out)
      *out << "\n" << err_name << " = " << err << " >= " << error_tol_name << " = " << error_tol << std::endl;
    *test_failed = true;
  }
  else if( err >= warning_tol ) {
    if(out)
      *out << "\n" << err_name << " = " << err << " >= " << warning_tol_name << " = " << warning_tol << std::endl;
  }
}

}	// end namespace

namespace ConstrainedOptPack {

// public

QPSolverRelaxedTester::QPSolverRelaxedTester(
  value_type   opt_warning_tol
  ,value_type  opt_error_tol
  ,value_type  feas_warning_tol
  ,value_type  feas_error_tol
  ,value_type  comp_warning_tol
  ,value_type  comp_error_tol
  )
  :opt_warning_tol_(opt_warning_tol)
  ,opt_error_tol_(opt_error_tol)
  ,feas_warning_tol_(feas_warning_tol)
  ,feas_error_tol_(feas_error_tol)
  ,comp_warning_tol_(comp_warning_tol)
  ,comp_error_tol_(comp_error_tol)
{}

bool QPSolverRelaxedTester::check_optimality_conditions(
  QPSolverStats::ESolutionType solution_type
  ,const value_type infinite_bound
  ,std::ostream* out, bool print_all_warnings, bool print_vectors
  ,const Vector& g, const MatrixSymOp& G
  ,value_type etaL
  ,const Vector& dL, const Vector& dU
  ,const MatrixOp& E, BLAS_Cpp::Transp trans_E, const Vector& b
  ,const Vector& eL, const Vector& eU
  ,const MatrixOp& F, BLAS_Cpp::Transp trans_F, const Vector& f
  ,const value_type* obj_d
  ,const value_type* eta, const Vector* d
  ,const Vector* nu
  ,const Vector* mu, const Vector* Ed
  ,const Vector* lambda, const Vector* Fd
  )
{
  return check_optimality_conditions(
    solution_type,infinite_bound,out,print_all_warnings,print_vectors
    ,g,G,etaL,&dL,&dU,&E,trans_E,&b,&eL,&eU,&F,trans_F,&f
    ,obj_d,eta,d,nu,mu,Ed,lambda,Fd);
}

bool QPSolverRelaxedTester::check_optimality_conditions(
  QPSolverStats::ESolutionType solution_type
  ,const value_type infinite_bound
  ,std::ostream* out, bool print_all_warnings, bool print_vectors
  ,const Vector& g, const MatrixSymOp& G
  ,value_type etaL
  ,const Vector& dL, const Vector& dU
  ,const MatrixOp& E, BLAS_Cpp::Transp trans_E, const Vector& b
  ,const Vector& eL, const Vector& eU
  ,const value_type* obj_d
  ,const value_type* eta, const Vector* d
  ,const Vector* nu
  ,const Vector* mu, const Vector* Ed
  )
{
  return check_optimality_conditions(
    solution_type,infinite_bound,out,print_all_warnings,print_vectors
    ,g,G,etaL,&dL,&dU,&E,trans_E,&b,&eL,&eU,NULL,BLAS_Cpp::no_trans,NULL
    ,obj_d,eta,d,nu,mu,Ed,NULL,NULL);
}

bool QPSolverRelaxedTester::check_optimality_conditions(
  QPSolverStats::ESolutionType solution_type
  ,const value_type infinite_bound
  ,std::ostream* out, bool print_all_warnings, bool print_vectors
  ,const Vector& g, const MatrixSymOp& G
  ,value_type etaL
  ,const Vector& dL, const Vector& dU
  ,const MatrixOp& F, BLAS_Cpp::Transp trans_F, const Vector& f
  ,const value_type* obj_d
  ,const value_type* eta, const Vector* d
  ,const Vector* nu
  ,const Vector* lambda, const Vector* Fd
  )
{
  return check_optimality_conditions(
    solution_type,infinite_bound,out,print_all_warnings,print_vectors
    ,g,G,etaL,&dL,&dU,NULL,BLAS_Cpp::no_trans,NULL,NULL,NULL,&F,trans_F,&f
    ,obj_d,eta,d,nu,NULL,NULL,lambda,Fd );
}

bool QPSolverRelaxedTester::check_optimality_conditions(
  QPSolverStats::ESolutionType solution_type
  ,const value_type infinite_bound
  ,std::ostream* out, bool print_all_warnings, bool print_vectors
  ,const Vector& g, const MatrixSymOp& G
  ,const Vector& dL, const Vector& dU
  ,const value_type* obj_d
  ,const Vector* d
  ,const Vector* nu
  )
{
  return check_optimality_conditions(
    solution_type,infinite_bound,out,print_all_warnings,print_vectors
    ,g,G,0.0,&dL,&dU,NULL,BLAS_Cpp::no_trans,NULL,NULL,NULL,NULL,BLAS_Cpp::no_trans,NULL
    ,obj_d,NULL,d,nu,NULL,NULL,NULL,NULL);
}

bool QPSolverRelaxedTester::check_optimality_conditions(
  QPSolverStats::ESolutionType solution_type
  ,const value_type infinite_bound
  ,std::ostream* out, bool print_all_warnings, bool print_vectors
  ,const Vector& g, const MatrixSymOp& G
  ,value_type etaL
  ,const Vector* dL, const Vector* dU
  ,const MatrixOp* E, BLAS_Cpp::Transp trans_E, const Vector* b
  ,const Vector* eL, const Vector* eU
  ,const MatrixOp* F, BLAS_Cpp::Transp trans_F, const Vector* f
  ,const value_type* obj_d
  ,const value_type* eta, const Vector* d
  ,const Vector* nu
  ,const Vector* mu, const Vector* Ed
  ,const Vector* lambda, const Vector* Fd
  )
{
  QPSolverRelaxed::validate_input(
    infinite_bound,g,G,etaL,dL,dU
    ,E,trans_E,b,eL,eU,F,trans_F,f
    ,obj_d,eta,d,nu,mu,Ed,lambda,Fd);

  return imp_check_optimality_conditions(
    solution_type,infinite_bound
    ,out,print_all_warnings,print_vectors,g,G,etaL,dL,dU
    ,E,trans_E,b,eL,eU,F,trans_F,f
    ,obj_d,eta,d,nu,mu,Ed,lambda,Fd);
}

// protected

bool QPSolverRelaxedTester::imp_check_optimality_conditions(
  QPSolverStats::ESolutionType solution_type
  ,const value_type infinite_bound
  ,std::ostream* out, bool print_all_warnings, bool print_vectors
  ,const Vector& g, const MatrixSymOp& G
  ,value_type etaL
  ,const Vector* dL, const Vector* dU
  ,const MatrixOp* E, BLAS_Cpp::Transp trans_E, const Vector* b
  ,const Vector* eL, const Vector* eU
  ,const MatrixOp* F, BLAS_Cpp::Transp trans_F, const Vector* f
  ,const value_type* obj_d
  ,const value_type* eta, const Vector* d
  ,const Vector* nu
  ,const Vector* mu, const Vector* Ed
  ,const Vector* lambda, const Vector* Fd
  )
{
  using std::endl;
  using BLAS_Cpp::trans_not;
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  using BLAS_Cpp::upper;
  using BLAS_Cpp::lower;
  using LinAlgOpPack::sum;
  using LinAlgOpPack::dot;
  using LinAlgOpPack::Vt_S;
  using LinAlgOpPack::V_VmV;
  using LinAlgOpPack::Vp_StV;
  using LinAlgOpPack::V_MtV;
  using LinAlgOpPack::Vp_MtV;
  using LinAlgOpPack::V_StMtV;
  using LinAlgOpPack::Vp_V;
  using AbstractLinAlgPack::max_element;
  typedef QPSolverStats qps_t;

  bool test_failed = false;

  const size_type
    nd = d->dim();

  const value_type
    really_big_error_tol = std::numeric_limits<value_type>::max();

  value_type opt_scale = 0.0;
  VectorSpace::vec_mut_ptr_t
    t_d = d->space().create_member(),
    u_d = d->space().create_member();

  value_type
    err = 0,
    d_norm_inf,	// holds ||d||inf
    e_norm_inf;	// holds ||e||inf

  if(out)
    *out
      << "\n*** Begin checking QP optimality conditions ***\n"
      << "\nThe solution type is " << solution_type_str(solution_type) << endl;

  bool force_opt_error_check
    = solution_type==qps_t::OPTIMAL_SOLUTION || solution_type==qps_t::DUAL_FEASIBLE_POINT;
  const bool force_inequality_error_check
    = solution_type==qps_t::OPTIMAL_SOLUTION || solution_type==qps_t::PRIMAL_FEASIBLE_POINT;
  const bool force_equality_error_check
    = solution_type!=qps_t::SUBOPTIMAL_POINT;
  const bool force_complementarity_error_check
    = solution_type!=qps_t::SUBOPTIMAL_POINT;

  const char sep_line[] = "\n--------------------------------------------------------------------------------\n";

  ////////////////////////////
  // Checking d(L)/d(d) = 0
  if(out)
    *out
      << sep_line
      << "\nChecking d(L)/d(d) = g + G*d + nu + op(E)'*mu - op(F)'*lambda == 0 ...\n";

  if(out && !force_opt_error_check)
    *out
      << "The optimality error tolerance will not be enforced ...\n";

  opt_scale = 0.0;

  *u_d = g;
  opt_scale += g.norm_inf();

  if(out) {
    *out << "||g||inf = " << g.norm_inf() << endl;
  }
  
  V_MtV( t_d.get(), G, no_trans, *d );
  Vp_V( u_d.get(), *t_d );
  opt_scale += t_d->norm_inf();

  if(out) {
    *out << "||G*d||inf = " << t_d->norm_inf() << endl;
    if(print_vectors)
      *out << "g + G*d =\n" << *u_d;
  }

  if( nu ) {
    Vp_V( u_d.get(), *nu );
    opt_scale += nu->norm_inf();
    if(out)
      *out << "||nu||inf = " << nu->norm_inf() << endl;
  }

  if(E) {
    V_MtV( t_d.get(), *E, trans_not(trans_E), *mu );
    Vp_V( u_d.get(), *t_d );
    opt_scale += t_d->norm_inf();
    if(out) {
      *out << "||op(E)'*mu||inf = " << t_d->norm_inf() << endl;
      if(print_vectors)
        *out << "op(E)'*mu =\n" << *t_d;
    }
  }			

  if(F) {
    V_MtV( t_d.get(), *F, trans_not(trans_F), *lambda );
    Vp_V( u_d.get(), *t_d );
    opt_scale += t_d->norm_inf();
    if(out) {
      *out << "||op(F)'*lambda||inf = " << t_d->norm_inf() << endl;
      if(print_vectors)
        *out << "op(F)'*lambda =\n" << *t_d;
    }
  }

  if( *eta > etaL ) { // opt_scale + |(eta - etaL) * (b'*mu + f'*lambda)|
    const value_type
      term = ::fabs( (*eta - etaL) * (E ? dot(*b,*mu) : 0.0) + (F ? dot(*f,*lambda) : 0.0 ) );
    if(out) {
      *out << "|(eta - etaL) * (b'*mu + f'*lambda)| = " << term << endl;
    }
    opt_scale += term;
  }

  if(out && print_vectors)
    *out
      << "g + G*d + nu + op(E)'*mu - op(F)'*lambda =\n" << *u_d;

  Vt_S( u_d.get(), 1.0/(1.0+opt_scale) );

  err = sum( *u_d );

  if(out)
    *out
      << "\nopt_scale = " << opt_scale << endl
      << "opt_err = sum( | g + G*d + nu + op(E)'*mu - op(F)'*lambda | / (1+opt_scale) ) / nd\n"
      << "        = " << err << " / " << nd << " = " << (err/nd) << endl;

  err *= nd;

  if( force_opt_error_check ) {
    if( err >= opt_error_tol() ) {
      if(out)
        *out << "\nopt_err = " << err << " >= opt_error_tol = " << opt_error_tol() << endl;
      test_failed = true;
    }
    else if( err >= opt_warning_tol() ) {
      if(out)
        *out << "\nopt_err = " << err << " >= opt_error_tol = " << opt_error_tol() << endl;
    }
  }

  if(out) {
    *out
      << sep_line
      << "\nTesting feasibility of the constraints and the complementarity conditions ...\n";
    if(!force_inequality_error_check)
      *out
        << "The inequality feasibility error tolerance will not be enforced ...\n";
    if(!force_equality_error_check)
      *out
        << "The equality feasibility error tolerance will not be enforced ...\n";
    if(!force_complementarity_error_check)
      *out
        << "The complementarity conditions error tolerance will not be enforced ...\n";
  }

  /////////////////////
  // etaL - eta
  if(out)
    *out
      << sep_line
      << "\nChecking etaL - eta <= 0 ...\n";
  if(out)
    *out
      << "etaL - eta = " << (etaL - (*eta)) << endl;
  if( etaL - (*eta) > feas_warning_tol() ) {
    if(out)
      *out
        << "Warning, etaL - eta = " << etaL << " - " << (*eta)
        << " = " << (etaL - (*eta)) << " >  feas_warning_tol = "
        <<  feas_warning_tol() << endl;
  }
  if( force_inequality_error_check && etaL - (*eta) > feas_error_tol() ) {
    if(out)
      *out
        << "Error, etaL - eta = " << etaL << " - " << (*eta)
        << " = " << (etaL - (*eta)) << " >  feas_error_tol = "
        <<  feas_error_tol() << endl;
    test_failed = true;
  } 

  d_norm_inf = d->norm_inf();

  if(dL) {

    ///////////////////////////////////
    // dL - d <= 0
    if(out)
      *out
        << sep_line
        << "\nChecking dL - d <= 0 ...\n";
    V_VmV( u_d.get(), *dL, *d );
    if(out && print_vectors)
      *out
        << "dL - d =\n" << *u_d;
    Vt_S( u_d.get(), 1.0/(1.0+d_norm_inf) );

    err = max_element(*u_d);
    if(out)
      *out
        << "\nmax(dL-d) = " << err << endl;
    if( force_inequality_error_check )
      handle_error(
        out,err,"max(dU-d)",feas_error_tol(),"feas_error_tol"
        ,feas_warning_tol(),"feas_error_tol",&test_failed
        );

    // ToDo: Update below code!
/*
    if(out)
      *out
        << sep_line
        << "\nChecking nuL(i) * (dL - d)(i) = 0  ...\n";
    set_complementarity( *nu, u(), *d, opt_scale, lower, &c );
    if(out && print_vectors)
      *out
        << "nuL(i) * (dL - d)(i) / ( 1 + |d(i)| + opt_scale ) =\n" << c();
    if(out) {
      *out
        << "Comparing:\n"
        << "u(i) = nuL(i) * (dL - d)(i) / ( 1 + |d(i)| + opt_scale ), v = 0 ...\n";
    }
    if(!comp_v.comp( c(), 0.0, opt_warning_tol()
      , force_complementarity_error_check ? comp_error_tol() : really_big_error_tol
      , print_all_warnings, out )) test_failed = true;
*/

    ///////////////////////////////////
    // d - dU <= 0
    if(out)
      *out
        << sep_line
        << "\nChecking d - dU <= 0 ...\n";
    V_VmV( u_d.get(), *d, *dU );
    if(out && print_vectors)
      *out
        << "d - dU =\n" << *u_d;
    Vt_S( u_d.get(), 1.0/(1.0+d_norm_inf) );

    err = max_element(*u_d);
    if(out)
      *out
        << "\nmax(d-dU) = " << err << endl;
    if( force_inequality_error_check )
      handle_error(
        out,err,"max(d-dU)",feas_error_tol(),"feas_error_tol"
        ,feas_warning_tol(),"feas_error_tol",&test_failed
        );
    // ToDo: Update below code!
/*
    if(out)
      *out
        << sep_line
        << "\nChecking nuU(i) * (d - dU)(i) = 0  ...\n";
    set_complementarity( *nu, u(), *d, opt_scale, upper, &c );
    if(out && print_vectors)
      *out
        << "nuU(i) * (d - dU)(i) / ( 1 + |d(i)| + opt_scale ) =\n" << c();
    if(out) {
      *out
        << "Comparing:\n"
        << "u(i) = nuU(i) * (dL - d)(i) / ( 1 + |d(i)| + opt_scale ), v = 0 ...\n";
    }
    if(!comp_v.comp( c(), 0.0, opt_warning_tol()
             , force_complementarity_error_check ? comp_error_tol() : really_big_error_tol
             , print_all_warnings, out )) test_failed = true;
*/
  }

  if( E ) {

    ///////////////////////////////////
    // e = op(E)*d + b*eta
    if(out)
      *out
        << sep_line
        << "\nComputing e = op(E)*d - b*eta ...\n";
    VectorSpace::vec_mut_ptr_t
      e   = ( trans_E == no_trans ? E->space_cols() : E->space_rows() ).create_member(),
      t_e = e->space().create_member();
    V_MtV( e.get(), *E, trans_E, *d );
    Vp_StV( e.get(), -(*eta), *b );
    e_norm_inf = e->norm_inf();
    if(out && print_vectors)
      *out
        << "e = op(E)*d - b*eta  =\n" << *e;

    ///////////////////////////////////
    // eL - e <= 0
    if(out)
      *out
        << sep_line
        << "\nChecking eL - e <= 0 ...\n";
    V_VmV( t_e.get(), *eL, *e );
    if(out && print_vectors)
      *out
        << "eL - e =\n" << *t_e;
    Vt_S( t_e.get(), 1.0/(1.0+e_norm_inf) );

    err = max_element(*t_e);
    if(out)
      *out
        << "\nmax(eL-e) = " << err << endl;
    if( force_inequality_error_check )
      handle_error(
        out,err,"max(eL-e)",feas_error_tol(),"feas_error_tol"
        ,feas_warning_tol(),"feas_error_tol",&test_failed
        );
    // ToDo: Update below code!
/*
    if(out)
      *out
        << sep_line
        << "\nChecking muL(i) * (eL - e)(i) = 0  ...\n";
    set_complementarity( *mu, u(), e(), opt_scale, lower, &c );
    if(out && print_vectors)
      *out
        << "muL(i) * (eL - e)(i) / ( 1 + |e(i)| + opt_scale ) =\n" << c();
    if(out) {
      *out
        << "Comparing:\n"
        << "u(i) = muL(i) * (eL - e)(i) / ( 1 + |e(i)| + opt_scale ), v = 0 ...\n";
    }
    if(!comp_v.comp( c(), 0.0, opt_warning_tol()
             , force_complementarity_error_check ? comp_error_tol() : really_big_error_tol
             , print_all_warnings, out )) test_failed = true;
*/
    
    ///////////////////////////////////
    // e - eU <= 0
    if(out)
      *out
        << sep_line
        << "\nChecking e - eU <= 0 ...\n";
    V_VmV( t_e.get(), *e, *eU );
    if(out && print_vectors)
      *out
        << "\ne - eU =\n" << *t_e;
    Vt_S( t_e.get(), 1.0/(1.0+e_norm_inf) );

    err = max_element(*t_e);
    if(out)
      *out
        << "\nmax(e-eU) = " << err << endl;
    if( force_inequality_error_check )
      handle_error(
        out,err,"max(e-eU)",feas_error_tol(),"feas_error_tol"
        ,feas_warning_tol(),"feas_error_tol",&test_failed
        );
    // ToDo: Update below code!
/*
    if(out)
      *out
        << sep_line
        << "\nChecking muU(i) * (e - eU)(i) = 0  ...\n";
    set_complementarity( *mu, u(), e(), opt_scale, upper, &c );
    if(out && print_vectors)
      *out
        << "\nmuU(i) * (e - eU)(i) / ( 1 + |e(i)| + opt_scale ) =\n" << c();
    if(out) {
      *out
        << "\nComparing:\n"
        << "u(i) = muU(i) * (e - eU)(i) / ( 1 + |e(i)| + opt_scale )\n"
        << "v = 0 ...\n";
    }
    if(!comp_v.comp( c(), 0.0, opt_warning_tol()
             , force_complementarity_error_check ? comp_error_tol() : really_big_error_tol
             , print_all_warnings, out )) test_failed = true;
*/
    
  }

  if( F ) {
    TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Update below code!
/*

    ///////////////////////////////////
    // r = - op(F)*d + eta * f 
    if(out)
      *out
        << sep_line
        << "\nComputing r = - op(F)*d + eta * f ...\n";
    V_StMtV( &r, -1.0, *F, trans_F, *d );
    Vp_StV( &r(), *eta, *f );
    if(out && print_vectors)
      *out
        << "\nr = - op(F)*d + eta * f =\n" << r();

    if(out) {
      *out
        << sep_line
        << "\nChecking r == f:\n"
        << "u = r, v = f ...\n";
    }
    if(!comp_v.comp( r(), *f, opt_warning_tol()
      , force_equality_error_check ? feas_error_tol() : really_big_error_tol
      , print_all_warnings, out )) test_failed = true;

*/
  }

  if(out) {
    *out
      << sep_line;
    if(solution_type != qps_t::SUBOPTIMAL_POINT) {
      if(test_failed) {
        *out
          << "\nDarn it!  At least one of the enforced QP optimality conditions were "
          << "not within the specified error tolerances!\n";
      }
      else {
        *out
          << "\nCongradulations!  All of the enforced QP optimality conditions were "
          << "within the specified error tolerances!\n";
      }
    }
    *out
      << "\n*** End checking QP optimality conditions ***\n";
  }

  return !test_failed;

}

} // end namespace ConstrainedOptPack
