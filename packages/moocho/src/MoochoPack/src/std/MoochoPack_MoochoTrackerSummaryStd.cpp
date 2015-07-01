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

#include <assert.h>

#include <iomanip>

#include "MoochoPack_MoochoTrackerSummaryStd.hpp"
#include "MoochoPack_NLPAlgoState.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "NLPInterfacePack_NLPFirstOrder.hpp"
#include "AbstractLinAlgPack_Vector.hpp"
#include "AbstractLinAlgPack_MatrixSymOp.hpp"
#include "Teuchos_dyn_cast.hpp"

using std::endl;
using std::setw;

namespace MoochoPack {

MoochoTrackerSummaryStd::MoochoTrackerSummaryStd(
  const ostream_ptr_t      &o
  ,const ostream_ptr_t     &journal_out
  ,EOptError               opt_error
  )
  :AlgorithmTracker(journal_out)
  ,o_(o)
  ,opt_error_(opt_error)
  ,num_total_qp_iter_(0)
{}	

void MoochoTrackerSummaryStd::set_output_stream(const ostream_ptr_t& o)
{	
  o_ = o;
}

const MoochoTrackerSummaryStd::ostream_ptr_t&
MoochoTrackerSummaryStd::get_output_stream() const
{
  return o_;
}

void MoochoTrackerSummaryStd::output_iteration(const Algorithm& algo) const
{

  const NLPAlgo            &_algo  = rsqp_algo(algo);
  const NLPAlgoState           &s      =_algo.rsqp_state();
  const NLP                 &nlp    = _algo.nlp(); 

  const size_type
    m = nlp.m();

  std::ostream& o = this->o();
  
  int w = 15;
  int prec = 6;
  o.precision(prec);

  // Output the table's header for the first iteration
  if(s.k() == 0) {
    print_header(s);
  }

  // ///////////////////////////////
  // Output a row for the iteration
  
  // Get active set and QP solver statistics.
  const ActSetStats		*act_stats =
    ( act_set_stats_.exists_in(s) && act_set_stats_(s).updated_k(0)
      ? &act_set_stats_(s).get_k(0)
      : NULL );
  const QPSolverStats		*qp_stats =
    ( qp_solver_stats_.exists_in(s) && qp_solver_stats_(s).updated_k(0)
      ? &qp_solver_stats_(s).get_k(0)
      : NULL );
  const QuasiNewtonStats	*quasi_newt_stats =
    ( quasi_newton_stats_.exists_in(s) && quasi_newton_stats_(s).updated_k(0)
      ? &quasi_newton_stats_(s).get_k(0)
      : NULL );

  // Get the norms of Ypy and Zpz
  value_type norm_2_Ypy = -1.0, norm_2_Zpz = -1.0;
  bool Ypy_exists, Zpz_exists;
  if( m && ( Ypy_exists = s.Ypy().updated_k(0) ) )
    norm_2_Ypy = s.Ypy().get_k(0).norm_2();
  if( Zpz_exists = s.Zpz().updated_k(0) )
    norm_2_Zpz = s.Zpz().get_k(0).norm_2();

  o	<< std::right
    << setw(5) << s.k();

  if( s.f().updated_k(0) )
    o << setw(w) << s.f().get_k(0);
  else
    o << setw(w) << "-";

  if( s.Gf().updated_k(0) )
    o << setw(w) << s.Gf().get_k(0).norm_inf();
  else
    o << setw(w) << "-";

  if( m && s.c().updated_k(0) )
    o << setw(w)
      << s.c().get_k(0).norm_inf();
  else
    o << setw(w) << "-";

  {
    const IterQuantityAccess<VectorMutable>
      &rGL_GL = ( opt_error_ == OPT_ERROR_REDUCED_GRADIENT_LAGR
              ? s.rGL() : s.GL()  );
    if( rGL_GL.updated_k(0) )
      o << setw(w) << rGL_GL.get_k(0).norm_inf();
    else
      o << setw(w) << "-";
  }

  if( quasi_newt_stats ) {
    o << setw(w);
    switch( quasi_newt_stats->updated() ) {
      case QuasiNewtonStats::UNKNOWN:
        o << "-";
        break;
      case QuasiNewtonStats:: REINITIALIZED:
        o << "initialized";
        break;
      case QuasiNewtonStats::DAMPENED_UPDATED:
        o << "damp.updated";
        break;
      case QuasiNewtonStats::UPDATED:
        o << "updated";
        break;
      case QuasiNewtonStats::SKIPED:
        o << "skiped";
        break;
      case QuasiNewtonStats::INDEF_SKIPED:
        o << "indef skiped";
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPT(true);
    }
  }
  else {
    o << setw(w) << "-";
  }

  if( act_stats ) {
    o << setw(7) << act_stats->num_active();
    // don't know num_add and num_drops on first iteration.
    if( act_stats->num_adds() == ActSetStats::NOT_KNOWN ) { 
      o	<< setw(7) << "-";
    }
    else {		
      o	<< setw(7) << act_stats->num_adds();
    }
    if( act_stats->num_drops() == ActSetStats::NOT_KNOWN ) {
      o	<< setw(7) << "-";
    }
    else {		
      o	<< setw(7) << act_stats->num_drops();
    }
  }
  else if( s.nu().updated_k(0) ) {
    o	<< setw(7) << s.nu().get_k(0).nz()
      << setw(7) << "-"
      << setw(7) << "-";
  }
  else {
    o	<< setw(7) << "-"
      << setw(7) << "-"
      << setw(7) << "-";
  }

  if( qp_stats ) {
    o	<< setw(7) << qp_stats->num_qp_iter()
      << setw(3) << ( qp_stats->warm_start() ? 'w' : 'c')
      << setw(2) << ( qp_stats->infeasible_qp() ? 'i' : 'f');
    num_total_qp_iter_ += qp_stats->num_qp_iter();
  }
  else {
  o	<< setw(7) << "-"
    << setw(3) << "-"
    << setw(2) << "-";
  }

  if(m && Ypy_exists)
    o	<< setw(w) << norm_2_Ypy;
  else
    o << setw(w) << "-";

  if(Zpz_exists)
    o	<< setw(w) << norm_2_Zpz;
  else
    o << setw(w) << "-";

  if( s.d().updated_k(0) )
    o << setw(w)
      << s.d().get_k(0).norm_inf();
  else
    o << setw(w) << "-";

  if( s.alpha().updated_k(0) )
    o << setw(w) << s.alpha().get_k(0);
  else
    o << setw(w) << "-";

  o << std::endl;
}

void MoochoTrackerSummaryStd::output_final(const Algorithm& algo
  , EAlgoReturn algo_return) const
{
  using Teuchos::dyn_cast;

  const NLPAlgo            &_algo  = rsqp_algo(algo);
  const NLPAlgoState           &s      =_algo.rsqp_state();
  const NLPObjGrad      &nlp    = dyn_cast<const NLPObjGrad>(_algo.nlp()); 
  const NLPFirstOrder  *nlp_foi = dynamic_cast<const NLPFirstOrder*>(&nlp); 

  const size_type
    m = nlp.m();

  std::ostream& o = this->o();

  int w = 15;
  int prec = 6;
  o.precision(prec);

  // Get active set, QP solver and quasi-newton statistics.
  const ActSetStats		*act_stats =
    ( act_set_stats_.exists_in(s) && act_set_stats_(s).updated_k(0)
      ? &act_set_stats_(s).get_k(0)
      : NULL );
  const QPSolverStats		*qp_stats =
    ( qp_solver_stats_.exists_in(s) && qp_solver_stats_(s).updated_k(0)
      ? &qp_solver_stats_(s).get_k(0)
      : NULL );
  const QuasiNewtonStats	*quasi_newt_stats =
    ( quasi_newton_stats_.exists_in(s) && quasi_newton_stats_(s).updated_k(0)
      ? &quasi_newton_stats_(s).get_k(0)
      : NULL );

  // Output the table's header for the first iteration
  if(s.k() == 0) {
    print_header(s);
  }
  else {
    o	<< " ----"
      << "   ------------"
      << "   ------------"
      << "   ------------"
      << "   ------------"
      << "   ------------"
      << " ------"
      << " ------"
      << " ------"
      << " ------"
      << " ----\n";
  }

  o	<< std::right
    << setw(5) << s.k();

  if( s.f().updated_k(0) )
    o << setw(w) << s.f().get_k(0);
  else
    o << setw(w) << "-";

  if( s.Gf().updated_k(0) )
    o << setw(w) << s.Gf().get_k(0).norm_inf();
  else
    o << setw(w) << "-";

  if( m && s.c().updated_k(0) )
    o << setw(w)
      << s.c().get_k(0).norm_inf();
  else
    o << setw(w) << "-";

  {
    const IterQuantityAccess<VectorMutable>
      &rGL_GL = ( opt_error_ == OPT_ERROR_REDUCED_GRADIENT_LAGR
              ? s.rGL() : s.GL()  );
    if( rGL_GL.updated_k(0) )
      o << setw(w) << rGL_GL.get_k(0).norm_inf();
    else
      o << setw(w) << "-";
  }

  o << setw(w);
  if( quasi_newt_stats ) {
    switch( quasi_newt_stats->updated() ) {
      case QuasiNewtonStats::UNKNOWN:
        o << "-";
        break;
      case QuasiNewtonStats:: REINITIALIZED:
        o << "initialized";
        break;
      case QuasiNewtonStats::DAMPENED_UPDATED:
        o << "damp.updated";
        break;
      case QuasiNewtonStats::UPDATED:
        o << "updated";
        break;
      case QuasiNewtonStats::SKIPED:
        o << "skiped";
        break;
      case QuasiNewtonStats::INDEF_SKIPED:
        o << "indef skiped";
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPT(true);
    }
  }
  else {
    o	<< setw(w) << "-";;
  }

  if( act_stats ) {
    o	<< setw(7) << act_stats->num_active()
      << setw(7) << act_stats->num_adds()
      << setw(7) << act_stats->num_drops();
  }
  else if( s.nu().updated_k(0) ) {
    o	<< setw(7) << s.nu().get_k(0).nz()
      << setw(7) << "-"
      << setw(7) << "-";
  }
  else {
    o	<< setw(7) << "-"
      << setw(7) << "-"
      << setw(7) << "-";
  }
  
  if( qp_stats ) {
    o << setw(7) << qp_stats->num_qp_iter()
      << setw(3) << ( qp_stats->warm_start() ? 'w' : 'c')
      << setw(2) << ( qp_stats->infeasible_qp() ? 'i' : 'f');
    num_total_qp_iter_ += qp_stats->num_qp_iter();
  }
  else {
    o << setw(7) << "-"
      << setw(3) << "-"
      << setw(2) << "-";
  }
  
  if(m && s.Ypy().updated_k(0))
    o	<< setw(w) 
      << s.Ypy().get_k(0).norm_2();
  else
    o << setw(w) << "-";
  
  if(s.Zpz().updated_k(0))
    o	<< setw(w)
      << s.Zpz().get_k(0).norm_2();
  else
    o << setw(w) << "-";
  
  if( s.d().updated_k(0) )
    o << setw(w)
      << s.d().get_k(0).norm_inf();
  else
    o << setw(w) << "-";

  o << endl;
  
  o	<< "\nNumber of function evaluations:\n"
    <<     "-------------------------------\n"
    << "f(x)  : " << nlp.num_f_evals() << endl
    << "c(x)  : " << ( m ? nlp.num_c_evals() : 0 ) << endl
    << "Gf(x) : " << nlp.num_Gf_evals() << endl
    << "Gc(x) : ";
  if(m){
    if( nlp_foi )
      o << nlp_foi->num_Gc_evals();
    else
      o << "?";
  }
  else {
    o << 0;
  }
  o << endl;

}

void MoochoTrackerSummaryStd::print_header(const NLPAlgoState &s) const
{
  // Reset the count of total QP iterations
  num_total_qp_iter_ = 0;

  std::ostream& o = this->o();

  NLPAlgoState::space_c_ptr_t
    space_c = s.get_space_c();

  o	<< "\n\n********************************\n"
    << "*** Start of rSQP Iterations ***\n"
    << "n = " << s.space_x().dim()
    << ", m = " << ( space_c.get() ? space_c->dim() : 0 )
    << ", nz = ";
  try {
    if(space_c.get()) {
      if( s.Gc().updated_k(0) )
        o	<< s.Gc().get_k(0).nz() << endl;
      else
        o	<< "?\n";
    }
    else {
      o	<< 0 << endl;
    }
  }
  catch( const AlgorithmState::DoesNotExist& ) {
      o	<< "?\n";
  }
  o
    << "\n k   "
    << "   f           "
    << "   ||Gf||inf   "
    << "   ||c||inf    ";
  switch(opt_error_) {
    case OPT_ERROR_REDUCED_GRADIENT_LAGR:
      o	<< "   ||rGL||inf  ";
      break;
    case OPT_ERROR_GRADIENT_LAGR:
      o	<< "   ||GL||inf  ";
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);
  }
  o 	<< "   quasi-Newton"
    << " #act  "
    << " #adds "
    << " #drops"
    << " #qpitr"
    << " wcfi  "
    << "   ||Ypy||2    "
    << "   ||Zpz||2    "
    << "   ||d||inf    "
    << "   alpha\n"
    << " ----"
    << "   ------------"
    << "   ------------"
    << "   ------------"
    << "   ------------"
    << "   ------------"
    << " ------"
    << " ------"
    << " ------"
    << " ------"
    << " ----"
    << "   ------------"
    << "   ------------"
    << "   ------------"
    << "   ------------\n";
}

} // end namespace MoochoPack
