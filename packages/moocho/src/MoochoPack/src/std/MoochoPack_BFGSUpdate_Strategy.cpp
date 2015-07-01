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

#include <limits>

#include "MoochoPack_BFGSUpdate_Strategy.hpp"
#include "MoochoPack_Exceptions.hpp"
#include "AbstractLinAlgPack_TestMatrixSymSecant.hpp"
#include "AbstractLinAlgPack_MatrixSymSecant.hpp"
#include "AbstractLinAlgPack_MatrixSymOp.hpp"
#include "AbstractLinAlgPack_MatrixOpOut.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_Assert.hpp"

namespace MoochoPack {

BFGSUpdate_Strategy::BFGSUpdate_Strategy(
  bool               rescale_init_identity
  ,bool              use_dampening
  ,ESecantTesting    secant_testing
  ,value_type        secant_warning_tol
  ,value_type        secant_error_tol
  )
  :rescale_init_identity_(rescale_init_identity)
  ,use_dampening_(use_dampening)
  ,secant_testing_(secant_testing)
  ,secant_warning_tol_(secant_warning_tol)
  ,secant_error_tol_(secant_error_tol)
{}

void BFGSUpdate_Strategy::perform_update(
  VectorMutable            *s_bfgs
  ,VectorMutable           *y_bfgs
  ,bool                    first_update
  ,std::ostream            &out
  ,EJournalOutputLevel     olevel
  ,bool                    check_results
  ,MatrixSymOp             *B
  ,QuasiNewtonStats        *quasi_newton_stats 
  )
{
  namespace rcp = MemMngPack;
  using Teuchos::dyn_cast;
  using AbstractLinAlgPack::dot;
  using AbstractLinAlgPack::Vt_S;
  using AbstractLinAlgPack::Vp_StV;
  using LinAlgOpPack::Vp_V;
  using LinAlgOpPack::V_StV;
  using LinAlgOpPack::V_VpV;
  using LinAlgOpPack::V_VmV;
  using LinAlgOpPack::V_MtV;

  const value_type
    NOT_CALCULATED = std::numeric_limits<value_type>::max();
  value_type
    sTy = NOT_CALCULATED,
    yTy = NOT_CALCULATED;
  bool used_dampening = false;

  // Make sure the update is defined!
  if( sTy == NOT_CALCULATED )
    sTy = dot( *s_bfgs, *y_bfgs );
  if( sTy == 0 ) {
    if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
      out << "\n(y'*s) == 0, skipping the update ...\n";
    }
    quasi_newton_stats->set_updated_stats( QuasiNewtonStats::INDEF_SKIPED );
    return;
  }

  MatrixSymSecant
      &B_updatable = dyn_cast<MatrixSymSecant>(*B);

  // /////////////////////////////////////////////////////////////
  // Consider rescaling the initial identity hessian before
  // the update.
  // 
  // This was taken from Nocedal & Wright, 1999, p. 200 
  // 
  // Bo = (y'*y)/(y'*s) * I
  // 
  if( first_update && rescale_init_identity() ) {
    if( sTy == NOT_CALCULATED )
      sTy = dot( *s_bfgs, *y_bfgs );
    if( yTy == NOT_CALCULATED )
      yTy = dot( *y_bfgs, *y_bfgs );
    const value_type
      Iscale = yTy/sTy;
    const value_type
      Iscale_too_small = 1e-5;  // ToDo: Make this adjustable
    if( Iscale >= Iscale_too_small ) {	
      if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
        out << "\nRescaling the initial identity matrix before the update as:\n"
          << "Iscale = (y'*y)/(y'*s) = ("<<yTy<<")/("<<sTy<<") = "<<Iscale<<"\n"
          << "B =  Iscale * eye(n-r) ...\n";
      }
      B_updatable.init_identity( y_bfgs->space(), Iscale );
      if( (int)olevel >= (int)PRINT_ITERATION_QUANTITIES ) {
        out << "\nB after rescaling = \n" << *B;
      }
    }
    else {
      if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
        out << "\nSkipping rescaling of the initial identity matrix since:\n"
          << "Iscale = (y'*y)/(y'*s) = ("<<yTy<<")/("<<sTy<<") = "<<Iscale
          << " < Iscale_too_small = " << Iscale_too_small << std::endl;
      }
    }
  }

  // ////////////////////////////////////////////////////
  // Modify the s_bfgs and y_bfgs vectors for dampening?
  VectorSpace::vec_mut_ptr_t
    Bs = Teuchos::null;
  if( use_dampening() ) {
    if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS )
    {
      out
        << "\nConsidering the dampened update ...\n";
    }
    // Bs = Bm1 * s_bfgs
    Bs = y_bfgs->space().create_member();
    V_MtV( Bs.get(), *B, BLAS_Cpp::no_trans, *s_bfgs );
    // sTy = s_bfgs' * y_bfgs
    if( sTy == NOT_CALCULATED )
      sTy = dot( *s_bfgs, *y_bfgs );
    // sTBs = s_bfgs' * Bs
    const value_type sTBs = dot( *s_bfgs, *Bs );
    // Select dampening parameter theta
    const value_type
      theta = ( sTy >= 0.2 * sTBs )
      ? 1.0
      : (0.8 * sTBs ) / ( sTBs - sTy );
    if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS )
    {
      out
        << "\ns_bfgs'*y_bfgs = " << sTy
        << ( theta == 1.0 ? " >= " : " < " )
        << " s_bfgs' * B * s_bfgs = " << sTBs << std::endl;
    }
    if( theta == 1.0 ) {
      if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS )
      {
        out
          << "Perform the undamped update ...\n";
      }
    }
    else {
      // y_bfgs = theta*y_bfgs + (1-theta)*B*s_bfgs
      Vt_S( y_bfgs, theta );
      Vp_StV( y_bfgs, (1.0-theta), *Bs );

      if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS )
      {
        out
          << "Dampen the update ...\n"
          << "theta = " << theta << std::endl
          << "y_bfgs = theta*y_bfgs + (1-theta)*B*s_bfgs ...\n"
          << "||y_bfgs||inf = " << y_bfgs->norm_inf() << std::endl;
      }

      if( (int)olevel >= (int)PRINT_VECTORS )
      {
        out
          << "y_bfgs =\n" << *y_bfgs;
      }

      used_dampening = true;
    }
  }

  // Perform the update if it is defined (s_bfgs' * y_bfgs > 0.0)
        
  VectorSpace::vec_mut_ptr_t
    s_bfgs_save = Teuchos::null,
    y_bfgs_save = Teuchos::null;
  if( check_results ) {
    // Save s and y since they may be overwritten in the update.
    s_bfgs_save = s_bfgs->clone();
    y_bfgs_save = y_bfgs->clone();
  }
  try {
    B_updatable.secant_update(
      s_bfgs
      ,y_bfgs
      ,Bs.get()
      );
  }
  catch( const MatrixSymSecant::UpdateFailedException& excpt ) {
    if( (int)olevel >= (int)PRINT_BASIC_ALGORITHM_INFO ) {
      out << "\nThe factorization of B failed in a major way!  Rethrow!\n";
    }
    throw;
  }					
  catch( const MatrixSymSecant::UpdateSkippedException& excpt ) {
    if( (int)olevel >= (int)PRINT_BASIC_ALGORITHM_INFO ) {
      out << excpt.what() << std::endl
        << "\nSkipping BFGS update.  B = B ...\n";
    }
    quasi_newton_stats->set_updated_stats(
      QuasiNewtonStats::INDEF_SKIPED );
    return;
  }					
    
  if( (int)olevel >= (int)PRINT_ITERATION_QUANTITIES ) {
    out << "\nB after the BFGS update = \n" << *B;
  }
  
  if( ( check_results && secant_testing() == SECANT_TEST_DEFAULT )
    || secant_testing() == SECANT_TEST_ALWAYS )
  {
    const bool result =
      AbstractLinAlgPack::TestMatrixSymSecant(
        *B, *s_bfgs_save, *y_bfgs_save
        , secant_warning_tol(), secant_error_tol()
        , (int)olevel >= (int)PRINT_VECTORS
        , (int)olevel >  (int)PRINT_NOTHING ? &out : NULL
        , (int)olevel >= (int)PRINT_ALGORITHM_STEPS
        );
    if( !result ) {
      const char
        msg[] =	"Error, the secant property for the BFGS update failed\n"
        "Stopping the algorithm ...\n";
      out << msg;
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, TestFailed
        ," BFGSUpdate_Strategy::perform_update(...) : " << msg );
    }
  }
  
  quasi_newton_stats->set_updated_stats(
    used_dampening
    ? QuasiNewtonStats::DAMPENED_UPDATED
    : QuasiNewtonStats::UPDATED );
}

void BFGSUpdate_Strategy::print_step( std::ostream& out, const std::string& L ) const
{
  out
    << L << "if use_dampening == true then\n"
    << L << "    if s'*y >= 0.2 * s'*B*s then\n"
    << L << "        theta = 1.0\n"
    << L << "    else\n"
    << L << "        theta = 0.8*s'*B*s / (s'*B*s - s'*y)\n"
    << L << "    end\n"
    << L << "    y = theta*y + (1-theta)*B*s\n"
    << L << "end\n"
    << L << "if first_update && rescale_init_identity and y'*s is sufficently positive then\n"
    << L << "    B = |(y'*y)/(y'*s)| * eye(size(s))\n"
    << L << "end\n"
    << L << "if s'*y is sufficently positive then\n"
    << L << "    *** Peform BFGS update\n"
    << L << "    (B, s, y ) -> B (through MatrixSymSecantUpdate interface)\n"
    << L << "    if ( check_results && secant_testing == SECANT_TEST_DEFAULT )\n"
    << L << "    or ( secant_testing == SECANT_TEST_ALWAYS ) then\n"
    << L << "        if B*s != y then\n"
    << L << "            *** The secant condition does not check out\n"
    << L << "            Throw TestFailed!\n"
    << L << "        end\n"
    << L << "    end\n"
    << L << "end\n"
    ;
}

}  // end namespace MoochoPack
