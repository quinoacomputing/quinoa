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

#include <assert.h>

#include <fstream>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "NLPInterfacePack_ExampleNLPDirectRun.hpp"
#include "AbstractLinAlgPack_exampleNLPDiagSetup.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_oblackholestream.hpp"

int main( int argc, char* argv[] ) {

  using std::endl;
  using std::setw;
  namespace mmp = MemMngPack;
  using Teuchos::RCP;
  typedef AbstractLinAlgPack::size_type size_type;
  typedef AbstractLinAlgPack::value_type value_type;
  namespace NLPIP = NLPInterfacePack;
  namespace rsqp = MoochoPack;
  using rsqp::MoochoSolver;

  using AbstractLinAlgPack::VectorSpace;

  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  int err = 0;

/*
  // Print out the input arguments
  printf("argc = %d\n",argc);
  {for( int i = 0; i < argc; ++i) {
    printf("argv[%d] = %s\n",i,argv[i]);
  }}
*/
  // Get an idea of what processors we have.
  MPI_Init(&argc,&argv);
  int num_proc, proc_rank;
  MPI_Comm_size( MPI_COMM_WORLD, &num_proc );
  MPI_Comm_rank( MPI_COMM_WORLD, &proc_rank );
/*
  // Print out the input arguments
  printf("\nproc_rank = %d\n",proc_rank);
  printf("argc = %d\n",argc);
  {for( int i = 0; i < argc; ++i) {
    printf("argv[%d] = %s\n",i,argv[i]);
  }}
*/

  // Define program return values
  const int
    PROG_SUCCESS				=  0,
    PROG_NLP_TEST_ERR			= -1,
    PROG_EXCEPTION				= -2,
    PROG_MAX_ITER_EXEEDED		= -3,
    PROG_MAX_TIME_EXEEDED		= -4;

  int prog_return = PROG_SUCCESS;

  // Set the output stream
  std::ostream &out  = std::cout;
  std::ostream &eout = std::cerr;
  Teuchos::oblackholestream  blackhole;

  try {
  
    //
    // Initialize stuff
    //

    size_type n;
    value_type xo;
    bool has_bounds;
    bool dep_bounded;
    
    VectorSpace::space_ptr_t    vec_space;
    const int err = AbstractLinAlgPack::exampleNLPDiagSetup(argc,argv,MPI_COMM_WORLD,&vec_space,&n,&xo,&has_bounds,&dep_bounded);
    if(err) return err;
    
    // Create and test the NLP using this vector space object
    const MoochoSolver::ESolutionStatus
      solve_return = NLPIP::ExampleNLPDirectRun(
        *vec_space, xo, has_bounds, dep_bounded
        ,proc_rank == 0 ? &out  : &blackhole  // console_out
        ,proc_rank == 0 ? &eout : &blackhole  // error_out
        ,proc_rank == 0 ? false : true        // throw_solve_exception
        ,proc_rank == 0 ? NULL  : &blackhole  // algo_out
        ,proc_rank == 0 ? NULL  : &blackhole  // summary_out
        ,proc_rank == 0 ? NULL  : &blackhole  // journal_out
        );
    
    switch(solve_return) {
      case MoochoSolver::SOLVE_RETURN_SOLVED:
        prog_return = PROG_SUCCESS;
        break;
      case MoochoSolver::SOLVE_RETURN_MAX_ITER:
        prog_return = PROG_MAX_ITER_EXEEDED;
        break;
      case MoochoSolver::SOLVE_RETURN_MAX_RUN_TIME:
        prog_return = PROG_MAX_TIME_EXEEDED;
        break;
      case MoochoSolver::SOLVE_RETURN_NLP_TEST_FAILED:
        prog_return = PROG_NLP_TEST_ERR;
        break;
      case MoochoSolver::SOLVE_RETURN_EXCEPTION:
        prog_return = PROG_EXCEPTION;
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPT(true);
    }
    
  }	// end try
  catch(const std::exception& excpt) {
    eout << "\nCaught a std::exception on process " << proc_rank<< ": " << excpt.what() << endl;
    prog_return = PROG_EXCEPTION;
  }
  catch(...) {
    eout << "\nCaught an unknown exception on process " << proc_rank<< "\n";
    prog_return = PROG_EXCEPTION;
  }

   MPI_Finalize();

  return prog_return;
}
