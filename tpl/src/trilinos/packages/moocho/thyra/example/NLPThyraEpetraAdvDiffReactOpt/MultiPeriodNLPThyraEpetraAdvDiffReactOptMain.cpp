/*
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
*/

#include "GLpApp_AdvDiffReactOptModelCreator.hpp"
#include "MoochoPack_MoochoThyraSolver.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_DefaultClusteredSpmdProductVectorSpace.hpp"
#include "Thyra_DefaultMultiPeriodModelEvaluator.hpp"
#include "Thyra_VectorSpaceTester.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_DefaultInverseModelEvaluator.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_OpaqueWrapper.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_Utils.hpp"
#include "Teuchos_DefaultComm.hpp"
#ifdef HAVE_MPI
#  include "Teuchos_DefaultMpiComm.hpp"
#  include "Epetra_MpiComm.h"
#else
#  include "Teuchos_DefaultSerialComm.hpp"
#  include "Epetra_SerialComm.h"
#endif

namespace {

typedef AbstractLinAlgPack::value_type  Scalar;

} // namespace

int main( int argc, char* argv[] )
{

  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcp_implicit_cast;
  using Teuchos::null;
  using Teuchos::RCP;
  using Teuchos::outArg;
  using Teuchos::Array;
  using Teuchos::tuple;
  using Teuchos::ParameterList;
  using Teuchos::OpaqueWrapper;
  using Teuchos::OSTab;
  using Teuchos::CommandLineProcessor;
  using Teuchos::toString;
  using Thyra::VectorBase;
  using Thyra::ProductVectorBase;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Thyra::Ordinal Ordinal;
  using Thyra::ModelEvaluator;
  using MoochoPack::MoochoSolver;
  using MoochoPack::MoochoThyraSolver;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  const int procRank = mpiSession.getRank();
  const int numProcs = mpiSession.getNProc();

  Teuchos::Time timer("");
  
  bool result, success = true;

  RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {
  
    Stratimikos::DefaultLinearSolverBuilder   lowsfCreator;
    GLpApp::AdvDiffReactOptModelCreator     epetraModelCreator;

    // Create the solver object
    MoochoThyraSolver solver;

    //
    // Get options from the command line
    //

    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);

    epetraModelCreator.setupCLP(&clp);
    lowsfCreator.setupCLP(&clp);
    solver.setupCLP(&clp);

    int numProcsPerCluster = -1;
    clp.setOption( "num-procs-per-cluster", &numProcsPerCluster,
      "Number of processes in a cluster (<=0 means only one cluster)." );
    int numPeriodsPerCluster = 1;
    clp.setOption( "num-periods-per-cluster", &numPeriodsPerCluster,
      "Number of periods in a cluster." );
    bool dumpAll = false;
    clp.setOption( "dump-all", "no-dump-all", &dumpAll,
      "Set to true, then a bunch of debugging output will be created for the clustered vector tests." );
    bool skipSolve = false;
    clp.setOption( "skip-solve", "no-skip-solve", &skipSolve,
      "Temporary flag for skip solve for testing." );
    double perturbedParamScaling = 1.0;
    clp.setOption( "p-perturb-scaling", &perturbedParamScaling,
      "Scaling for perturbed paramters from the initial forward solve." );
    bool doMultiPeriod = true;
    clp.setOption( "multi-period", "no-multi-period", &doMultiPeriod,
      "Do a mulit-period solve or not." );
    bool useOuterInverse = true;
    clp.setOption( "use-outer-inverse", "use-inner-inverse", &useOuterInverse,
      "Determines if the outer inverse model will be used or the inner inverse." );
    double periodParamScale = 1.0;
    clp.setOption( "period-param-scale", &periodParamScale,
      "Sets the scaling factor to scale z[i] from one period to the next." );
    bool initialSolveContinuation = false;
    clp.setOption( "init-solve-continuation", "init-solve-all-at-once", &initialSolveContinuation,
      "Determines if the inital solve is done using continuation or all at once." );
    bool useStatelessPeriodModel = false;
    clp.setOption( "use-stateless-period-model", "use-statefull-period-model", &useStatelessPeriodModel,
      "Determines if a stateless or a statefull period model should be used or not." );
    double stateInvError = 1e-8;
    clp.setOption( "state-inv-error", &stateInvError,
      "The error in the l2 norm of the state inverse solution error." );
    double paramInvError = 1e-8;
    clp.setOption( "param-inv-error", &paramInvError,
      "The error in the l2 norm of the parameter inverse solution error." );

    CommandLineProcessor::EParseCommandLineReturn
      parse_return = clp.parse(argc,argv,&std::cerr);

    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL )
      return parse_return;

    lowsfCreator.readParameters(out.get());
    solver.readParameters(out.get());

    *out
      << "\n***"
      << "\n*** NLPThyraEpetraAdvDiffReactOptMain, Global numProcs = "<<numProcs
      << "\n***\n";

    int clusterRank = -1;
    int numClusters = -1;

#ifdef HAVE_MPI

    RCP<OpaqueWrapper<MPI_Comm> >
      intraClusterMpiComm = Teuchos::opaqueWrapper<MPI_Comm>(MPI_COMM_WORLD),
      interClusterMpiComm = Teuchos::null;
    
    {
      if ( numProcsPerCluster <= 0 ) {
        *out
          << "\nnumProcsPerCluster = " << numProcsPerCluster
          << " <= 0: Setting to " << numProcs << "...\n";
        numProcsPerCluster = numProcs;
      }
      *out << "\nCreating communicator for local cluster of "<<numProcsPerCluster<<" processes ...\n";
      numClusters = numProcs/numProcsPerCluster;
      const int remainingProcs = numProcs%numProcsPerCluster;
      TEUCHOS_TEST_FOR_EXCEPTION(
        remainingProcs!=0,std::logic_error
        ,"Error, The number of processes per cluster numProcsPerCluster="<<numProcsPerCluster
        << " does not divide into the global number of processes numProcs="<<numProcs
        << " and instead has remainder="<<remainingProcs<<"!"
        );
      // Determine which cluster this process is part of and what the global
      // process ranges are.
      clusterRank = procRank / numProcsPerCluster; // Integer division!
      *out << "\nclusterRank = " << clusterRank << "\n";
      const int firstClusterProcRank = clusterRank * numProcsPerCluster;
      const int lastClusterProcRank = firstClusterProcRank + numProcsPerCluster - 1;
      *out << "\nclusterProcRange = ["<<firstClusterProcRank<<","<<lastClusterProcRank<<"]\n";
      // Create the communicator for this cluster of processes
      *out << "\nCreating intraClusterMpiComm ...";
      MPI_Comm rawIntraClusterMpiComm = MPI_COMM_NULL;
      MPI_Comm_split(
        MPI_COMM_WORLD        // comm
        ,clusterRank          // color (will all be put in the same output comm)
        ,0                    // key (not important here)
        ,&rawIntraClusterMpiComm // newcomm
        );
      intraClusterMpiComm = Teuchos::opaqueWrapper(rawIntraClusterMpiComm,MPI_Comm_free);
      {
        *out << "\nintraClusterMpiComm:";
        Teuchos::OSTab tab(out);
        int rank, size;
        MPI_Comm_size(*intraClusterMpiComm,&size);
        MPI_Comm_rank(*intraClusterMpiComm,&rank);
        *out << "\nsize="<<size;
        *out << "\nrank="<<rank;
        *out << "\n";
      }
      // Create the communicator for just the root process in each cluster
      *out << "\nCreating interClusterMpiComm ...";
      MPI_Comm rawInterClusterMpiComm = MPI_COMM_NULL;
      MPI_Comm_split(
        MPI_COMM_WORLD                                  // comm
        ,procRank==firstClusterProcRank?0:MPI_UNDEFINED // color
        ,0                                              // key
        ,&rawInterClusterMpiComm                           // newcomm
        );
      if(rawInterClusterMpiComm!=MPI_COMM_NULL)
        interClusterMpiComm = Teuchos::opaqueWrapper(rawInterClusterMpiComm,MPI_Comm_free);
      else
        interClusterMpiComm = Teuchos::opaqueWrapper(rawInterClusterMpiComm);
      {
        *out << "\ninterClusterMpiComm:";
        Teuchos::OSTab tab(out);
        if(*interClusterMpiComm==MPI_COMM_NULL) {
          *out << " NULL\n";
        }
        else {
          int rank, size;
          MPI_Comm_size(*interClusterMpiComm,&size);
          MPI_Comm_rank(*interClusterMpiComm,&rank);
          *out << "\nsize="<<size;
          *out << "\nrank="<<rank;
          *out << "\n";
        }
      }
    }

#endif

    RCP<Epetra_Comm> comm = Teuchos::null;
#ifdef HAVE_MPI
    comm = Teuchos::rcp(new Epetra_MpiComm(*intraClusterMpiComm));
    Teuchos::set_extra_data(intraClusterMpiComm, "mpiComm", outArg(comm));
#else
    comm = Teuchos::rcp(new Epetra_SerialComm());
#endif
    
    //
    // Create the Thyra::ModelEvaluator object
    //
    
    *out << "\nCreate the GLpApp::AdvDiffReactOptModel wrapper object ...\n";
    
    RCP<GLpApp::AdvDiffReactOptModel>
      epetraModel = epetraModelCreator.createModel(comm);

    *out << "\nCreate the Thyra::LinearOpWithSolveFactory object ...\n";

    RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >
      lowsFactory = lowsfCreator.createLinearSolveStrategy("");
    // ToDo: Set the output stream before calling above!
    
    *out << "\nCreate the Thyra::EpetraModelEvaluator wrapper object ...\n";
    
    RCP<Thyra::EpetraModelEvaluator>
      epetraThyraModel = rcp(new Thyra::EpetraModelEvaluator()); // Sets default options!
    epetraThyraModel->setOStream(out);
    epetraThyraModel->initialize(epetraModel,lowsFactory);

    *out
      << "\nnx = " << epetraThyraModel->get_x_space()->dim()
      << "\nnp = " << epetraThyraModel->get_p_space(0)->dim() << "\n";

    //
    // Create the parallel product spaces for x and f
    //

    RCP<const Thyra::ProductVectorSpaceBase<Scalar> >
      x_bar_space, f_bar_space;
    
#ifdef HAVE_MPI

    // For now just build and test these vector spaces if we are not doing a
    // solve!  We have a lot more work to do to the "Clustered" support
    // software before this will work for a solve.
    
    if (skipSolve) {
      
      *out << "\nCreate block parallel vector spaces for multi-period model.x and model.f ...\n";
      RCP<const Teuchos::Comm<Ordinal> >
        intraClusterComm = rcp(new Teuchos::MpiComm<Ordinal>(intraClusterMpiComm)),
        interClusterComm = Teuchos::createMpiComm<Ordinal>(interClusterMpiComm);
      x_bar_space = Teuchos::rcp(
        new Thyra::DefaultClusteredSpmdProductVectorSpace<Scalar>(
          intraClusterComm
          ,0 // clusterRootRank
          ,interClusterComm
          ,1 // numBlocks
          ,tuple<RCP<const Thyra::VectorSpaceBase<Scalar> > >(
            epetraThyraModel->get_x_space()
            ).getRawPtr()
          )
        );
      f_bar_space = Teuchos::rcp(
        new Thyra::DefaultClusteredSpmdProductVectorSpace<Scalar>(
          intraClusterComm
          ,0 // clusterRootRank
          ,interClusterComm
          ,1 // numBlocks
          ,tuple<RCP<const Thyra::VectorSpaceBase<Scalar> > >(
            epetraThyraModel->get_f_space()
            ).getRawPtr()
          )
        );
      
      Thyra::VectorSpaceTester<Scalar> vectorSpaceTester;
      vectorSpaceTester.show_all_tests(true);
      vectorSpaceTester.dump_all(dumpAll);

#ifdef RTOPPACK_SPMD_APPLY_OP_DUMP
      RTOpPack::show_mpi_apply_op_dump = dumpAll;
#endif
#ifdef THYRA_SPMD_VECTOR_BASE_DUMP
      Thyra::SpmdVectorBase<Scalar>::show_dump = dumpAll;
#endif

      *out << "\nTesting the vector space x_bar_space ...\n";
      result = vectorSpaceTester.check(*x_bar_space,OSTab(out).get());
      if(!result) success = false;

      *out << "\nTesting the vector space f_bar_space ...\n";
      result = vectorSpaceTester.check(*f_bar_space,OSTab(out).get());
      if(!result) success = false;
      
      RCP<const VectorBase<Scalar> >
        x0 = epetraThyraModel->getNominalValues().get_x();
      
/* 2008/02/21: rabartl: I have commented this out since it is causing an MPI error?

      double nrm_x0;

      *out << "\nTiming a global reduction across just this cluster: ||x0||_1 = ";
      timer.start(true);
      nrm_x0 = Thyra::norm_1(*x0);
      *out << nrm_x0 << "\n";
      timer.stop();
      *out << "\n    time = " << timer.totalElapsedTime() << " seconds\n";
      
      *out << "\nTiming a global reduction across the entire set of processes: ||x0||_1 = ";
      timer.start(true);
      RTOpPack::ROpNorm1<Scalar> norm_1_op;
      RCP<RTOpPack::ReductTarget> norm_1_targ = norm_1_op.reduct_obj_create();
      const Teuchos::RCP<const Teuchos::Comm<Teuchos_Ordinal> > comm
        = Teuchos::DefaultComm<Ordinal>::getComm();
      const Teuchos::ArrayView<const Teuchos::Ptr<const VectorBase<Scalar> > >
        vecs = Teuchos::tuple(Teuchos::ptrInArg(*x0));
      Teuchos::dyn_cast<const Thyra::SpmdVectorBase<Scalar> >(*x0).applyOpImplWithComm(
        comm.ptr(),
        norm_1_op, vecs, Teuchos::null, norm_1_targ.ptr(),
        0, -1, 0
        );
      nrm_x0 = norm_1_op(*norm_1_targ);
      *out << nrm_x0 << "\n";
      timer.stop();
      *out << "\n    time = " << timer.totalElapsedTime() << " seconds\n";

*/

#ifdef RTOPPACK_SPMD_APPLY_OP_DUMP
      RTOpPack::show_mpi_apply_op_dump = false;
#endif
#ifdef THYRA_SPMD_VECTOR_BASE_DUMP
      Thyra::SpmdVectorBase<Scalar>::show_dump = false;
#endif

    }

#endif // HAVE_MPI
    
    if(skipSolve) {

      if(success)
        *out << "\nEnd Result: TEST PASSED" << std::endl;
      else
        *out << "\nEnd Result: TEST FAILED" << std::endl;

      return ( success ? 0 : 1 );

    }

    const int N = numPeriodsPerCluster;

    Array<RCP<Thyra::ModelEvaluator<Scalar> > >
      inverseThyraModels(N);
    if (useOuterInverse) {
      *out << "\nUsing Thyra::DefaultInverseModelEvaluator for the objective function ...\n";
      if ( useStatelessPeriodModel ) {
        *out << "\nBuilding a single Thyra::DefaultInverseModelEvaluator object where the matching vector will be maintained externally ...\n";
      }
      else {
        *out << "\nBuilding multiple Thyra::DefaultInverseModelEvaluator objects where the matching vector is held internally ...\n";
      }
      RCP<ParameterList> invMEPL = Teuchos::parameterList();
      invMEPL->set( "Observation Multiplier", 1.0 );
      invMEPL->set( "Parameter Multiplier", epetraModel->getDataPool()->getbeta() );
      if ( useStatelessPeriodModel )
        invMEPL->set( "Observation Target as Parameter", true );
      RCP<const Thyra::EpetraLinearOp>
        H = Thyra::epetraLinearOp(epetraModel->getDataPool()->getH(),"H"),
        R = Thyra::epetraLinearOp(epetraModel->getDataPool()->getR(),"R");
      RCP<const Thyra::MultiVectorBase<Scalar> >
        B_bar = Thyra::create_MultiVector(
          epetraModel->get_B_bar(),
          R->spmdRange()
          );
      RCP<const Thyra::LinearOpBase<Scalar> >
        R_bar = Thyra::multiply<Scalar>(Thyra::adjoint<Scalar>(B_bar),R,B_bar);
      for ( int i = 0; i < N; ++i ) {
        if ( ( useStatelessPeriodModel && i==0 ) || !useStatelessPeriodModel ) {
          RCP<Thyra::DefaultInverseModelEvaluator<Scalar> >
            _inverseThyraModel = Thyra::defaultInverseModelEvaluator<Scalar>(
              epetraThyraModel );
          _inverseThyraModel->setParameterList(invMEPL);
          _inverseThyraModel->set_observationMatchWeightingOp(H);
          _inverseThyraModel->set_parameterRegularizationWeightingOp(R_bar);
          inverseThyraModels[i] = _inverseThyraModel;
        }
        else {
#ifdef TEUCHOS_DEBUG
          TEUCHOS_TEST_FOR_EXCEPT( ! ( useStatelessPeriodModel && i > 0 ) );
#endif
          inverseThyraModels[i] = inverseThyraModels[0];
        }
      }
    }
    else {
      *out << "\nUsing built-in inverse objective function ...\n";
      TEUCHOS_TEST_FOR_EXCEPTION(
        N != 1, std::logic_error,
        "Error, you can't have N = "<<N<<" > 1\n"
        "and be using an internal inverse objective!" );
      inverseThyraModels[0] = epetraThyraModel;
    }

    const int p_index = 0;
    const int z_index = 1;
    const int z_p_index = 1; // Ordinal of the reaction rate parameter parameter subvector
    const int z_x_index = 2; // Ordinal of the state matching subvector parameter
    Array<int> z_indexes;
    if (useStatelessPeriodModel)
      z_indexes = tuple<int>(z_p_index, z_x_index);
    else
      z_indexes = tuple<int>(z_p_index);
    Array<Array<RCP<const VectorBase<Scalar> > > > z;
    const int g_index = ( useOuterInverse ? 1 : 0 );
    Array<Scalar> weights;
    RCP<VectorBase<Scalar> >
      z_base = inverseThyraModels[0]->getNominalValues().get_p(z_index)->clone_v();
    *out << "\nz_base =\n" << Teuchos::describe(*z_base,Teuchos::VERB_EXTREME);
    Scalar scale_z_i = 1.0;
    for ( int i = 0; i < N; ++i ) {
      weights.push_back(1.0);
      const RCP<VectorBase<Scalar> > z_i = z_base->clone_v();
      Vt_S( z_i.ptr(), scale_z_i );
      *out << "\nz["<<i<<"] =\n" << Teuchos::describe(*z_i,Teuchos::VERB_EXTREME);
      if ( useStatelessPeriodModel ) {
        z.push_back(
          tuple<RCP<const VectorBase<Scalar> > >(
            z_i,
            null // We will set this again later!
            )
          );
      }
      else {
        z.push_back(
          tuple<RCP<const VectorBase<Scalar> > >(z_i)
          );
      }
      scale_z_i *= periodParamScale;
    }

    RCP<Thyra::ModelEvaluator<Scalar> >
      thyraModel = inverseThyraModels[0];

    if (doMultiPeriod) {
      thyraModel =
        rcp(
          new Thyra::DefaultMultiPeriodModelEvaluator<Scalar>(
            N, inverseThyraModels, z_indexes, z, g_index, weights,
            x_bar_space, f_bar_space
            )
          );
    }
    
    MoochoSolver::ESolutionStatus solution_status;

    //
    *out << "\n***\n*** Solving the initial forward problem\n***\n";
    //

    // Set the solve mode to solve the forward problem
    solver.setSolveMode(MoochoThyraSolver::SOLVE_MODE_FORWARD);

    // Save the solution for model.x and model.p to be used later
    RCP<const VectorBase<Scalar> >
      x_opt, // Will be set below
      x_init, // Will be set below
      p_opt = inverseThyraModels[0]->getNominalValues().get_p(0)->clone_v();

    *out << "\np_opt =\n" << Teuchos::describe(*p_opt,Teuchos::VERB_EXTREME);
    
    if ( initialSolveContinuation ) {
      
      *out << "\nSolving individual period problems one at time using continuation ...\n";

      RCP<ProductVectorBase<Scalar> >
        x_opt_prod = rcp_dynamic_cast<ProductVectorBase<Scalar> >(
         createMember( thyraModel->get_x_space() ), true );

      RCP<const VectorBase<Scalar> > period_x;

      for ( int i = 0; i < N; ++i ) {

        *out << "\nSolving period i = " << i << " using guess from last period ...\n";
      
        // Set the deliminator for the output files!
        solver.getSolver().set_output_context("fwd-init-"+toString(i));

        // Set the period model
        solver.setModel(inverseThyraModels[i]);
        
        // Set the initial guess and the parameter values
        MEB::InArgs<Scalar> initialGuess = inverseThyraModels[i]->createInArgs();
        initialGuess.set_p(z_index,z[i][0]->clone_v());
        initialGuess.set_p(p_index,p_opt->clone_v());
        if ( i == 0 ) {
          // For the first period just use whatever initial guess is built
          // into the model
        }
        else {
          // Set the final solution for x from the last period!
          initialGuess.set_x(period_x);
        }
        solver.setInitialGuess(initialGuess);

        // Solve the period model
        solution_status = solver.solve();
        TEUCHOS_TEST_FOR_EXCEPT( solution_status != MoochoSolver::SOLVE_RETURN_SOLVED );

        // Save the final solution for the next period!
        period_x = solver.getFinalPoint().get_x()->clone_v();
        assign( x_opt_prod->getNonconstVectorBlock(i).ptr(), *period_x );
        if ( useStatelessPeriodModel )
          z[i][1] = period_x->clone_v(); // This is our matching vector!
        
      }

      x_opt = x_opt_prod;
      x_init = x_opt->clone_v();

      if ( useStatelessPeriodModel ) {
        rcp_dynamic_cast<Thyra::DefaultMultiPeriodModelEvaluator<Scalar> >(
          thyraModel
          )->reset_z(z);
      }

    }
    else {

      *out << "\nSolving all periods simultaniously ...\n";
      
      // Set the deliminator for the output files!
      solver.getSolver().set_output_context("fwd-init");
      
      // Set the model
      solver.setModel(thyraModel);
      
      // Set the initial guess from files (if specified on commandline)
      solver.readInitialGuess(out.get());
      
      // Solve the initial forward problem
      solution_status = solver.solve();
      TEUCHOS_TEST_FOR_EXCEPT( solution_status != MoochoSolver::SOLVE_RETURN_SOLVED );

      // Save the solution for model.x and model.p to be used later
      x_opt = solver.getFinalPoint().get_x()->clone_v();
      x_init = solver.getFinalPoint().get_x()->clone_v();
      
    }
    
    //
    *out << "\n***\n*** Solving the perturbed forward problem\n***\n";
    //
    
    // Set the deliminator for the output files!
    solver.getSolver().set_output_context("fwd");
    
    // Set the solve mode to solve the forward problem
    solver.setSolveMode(MoochoThyraSolver::SOLVE_MODE_FORWARD);
    
    // Set the model
    solver.setModel(thyraModel);
    
    // Set the initial guess and the perturbed parameters
    RCP<VectorBase<Scalar> >
      p_init = p_opt->clone_v();
    {
      MEB::InArgs<Scalar> initialGuess = thyraModel->createInArgs();
      initialGuess.setArgs(thyraModel->getNominalValues());
      initialGuess.set_x(x_init);
      Thyra::Vt_S(p_init.ptr(), perturbedParamScaling);
      initialGuess.set_p(0,p_init);
      //*out << "\nInitial Guess:\n" << Teuchos::describe(initialGuess,Teuchos::VERB_EXTREME);
      solver.setInitialGuess(initialGuess);
    }

    // Solve the perturbed forward problem
    solution_status = solver.solve();
    
    // Save the solution for model.x and model.p to be used later
    x_init = solver.getFinalPoint().get_x()->clone_v();
    p_init = solver.getFinalPoint().get_p(0)->clone_v();

    *out
      << "\nrelVectorErr(x_perturb,x_opt) = " << Thyra::relVectorErr(*x_init,*x_opt)
      << "\nrelVectorErr(p_perturb,p_opt) = " << Thyra::relVectorErr(*p_init,*p_opt)
      << "\n";
    
    //
    *out << "\n***\n*** Solving the perturbed inverse problem\n***\n";
    //

    // Set the deliminator for the output files!
    solver.getSolver().set_output_context("inv");

    //TEUCHOS_TEST_FOR_EXCEPT("ToDo: We need to use the DefaultInverseModelEvaluator to set the matching vector correctly!");

    // Set the matching vector
    if ( N > 1 ) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        !useOuterInverse, std::logic_error,
        "Error, if N > 1, you have to use the outer inverse objective function\n"
        "since each target vector will be different!" );
      RCP<const ProductVectorBase<Scalar> >
        x_opt_prod = rcp_dynamic_cast<const ProductVectorBase<Scalar> >(
          rcp_implicit_cast<const VectorBase<Scalar> >(x_opt), true
          ); // This cast can *not* fail!
      for ( int i = 0; i < N; ++i ) {
        rcp_dynamic_cast<Thyra::DefaultInverseModelEvaluator<Scalar> >(
          inverseThyraModels[i], true
          )->set_observationTarget(x_opt_prod->getVectorBlock(i));
      }
    }
    else if ( 1 == N ) {
      RCP<const ProductVectorBase<Scalar> >
        x_opt_prod = rcp_dynamic_cast<const ProductVectorBase<Scalar> >(
          rcp_implicit_cast<const VectorBase<Scalar> >(x_opt)
          ); // This cast can fail!
      RCP<const VectorBase<Scalar> >
        x_opt_i =
        ( !is_null(x_opt_prod)
          ? x_opt_prod->getVectorBlock(0)
          : rcp_implicit_cast<const VectorBase<Scalar> >(x_opt)
          );
      if (useOuterInverse) {
        rcp_dynamic_cast<Thyra::DefaultInverseModelEvaluator<Scalar> >(
          inverseThyraModels[0], true
          )->set_observationTarget(x_opt_i);
      }
      else {
        epetraModel->set_q(
          Thyra::get_Epetra_Vector(*epetraModel->get_x_map(), x_opt_i)
          );
      }
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPT("Error, should not get here!");
    }
    
    // Set the solve mode to solve the inverse problem
    solver.setSolveMode(MoochoThyraSolver::SOLVE_MODE_OPTIMIZE);
   
    // Set the model
    solver.setModel(thyraModel);

    // Set the initial guess for model.x and model.p
    {
      MEB::InArgs<Scalar> initialGuess = thyraModel->createInArgs();
      initialGuess.setArgs(thyraModel->getNominalValues());
      initialGuess.set_x(x_init);
      initialGuess.set_p(0,p_init);
      //*out << "\nInitial Guess:\n" << Teuchos::describe(initialGuess,Teuchos::VERB_EXTREME);
      solver.setInitialGuess(initialGuess);
    }
    
    // Solve the inverse problem
    solution_status = solver.solve();
    TEUCHOS_TEST_FOR_EXCEPT( solution_status != MoochoSolver::SOLVE_RETURN_SOLVED );

    //
    *out << "\n***\n*** Testing the error in the inversion\n***\n";
    //
    
    // Get the inverted for solution and compare it to the optimal solution

    RCP<const VectorBase<Scalar> >
      x_inv = solver.getFinalPoint().get_x(),
      p_inv = solver.getFinalPoint().get_p(0);

    *out << "\np_opt =\n" << Teuchos::describe(*p_opt,Teuchos::VERB_EXTREME);
    *out << "\np_inv =\n" << Teuchos::describe(*p_inv,Teuchos::VERB_EXTREME);

    const Scalar
      x_err = Thyra::relVectorErr( *x_inv, *x_opt ),
      p_err = Thyra::relVectorErr( *p_inv, *p_opt );

    const bool
      x_test_passed = ( x_err <= stateInvError ),
      p_test_passed = ( p_err <= paramInvError );

    *out
      << "\nrelVectorErr(x_inv,x_opt) = " << x_err << " <= " << stateInvError
      << " : " << Thyra::passfail(x_test_passed)
      << "\nrelVectorErr(p_inv,p_opt) = " << p_err << " <= " << paramInvError
      << " : " << Thyra::passfail(p_test_passed)
      << "\n";
    
    // Write the final solution
    solver.writeFinalSolution(out.get());
    
    // Write the final parameters to file
    lowsfCreator.writeParamsFile(*lowsFactory);
    solver.writeParamsFile();
    
    TEUCHOS_TEST_FOR_EXCEPT(
      solution_status != MoochoSolver::SOLVE_RETURN_SOLVED
      );
    
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,success)

  if(success)
    *out << "\nEnd Result: TEST PASSED" << std::endl;
  else
    *out << "\nEnd Result: TEST FAILED" << std::endl;
  
  return ( success ? 0 : 1 );
  
}
