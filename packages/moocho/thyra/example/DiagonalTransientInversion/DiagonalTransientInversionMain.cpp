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


#include "EpetraExt_DiagonalTransientModel.hpp"
#include "MoochoPack_MoochoThyraSolver.hpp"
#include "Rythmos_BackwardEulerStepper.hpp"
#include "Rythmos_ImplicitBDFStepper.hpp"
#include "Rythmos_ForwardSensitivityStepper.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"
#include "Rythmos_DefaultIntegrator.hpp"
#include "Rythmos_SimpleIntegrationControlStrategy.hpp"
#include "Rythmos_StepperAsModelEvaluator.hpp"
#include "Rythmos_ForwardSensitivityIntegratorAsModelEvaluator.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_DefaultSpmdMultiVectorFileIO.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_DefaultInverseModelEvaluator.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_DefaultFiniteDifferenceModelEvaluator.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_VerbosityLevelCommandLineProcessorHelpers.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif


namespace {


typedef AbstractLinAlgPack::value_type  Scalar;


} // namespace


int main( int argc, char* argv[] )
{

  using std::endl;
  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::OSTab;
  using Teuchos::Array;
  using Teuchos::includesVerbLevel;
  using Teuchos::ParameterList;
  using Teuchos::CommandLineProcessor;
  typedef Teuchos::ParameterList::PrintOptions PLPO;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Thyra::DefaultMultiVectorProductVectorSpace<Scalar> DMVPVS;
  using MoochoPack::MoochoSolver;
  using MoochoPack::MoochoThyraSolver;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  RCP<Epetra_Comm> epetra_comm;
#ifdef HAVE_MPI
  epetra_comm = rcp( new Epetra_MpiComm(MPI_COMM_WORLD) );
#else
  epetra_comm = rcp( new Epetra_SerialComm );
#endif // HAVE_MPI
  
  bool success = true, result;

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  
  // Rythmos uses long names so enlarge the size of this!
  out->setMaxLenLinePrefix(20);

  try {
  
    // Create the solver objects
    Stratimikos::DefaultLinearSolverBuilder lowsfCreator;
    MoochoThyraSolver moochoThyraSolver;

    //
    // A) Get options from the command line and get the parameter list
    //

    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);

    moochoThyraSolver.setupCLP(&clp);

    std::string paramsFileName = "";
    clp.setOption( "params-file", &paramsFileName,
      "File name for varaious parameter lists in xml format.",
      true );

    bool useBDF = false;
    clp.setOption( "use-BDF", "use-BE", &useBDF,
      "Use BDF or Backward Euler (BE)" );

    double finalTime = 1e-3;
    clp.setOption( "final-time", &finalTime,
      "Final integration time (initial time is 0.0)" );
    
    int numObservationPoints= 10;
    clp.setOption( "num-observations", &numObservationPoints,
      "Number of state observations to take as basis for inversion" );

    double scaleParamsBy = 1e-2;
    clp.setOption( "scale-params-by", &scaleParamsBy,
      "Amount to scale the parameters by to perturb them" );

    Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_DEFAULT;
    setVerbosityLevelOption( "verb-level", &verbLevel,
      "Top-level verbosity level.  By default, this gets deincremented"
      " as you go deeper into numerical objects.",
      &clp );

    bool printValidOptions = false;
    clp.setOption( "print-valid-options", "no-print-valid-options",
      &printValidOptions, "Print the valid options or not" );

    bool testTransientModelEvaluator = true;
    clp.setOption( "test-transient-model", "no-test-transient-model", &testTransientModelEvaluator,
      "Test the Rythmos::ForwardSensitivityIntegratorAsModelEvaluator object or not." );

    double maxSensError = 1e-4;
    clp.setOption( "max-sens-error", &maxSensError,
      "The maximum allowed error in the integrated sensitivity in relation to"
      " the finite-difference sensitivity" );

    bool doInverseProblem = true;
    clp.setOption( "do-inv-prob", "no-do-inv-prob", &doInverseProblem,
      "Run and test the inverse problem or not" );

    double p_inv_err_tol = 1e-8;
    clp.setOption( "max-p-inv-error", &p_inv_err_tol,
      "The maximum allowed error in the inverted for parameters" );
    
    CommandLineProcessor::EParseCommandLineReturn
      parse_return = clp.parse(argc,argv,&std::cerr);

    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL )
      return parse_return;

    moochoThyraSolver.readParameters(&*out);

    RCP<ParameterList>
      paramList = Teuchos::parameterList();
    if (paramsFileName.length())
      updateParametersFromXmlFile( paramsFileName, paramList.ptr() );
    *out << "\nRead in parameters list:\n";
    paramList->print(*out, PLPO().indent(2).showTypes(true).showDoc(true));
    
    //
    // B) Create the Thyra::ModelEvaluator object for the ODE
    //
    
    *out << "\nCreate the EpetraExt::DiagonalTransientModel wrapper object ...\n";

    RCP<EpetraExt::DiagonalTransientModel>
      epetraStateModel = EpetraExt::diagonalTransientModel(
        epetra_comm,
        sublist(paramList,"DiagonalTransientModel",true)
        );

    if (printValidOptions) {
      *out <<"\nepetraStateModel valid options:\n";
      epetraStateModel->getValidParameters()->print(
        *out, PLPO().indent(2).showTypes(true).showDoc(true)
        );
    }

    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
    linearSolverBuilder.setParameterList(sublist(paramList,"Stratimikos",true));
    RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >
      W_factory = createLinearSolveStrategy(linearSolverBuilder);

    RCP<Thyra::ModelEvaluator<Scalar> >
      stateModel = epetraModelEvaluator(epetraStateModel,W_factory);
    
    //
    // C) Create and initialize the state stepper and state integrator objects
    //

    RCP<ParameterList>
      nonlinearSolverPL = sublist(paramList,"TimeStepNonlinearSolver",true);
    //nonlinearSolverPL->get("Default Tol",1e-3*maxStateError); // Set default if not set
    RCP<Rythmos::TimeStepNonlinearSolver<Scalar> >
      nonlinearSolver = rcp(new Rythmos::TimeStepNonlinearSolver<Scalar>());
    nonlinearSolver->setParameterList(nonlinearSolverPL);

    RCP<Rythmos::StepperBase<Scalar> > stateStepper;
    if (useBDF) {
      stateStepper = rcp(
        new Rythmos::ImplicitBDFStepper<Scalar>(
          stateModel, nonlinearSolver
          )
        );
    }
    else {
      stateStepper = rcp(
        new Rythmos::BackwardEulerStepper<Scalar>(
          stateModel, nonlinearSolver
          )
        );
    }
    stateStepper->setParameterList(sublist(paramList,"Rythmos Stepper",true));

    *out <<"\nstateStepper: " << describe(*stateStepper,verbLevel);
    if (printValidOptions) {
      *out <<"\nstateStepper valid options:\n";
      stateStepper->getValidParameters()->print(
        *out, PLPO().indent(2).showTypes(true).showDoc(true)
        );
    }
      
    RCP<Rythmos::IntegratorBase<Scalar> > integrator;
    {
      integrator = Rythmos::controlledDefaultIntegrator<Scalar>(
        Rythmos::simpleIntegrationControlStrategy<Scalar>(
          sublist(paramList,"Rythmos Integration Control",true)
          )
        );
    }
    
    const MEB::InArgs<Scalar>
      state_ic = stateModel->getNominalValues();
    *out << "\nstate_ic: " << describe(state_ic,verbLevel);

    RCP<const Thyra::VectorBase<Scalar> >
      p = state_ic.get_p(0);
    
    stateStepper->setInitialCondition(state_ic);
    integrator->setStepper(stateStepper,finalTime);
      
    if (printValidOptions) {
      *out <<"\nintegrator valid options:\n";
      integrator->getValidParameters()->print(
        *out, PLPO().indent(2).showTypes(true).showDoc(true)
        );
    }

    //
    // D) Generate the matching vectors by running the forward problem through
    //    the state integrator
    //

    Array<Scalar> observationTimes;
    Array<RCP<const Thyra::VectorBase<Scalar> > > stateObservations;
    {
      const double obs_dt = (finalTime - state_ic.get_t())/numObservationPoints;
      *out << "\nCollecting observations at intervals of obs_dt = "
           << obs_dt << " ...\n";
      OSTab tab(out);
      int obs_idx;
      double curr_t;
      for (
        obs_idx = 0, curr_t = state_ic.get_t() + obs_dt; 
        obs_idx < numObservationPoints;
        ++obs_idx, curr_t += obs_dt
        )
      {
        // Make sure we don't step over final time
        curr_t = std::min(curr_t,finalTime);

        *out << "\nCollecting state observation at curr_t = " << curr_t << " ...\n";
        observationTimes.push_back(curr_t);
        stateObservations.push_back(get_fwd_x(*integrator,curr_t)->clone_v());
        if (includesVerbLevel(verbLevel,Teuchos::VERB_EXTREME)) {
          *out << "\nx("<<curr_t<<") = "
               << Teuchos::describe(*stateObservations.back(),Teuchos::VERB_EXTREME);
        }
      }
    }
    
    //
    // F) Create the forward sensitivity stepper and reduced sensitivity model
    //    evaluator
    //
    
    //
    // F.1) Create the forward sensitivity stepper
    //
    
    RCP<Rythmos::ForwardSensitivityStepper<Scalar> >
      stateAndSensStepper = Rythmos::forwardSensitivityStepper<Scalar>(
        stateModel, 0, stateModel->getNominalValues(),
        stateStepper, nonlinearSolver
        );
    // The above call will result in stateStepper and nonlinearSolver being
    // cloned.  This helps to ensure consistency between the state and
    // sensitivity computations!

    const RCP<const DMVPVS>
      s_bar_space = rcp_dynamic_cast<const DMVPVS>(
        stateAndSensStepper->getFwdSensModel()->get_x_space()
        );

    //
    // F.2) Set the initial condition for the state and forward sensitivities
    //

    RCP<Thyra::VectorBase<Scalar> > s_bar_init
      = createMember(stateAndSensStepper->getFwdSensModel()->get_x_space());
    assign( s_bar_init.ptr(), 0.0 );
    RCP<Thyra::VectorBase<Scalar> > s_bar_dot_init
      = createMember(stateAndSensStepper->getFwdSensModel()->get_x_space());
    assign( s_bar_dot_init.ptr(), 0.0 );
    // Above, I believe that these are the correct initial conditions for
    // s_bar and s_bar_dot given how the EpetraExt::DiagonalTransientModel
    // is currently implemented!

    RCP<const Rythmos::StateAndForwardSensitivityModelEvaluator<Scalar> >
      stateAndSensModel = stateAndSensStepper->getStateAndFwdSensModel();

    MEB::InArgs<Scalar>
      state_and_sens_ic = stateAndSensStepper->getModel()->createInArgs();

    // Copy time, parameters etc.
    state_and_sens_ic.setArgs(state_ic);
    // Set initial condition for x_bar = [ x; s_bar ]
    state_and_sens_ic.set_x(
      stateAndSensModel->create_x_bar_vec(state_ic.get_x(),s_bar_init)
      );
    // Set initial condition for x_bar_dot = [ x_dot; s_bar_dot ]
    state_and_sens_ic.set_x_dot(
      stateAndSensModel->create_x_bar_vec(state_ic.get_x_dot(),s_bar_dot_init)
      );
    
    *out << "\nstate_and_sens_ic:\n" << describe(state_and_sens_ic,verbLevel);
    
    stateAndSensStepper->setInitialCondition(state_and_sens_ic);

    //
    // F.3) Setup the state-matching response functions
    //

    Array<RCP<const Thyra::ModelEvaluator<Scalar> > > matchingFuncs;
    Array<Thyra::ModelEvaluatorBase::InArgs<Scalar> > matchingFuncBasePoints;

    for (int obs_idx = 0; obs_idx < numObservationPoints; ++obs_idx ) {
      RCP<Thyra::DefaultInverseModelEvaluator<Scalar> >
        matchingFunc = Thyra::defaultInverseModelEvaluator(stateModel);
      matchingFunc->setParameterList(
        sublist(paramList,"DefaultInverseModelEvalautor",true));
      matchingFunc->set_observationTarget(stateObservations[obs_idx]);
      matchingFuncs.push_back(matchingFunc);
      matchingFuncBasePoints.push_back(state_ic);
    }

    //
    // F.4) Setup the final integrator that will evaluate the reduced response
    //      function and it's derivatives
    
    namespace FSIAMET = Rythmos::ForwardSensitivityIntegratorAsModelEvaluatorTypes;
    const RCP<Rythmos::IntegratorBase<Scalar> > nullIntegrator;
    RCP<Rythmos::ForwardSensitivityIntegratorAsModelEvaluator<Scalar> >
      sensIntegatorAsModelEvaluator =
      Rythmos::forwardSensitivityIntegratorAsModelEvaluator<Scalar>(
        stateStepper, integrator, stateAndSensStepper, nullIntegrator,
        state_and_sens_ic, observationTimes, matchingFuncs,
        matchingFuncBasePoints, FSIAMET::RESPONSE_TYPE_SUM
        );
    sensIntegatorAsModelEvaluator->setParameterList(
      sublist(paramList,"ForwardSensitivityIntegratorAsModelEvaluator",true) );
    
    //
    // G) Test the ForwardSensitivityIntegratorAsModelEvaluator object
    //

    if (testTransientModelEvaluator) {

      //
      // G.1) Test that the base evaluation computes a zero response function.
      //

      // Here, we check that if we evaluation p -> g(p) for the base p, then we
      // should get zero since the state should match the read observations
      // exactly.  We also check that the computed state is exactly the same
      // also.  We can cheat and get the final state from stateStepper which is
      // used internally.
      //

      *out << "\nTest that sensIntegatorAsModelEvaluator computes the right base point ... \n";
      {
      
        RCP<Thyra::VectorBase<Scalar> >
          d_hat = createMember(sensIntegatorAsModelEvaluator->get_g_space(0));

        eval_g(
          *sensIntegatorAsModelEvaluator,
          0, *state_ic.get_p(0),
          0, d_hat.ptr()
          );

        *out
          << "\nd_hat(p) evaluated using sensIntegatorAsModelEvaluator:\n"
          << *d_hat;

        // Test that d_hat == 0.0 to all digits since 

        const Scalar d_hat_nrm = norm(*d_hat);

        const bool d_hat_norm_is_zero = (d_hat_nrm == 0.0);

        *out << "\nCheck: ||d_hat|| = " << d_hat_nrm << " == 0 : "
             << Thyra::passfail(d_hat_norm_is_zero) << endl;

        if (!d_hat_norm_is_zero)
          success = false;

        // Test that:
        //
        //  get_x(stateStepper,finalTime) == stateObservations[numObservationPoints-1]
        //
        // to the given precision.

        *out
          << "\nChecking that get_x(stateStepper,finalTime) == stateObservations[numObservationPoints-1] ... \n";

        RCP<const Thyra::VectorBase<Scalar> >
          x_final = get_x(*stateStepper,observationTimes[numObservationPoints-1]);
      
        *out << endl;
        result = Thyra::testRelNormDiffErr(
          "get_x(stateStepper,finalTime)", *x_final,
          "stateObservations[numObservationPoints-1]", *stateObservations[numObservationPoints-1],
          "errorTol", 0.0, "warningTol", 1e-8,
          &OSTab(out).o(), verbLevel
          );
        if (!result) success = false;

      }

      *out << "\nBase parameters p = " << *p;

      RCP<Thyra::VectorBase<Scalar> >
        p_perturb = createMember(p->space());

      V_StV( p_perturb.ptr(), scaleParamsBy, *p );

      *out << "\nPerturbed parameters p_perturb = " << *p_perturb;

      //
      // G.2) Test the sensitivity D(d_hat)/d(p) by finite differences
      //

      //
      *out << "\nCompute analytic D(d_hat)/d(p) at p_perturb ...\n";
      //

      MEB::DerivativeMultiVector<Scalar> D_d_hat_D_p_exact;
      {

        D_d_hat_D_p_exact = Thyra::create_DgDp_mv(
          *sensIntegatorAsModelEvaluator, 0, 0, MEB::DERIV_TRANS_MV_BY_ROW );

        MEB::InArgs<Scalar>
          inArgs = sensIntegatorAsModelEvaluator->createInArgs();
        inArgs.set_p(0,p_perturb);

        MEB::OutArgs<Scalar>
          outArgs = sensIntegatorAsModelEvaluator->createOutArgs();
        outArgs.set_DgDp(0,0,D_d_hat_D_p_exact);

        sensIntegatorAsModelEvaluator->evalModel(inArgs,outArgs);
      
        *out
          << "\nAnalytic D(d_hat)/d(p)^T = "
          << describe(*D_d_hat_D_p_exact.getMultiVector(),verbLevel);

      }
 
      //
      *out << "\nCompute finite-diff D(d_hat)/d(p) at p_perturb ...\n";
      //

      MEB::DerivativeMultiVector<Scalar> D_d_hat_D_p_fd;
      {

        RCP<Thyra::DirectionalFiniteDiffCalculator<Scalar> > fdCalc =
          Thyra::directionalFiniteDiffCalculator<Scalar>(
            sublist(paramList,"DirectionalFiniteDiffCalculator",true)
            );

        RCP<Thyra::DefaultFiniteDifferenceModelEvaluator<Scalar> > fdModelEval =
          Thyra::defaultFiniteDifferenceModelEvaluator<Scalar>(
            sensIntegatorAsModelEvaluator,
            fdCalc
            );

        MEB::InArgs<Scalar>
          inArgs = fdModelEval->createInArgs();
        inArgs.set_p(0,p_perturb);
      
        D_d_hat_D_p_fd = Thyra::create_DgDp_mv(
          *fdModelEval, 0, 0, MEB::DERIV_TRANS_MV_BY_ROW
          );
      
        MEB::OutArgs<Scalar>
          outArgs = fdModelEval->createOutArgs();
        outArgs.set_DgDp( 0, 0, D_d_hat_D_p_fd );
      
        // Silence the model evaluators that are called.  The fdCal object
        // will show all of the inputs and outputs for each call.
        stateStepper->setVerbLevel(Teuchos::VERB_NONE);
        sensIntegatorAsModelEvaluator->setVerbLevel(Teuchos::VERB_NONE);

        fdModelEval->evalModel(inArgs,outArgs);
      
        *out
          << "\nFinite difference D(d_hat)/d(p)^T = "
          << describe(*D_d_hat_D_p_fd.getMultiVector(),verbLevel);
      
      }
    
      //
      *out << "\nCompare analytic and finite-diff D(d_hat)/d(p)^T ...\n";
      //

      RCP<const Thyra::VectorBase<Scalar> >
        D_d_hat_D_p_exact_vec = D_d_hat_D_p_exact.getMultiVector()->col(0),
        D_d_hat_D_p_fd_vec = D_d_hat_D_p_fd.getMultiVector()->col(0);
    
      result = Thyra::testRelNormDiffErr(
        "D_d_hat_D_p_exact_vec", *D_d_hat_D_p_exact_vec,
        "D_d_hat_D_p_fd_vec", *D_d_hat_D_p_fd_vec,
        "maxSensError", maxSensError,
        "warningTol", 1.0, // Don't warn
        &*out, verbLevel
        );
      if (!result) success = false;

    }

    if (doInverseProblem) {
    
      //
      // H) Setup the MoochoThyraSolver object
      //
      
      *out << "\nSetup the MoochoThyraSolver object ...\n";
      
      moochoThyraSolver.setParameterList(
        sublist(paramList,"MoochoThyraSolver",true) );
      
      moochoThyraSolver.setModel(sensIntegatorAsModelEvaluator);
      
      //
      // I) Solve the transient inverse problem using MOOCHO.
      //
      
      *out << "\nSetup the transient inverse problem ...\n";
      
      moochoThyraSolver.readInitialGuess(&*out);
      
      const MoochoSolver::ESolutionStatus
        solution_status = moochoThyraSolver.solve();
      
      if ( solution_status != MoochoSolver::SOLVE_RETURN_SOLVED )
        success = false;
      
      moochoThyraSolver.writeFinalSolution(&*out);
      
      //
      // J) Check the inverse solution against the true solution
      //
      
      RCP<const Thyra::VectorBase<Scalar> >
        p_inv = moochoThyraSolver.getFinalPoint().get_p(0);
      
      result = Thyra::testRelNormDiffErr(
        "p_inv", *p_inv,
        "p", *p,
        "errorTol", p_inv_err_tol,
        "warningTol", p_inv_err_tol * 1e+2,
        &OSTab(out).o(), Teuchos::VERB_EXTREME
        );
      if (!result) success = false;

    }
      
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,success);

  if (success)
    *out << "\nEnd Result: TEST PASSED\n";
  else
    *out << "\nEnd Result: TEST FAILED\n";
  
  return ( success ? 0 : 1 );

}
