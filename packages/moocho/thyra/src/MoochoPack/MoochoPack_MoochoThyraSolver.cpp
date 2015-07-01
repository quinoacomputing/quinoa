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

#include "MoochoPack_MoochoThyraSolver.hpp"
#include "NLPInterfacePack_NLPDirectThyraModelEvaluator.hpp"
#include "NLPInterfacePack_NLPFirstOrderThyraModelEvaluator.hpp"
#include "Thyra_DefaultFiniteDifferenceModelEvaluator.hpp"
#include "Thyra_DefaultStateEliminationModelEvaluator.hpp"
#include "Thyra_DefaultEvaluationLoggerModelEvaluator.hpp"
#include "Thyra_DefaultInverseModelEvaluator.hpp"
#include "Thyra_DefaultLumpedParameterModelEvaluator.hpp"
#include "Thyra_DefaultSpmdMultiVectorFileIO.hpp"
#include "Thyra_DampenedNewtonNonlinearSolver.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

namespace {

//
// ParameterList parameters and sublists
//

const std::string SolveMode_name = "Solve Mode";
const Teuchos::RCP<
  Teuchos::StringToIntegralParameterEntryValidator<
    MoochoPack::MoochoThyraSolver::ESolveMode
  >
>
solveModeValidator = Teuchos::rcp(
  new Teuchos::StringToIntegralParameterEntryValidator<MoochoPack::MoochoThyraSolver::ESolveMode>(
    Teuchos::tuple<std::string>(
      "Forward Solve"
      ,"Optimize"
      )
    ,Teuchos::tuple<std::string>(
      "Only solve state equaitons f(x,p)=0 for states x\n"
      "given fixed parameters values p."
      ,"Solve the simulation constrained optimization problem\n"
      "  min  g(x,p)\n"
      "  s.t. f(x,p)=0\n"
      "for the state varaibles x and parameters p."
      )
    ,Teuchos::tuple<MoochoPack::MoochoThyraSolver::ESolveMode>(
      MoochoPack::MoochoThyraSolver::SOLVE_MODE_FORWARD
      ,MoochoPack::MoochoThyraSolver::SOLVE_MODE_OPTIMIZE
      )
    ,""
    )
  );
const std::string SolveMode_default = "Optimize";

const std::string NLPType_name = "NLP Type";
const Teuchos::RCP<
  Teuchos::StringToIntegralParameterEntryValidator<
    MoochoPack::MoochoThyraSolver::ENLPType
    >
  >
nlpTypeValidator = Teuchos::rcp(
  new Teuchos::StringToIntegralParameterEntryValidator<MoochoPack::MoochoThyraSolver::ENLPType>(
    Teuchos::tuple<std::string>(
      "First Order"
      ,"Direct"
      )
    ,Teuchos::tuple<std::string>(
      "Support the NLPInterfacePack::NLPFirstOrder interface which assumes\n"
      "that full adjoints for the objective and constraint derivatives are\n"
      "available."
      ,"Support the NLPInterfacePack::NLPDirect interface which only assumes\n"
      "that forward or direct sensitivities and state solves are supported."
      )
    ,Teuchos::tuple<MoochoPack::MoochoThyraSolver::ENLPType>(
      MoochoPack::MoochoThyraSolver::NLP_TYPE_FIRST_ORDER
      ,MoochoPack::MoochoThyraSolver::NLP_TYPE_DIRECT
      )
    ,""
    )
  );
const std::string NLPType_default = "First Order";

const std::string NonlinearlyEliminateStates_name = "Nonlinearly Eliminate States";
const bool NonlinearlyEliminateStates_default = false;

const std::string UseFiniteDifferencesForObjective_name = "Use Finite Differences For Objective";
const bool UseFiniteDifferencesForObjective_default = false;

const std::string ObjectiveFiniteDifferenceSettings_name = "Objective Finite Difference Settings";

const std::string UseFiniteDifferencesForConstraints_name = "Use Finite Differences For Constraints";
const bool UseFiniteDifferencesForConstraints_default = false;

const std::string ConstraintsFiniteDifferenceSettings_name = "Constraints Finite Difference Settings";

const std::string FwdNewtonTol_name = "Forward Newton Tolerance";
const double FwdNewtonTol_default = -1.0;

const std::string FwdNewtonMaxIters_name = "Forward Newton Max Iters";
const int FwdNewtonMaxIters_default = 20;

const std::string ForwardNewtonDampening_name = "Forward Newton Dampening";
const bool ForwardNewtonDampening_default = true;

const std::string FwdNewtonMaxLineSearchIters_name = "Forward Newton Max Line Search Iters";
const int FwdNewtonMaxLineSearchIters_default = 20;

const std::string UseBuiltInInverseObjectiveFunction_name = "Use Built-in Inverse Objective Function";
const bool UseBuiltInInverseObjectiveFunction_default = false;

const std::string InverseObjectiveFunctionSettings_name = "Inverse Objective Function Settings";

const std::string UseParameterLumping_name = "Use Parameter Lumping";
const bool UseParameterLumping_default = false;

const std::string LumpedParametersSettings_name = "Lumped Parameters Settings";

const std::string OutputFileTag_name = "Output File Tag";
const std::string OutputFileTag_default = "";

const std::string OutputOnAllProcesses_name = "Output On All Processes";
const bool OutputOnAllProcesses_default = false;

const std::string ShowModelEvaluatorTrace_name = "Show Model Evaluator Trace";
const bool ShowModelEvaluatorTrace_default = "false";

const std::string StateGuess_name = "State Guess";

const std::string ParamGuess_name = "Parameter Guess";

const std::string ParamLowerBounds_name = "Parameter Lower Bounds";

const std::string ParamUpperBounds_name = "Parameter Upper Bounds";

const std::string StateSoluFileBaseName_name = "State Solution File Base Name";
const std::string StateSoluFileBaseName_default = "";

const std::string ParamSoluFileBaseName_name = "Parameters Solution File Base Name";
const std::string ParamSoluFileBaseName_default = "";

} // namespace

namespace MoochoPack {

// Constructors/initialization

MoochoThyraSolver::MoochoThyraSolver(
  const std::string    &paramsXmlFileName
  ,const std::string   &extraParamsXmlString
  ,const std::string   &paramsUsedXmlOutFileName
  ,const std::string   &paramsXmlFileNameOption
  ,const std::string   &extraParamsXmlStringOption
  ,const std::string   &paramsUsedXmlOutFileNameOption
  )
  :paramsXmlFileName_(paramsXmlFileName)
  ,extraParamsXmlString_(extraParamsXmlString)
  ,paramsUsedXmlOutFileName_(paramsUsedXmlOutFileName)
  ,paramsXmlFileNameOption_(paramsXmlFileNameOption)
  ,extraParamsXmlStringOption_(extraParamsXmlStringOption)
  ,paramsUsedXmlOutFileNameOption_(paramsUsedXmlOutFileNameOption)
  ,stateVectorIO_(Teuchos::rcp(new Thyra::DefaultSpmdMultiVectorFileIO<value_type>))
  ,parameterVectorIO_(Teuchos::rcp(new Thyra::DefaultSpmdMultiVectorFileIO<value_type>))
  ,solveMode_(SOLVE_MODE_OPTIMIZE)
  ,nlpType_(NLP_TYPE_FIRST_ORDER)
  ,nonlinearlyElimiateStates_(false)
  ,use_finite_diff_for_obj_(false)
  ,use_finite_diff_for_con_(false)
  ,fwd_newton_tol_(-1.0)
  ,fwd_newton_max_iters_(20)
  ,fwd_newton_dampening_(false)
  ,fwd_newton_max_ls_iters_(20)
  ,useInvObjFunc_(false)
  ,useParameterLumping_(false)
  ,outputFileTag_("")
  ,showModelEvaluatorTrace_(false)
  ,stateSoluFileBase_("")
  ,paramSoluFileBase_("")
{}

MoochoThyraSolver::~MoochoThyraSolver()
{}

void MoochoThyraSolver::setupCLP(
  Teuchos::CommandLineProcessor *clp
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(0==clp);
  solver_.setup_commandline_processor(clp);
  clp->setOption(
    paramsXmlFileNameOption().c_str(),&paramsXmlFileName_
    ,"Name of an XML file containing parameters for linear solver options to be appended first."
    );
  clp->setOption(
    extraParamsXmlStringOption().c_str(),&extraParamsXmlString_
    ,"An XML string containing linear solver parameters to be appended second."
    );
  clp->setOption(
    paramsUsedXmlOutFileNameOption().c_str(),&paramsUsedXmlOutFileName_
    ,"Name of an XML file that can be written with the parameter list after it has been used on completion of this program."
    );
}

void MoochoThyraSolver::readParameters( std::ostream *out_arg )
{
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::getFancyOStream(Teuchos::rcp(out_arg,false));
  Teuchos::OSTab tab(out);
  if(out.get()) *out << "\nMoochoThyraSolver::readParameters(...):\n";
  Teuchos::OSTab tab2(out);
  Teuchos::RCP<Teuchos::ParameterList>
    paramList = this->getNonconstParameterList();
  if(!paramList.get()) {
    if(out.get()) *out << "\nCreating a new Teuchos::ParameterList ...\n";
    paramList = Teuchos::rcp(new Teuchos::ParameterList("MoochoThyraSolver"));
  }
  if(paramsXmlFileName().length()) {
    if(out.get()) *out << "\nReading parameters from XML file \""<<paramsXmlFileName()<<"\" ...\n";
    Teuchos::updateParametersFromXmlFile(paramsXmlFileName(), paramList.ptr());
  }
  if(extraParamsXmlString().length()) {
    if(out.get())
      *out << "\nAppending extra parameters from the XML string \""<<extraParamsXmlString()<<"\" ...\n";
    Teuchos::updateParametersFromXmlString(extraParamsXmlString(), paramList.ptr());
  }
  if( paramsXmlFileName().length() || extraParamsXmlString().length() ) {
    typedef Teuchos::ParameterList::PrintOptions PLPrintOptions;
    if(out.get()) {
      *out  << "\nUpdated parameter list:\n";
      paramList->print(
        *out,PLPrintOptions().indent(2).showTypes(true)
        );
    }
    this->setParameterList(paramList);
  }
}

// Overridden from ParameterListAcceptor

void MoochoThyraSolver::setParameterList(
  Teuchos::RCP<Teuchos::ParameterList> const& paramList
  )
{
  typedef MoochoPack::MoochoSolver MS;
  TEUCHOS_TEST_FOR_EXCEPT(!paramList.get());
  paramList->validateParameters(*getValidParameters(),0); // Just validate my level!
  paramList_ = paramList;
  solveMode_ = solveModeValidator->getIntegralValue(
    *paramList_,SolveMode_name,SolveMode_default);
  nlpType_ = nlpTypeValidator->getIntegralValue(
    *paramList_,NLPType_name,NLPType_default);
  nonlinearlyElimiateStates_ = paramList_->get(
    NonlinearlyEliminateStates_name,NonlinearlyEliminateStates_default);
  use_finite_diff_for_obj_ = paramList_->get(
    UseFiniteDifferencesForObjective_name,UseFiniteDifferencesForObjective_default);
  use_finite_diff_for_con_ = paramList_->get(
    UseFiniteDifferencesForConstraints_name,UseFiniteDifferencesForConstraints_default);
  fwd_newton_tol_ = paramList_->get(
    FwdNewtonTol_name,FwdNewtonTol_default);
  fwd_newton_max_iters_ = paramList_->get(
    FwdNewtonMaxIters_name,FwdNewtonMaxIters_default);
  fwd_newton_dampening_ = paramList_->get(
    ForwardNewtonDampening_name,ForwardNewtonDampening_default);
  fwd_newton_max_ls_iters_ = paramList_->get(
    FwdNewtonMaxLineSearchIters_name,FwdNewtonMaxLineSearchIters_default);
  useInvObjFunc_ = paramList_->get(
    UseBuiltInInverseObjectiveFunction_name,UseBuiltInInverseObjectiveFunction_default);
  useParameterLumping_ = paramList_->get(
    UseParameterLumping_name, UseParameterLumping_default);
  outputFileTag_ = paramList->get(
    OutputFileTag_name,OutputFileTag_default);
  solver_.set_output_file_tag(outputFileTag_);
  solver_.set_output_context("",
    paramList_->get(OutputOnAllProcesses_name, OutputOnAllProcesses_default)
    ? MS::OUTPUT_TO_BLACK_HOLE_FALSE
    : MS::OUTPUT_TO_BLACK_HOLE_DEFAULT
    );
  showModelEvaluatorTrace_ = paramList->get(
    ShowModelEvaluatorTrace_name,ShowModelEvaluatorTrace_default);
  x_reader_.setParameterList(
    sublist(paramList_,StateGuess_name)
    );
  p_reader_.setParameterList(
    sublist(paramList_,ParamGuess_name)
    );
  p_l_reader_.setParameterList(
    sublist(paramList_,ParamLowerBounds_name)
    );
  p_u_reader_.setParameterList(
    sublist(paramList_,ParamUpperBounds_name)
    );
  stateSoluFileBase_ = paramList_->get(
    StateSoluFileBaseName_name,StateSoluFileBaseName_default);
  paramSoluFileBase_ = paramList_->get(
    ParamSoluFileBaseName_name,ParamSoluFileBaseName_default);
#ifdef TEUCHOS_DEBUG
  paramList->validateParameters(*getValidParameters(),0); // Just validate my level!
#endif
}

Teuchos::RCP<Teuchos::ParameterList>
MoochoThyraSolver::getNonconstParameterList()
{
  return paramList_;
}

Teuchos::RCP<Teuchos::ParameterList>
MoochoThyraSolver::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}

Teuchos::RCP<const Teuchos::ParameterList>
MoochoThyraSolver::getParameterList() const
{
  return paramList_;
}

Teuchos::RCP<const Teuchos::ParameterList>
MoochoThyraSolver::getValidParameters() const
{
  static Teuchos::RCP<Teuchos::ParameterList> pl;
  if(pl.get()==NULL) {
    pl = Teuchos::rcp(new Teuchos::ParameterList());
    pl->set(
      SolveMode_name,SolveMode_default
      ,"The type of solve to perform."
      ,solveModeValidator
      );
    pl->set(
      NLPType_name,NLPType_default
      ,"The type of MOOCHO NLP subclass to use."
      ,nlpTypeValidator
      );
    pl->set(
      NonlinearlyEliminateStates_name,NonlinearlyEliminateStates_default
      ,"If true, then the model's state equations and state variables\n"
      "are nonlinearlly eliminated using a forward solver."
      );
    pl->set(
      UseFiniteDifferencesForObjective_name,UseFiniteDifferencesForObjective_default
      ,"Use finite differences for missing objective function derivatives (Direct NLP only).\n"
      "See the options in the sublist \"" + ObjectiveFiniteDifferenceSettings_name + "\"."
      );
    {
      Thyra::DirectionalFiniteDiffCalculator<Scalar> dfdcalc;
      {
        Teuchos::ParameterList
          &fdSublist = pl->sublist(ObjectiveFiniteDifferenceSettings_name);
        fdSublist.setParameters(*dfdcalc.getValidParameters());
      }
      pl->set(
        UseFiniteDifferencesForConstraints_name,UseFiniteDifferencesForConstraints_default
        ,"Use  finite differences for missing constraint derivatives (Direct NLP only).\n"
        "See the   options in the sublist \"" + ConstraintsFiniteDifferenceSettings_name + "\"."
        );
      {
        Teuchos::ParameterList
          &fdSublist = pl->sublist(ConstraintsFiniteDifferenceSettings_name);
        fdSublist.setParameters(*dfdcalc.getValidParameters());
      }
    }
    pl->set(
      FwdNewtonTol_name,FwdNewtonTol_default
      ,"Tolarance used for the forward state solver in eliminating\n"
      "the state equations/variables."
      );
    pl->set(
      FwdNewtonMaxIters_name,FwdNewtonMaxIters_default
      ,"Maximum number of iterations allowed for the forward state\n"
      "solver in eliminating the state equations/variables."
      );
    pl->set(
      ForwardNewtonDampening_name, ForwardNewtonDampening_default,
      "If true, then the state elimination nonlinear solver will\n"
      "use a dampened line search.  Otherwise, it will just take fulls steps."
      );
    pl->set(
      FwdNewtonMaxLineSearchIters_name, FwdNewtonMaxLineSearchIters_default,
      "Maximum number of linea search iterations per newton iteration\n"
      "allowed for the forward state solver in eliminating the state equations/variables."
      );
    pl->set(
      UseBuiltInInverseObjectiveFunction_name,UseBuiltInInverseObjectiveFunction_default
      ,"Use a built-in form of a simple inverse objection function instead\n"
      "of a a response function contained in the underlying model evaluator\n"
      "object itself.  The settings are contained in the sublist\n"
      "\""+InverseObjectiveFunctionSettings_name+"\".\n"
      "Note that this feature allows the client to form a useful type\n"
      "of optimization problem just with a model that supports only the\n"
      "parameterized state function f(x,p)=0."
      );
    {
      Teuchos::RCP<Thyra::DefaultInverseModelEvaluator<Scalar> >
        inverseModel = rcp(new Thyra::DefaultInverseModelEvaluator<Scalar>());
      pl->sublist(
        InverseObjectiveFunctionSettings_name,false
        ,"Settings for the built-in inverse objective function.\n"
        "See the outer parameter \""+UseBuiltInInverseObjectiveFunction_name+"\"."
        ).setParameters(*inverseModel->getValidParameters());
    }
    pl->set(
      UseParameterLumping_name,UseParameterLumping_default,
      "Use a reduced basis to lump optimization parameters as\n"
      "p_orig = P_basis * p.  If set to true, then the settings\n"
      "in \""+LumpedParametersSettings_name+"\" determine how the\n"
      "parameter basis is set.  This feature can be used to safely\n"
      "regularize a problem if there are linearly dependent parameters\n"
      "and will generally speed up the optimiztation algorithms."
      );
    {
      Teuchos::RCP<Thyra::DefaultLumpedParameterModelEvaluator<Scalar> >
        lumpedParamModel = rcp(new Thyra::DefaultLumpedParameterModelEvaluator<Scalar>);
      pl->sublist(
        LumpedParametersSettings_name,false
        ,"Settings for parameter lumping.\n"
        "See the outer parameter \""+UseParameterLumping_name+"\"."
        ).setParameters(*lumpedParamModel->getValidParameters());
    }
    pl->set(OutputFileTag_name,OutputFileTag_default,
      "A tag that is attached to every output file that is created by the\n"
      "solver.  If empty \"\", then no tag is used.  This option simply is\n"
      "passed into the set_output_context(...) function on the underlying\n"
      "MoochoPack::MoochoSolver object.  Therefore, this same parameter\n"
      "can be set in code as well without going through the parameter list.");
    pl->set(OutputOnAllProcesses_name, OutputOnAllProcesses_default,
      "If set to true, then the console output and all MOOCHO files will be\n"
      "written for every process.  This helps in debugging very hard\n"
      "parallel problems.");
    pl->set(ShowModelEvaluatorTrace_name,ShowModelEvaluatorTrace_default
      ,"Determine if a trace of the objective function will be shown or not\n"
      "when the NLP is evaluated."
      );
    if(this->get_stateVectorIO().get())
      x_reader_.set_fileIO(this->get_stateVectorIO());
    pl->sublist(StateGuess_name).setParameters(*x_reader_.getValidParameters());
    if(this->get_parameterVectorIO().get()) {
      p_reader_.set_fileIO(this->get_parameterVectorIO());
      p_l_reader_.set_fileIO(this->get_parameterVectorIO());
      p_u_reader_.set_fileIO(this->get_parameterVectorIO());
      pl->sublist(ParamGuess_name).setParameters(*p_reader_.getValidParameters());
      pl->sublist(ParamLowerBounds_name).setParameters(*p_l_reader_.getValidParameters());
      pl->sublist(ParamUpperBounds_name).setParameters(*p_u_reader_.getValidParameters());
    }
    pl->set(
      StateSoluFileBaseName_name,StateSoluFileBaseName_default
      ,"If specified, a file with this basename will be written to with\n"
      "the final value of the state variables.  A different file for each\n"
      "process will be created.  Note that these files can be used for the\n"
      "initial guess for the state variables."
      );
    pl->set(
      ParamSoluFileBaseName_name,ParamSoluFileBaseName_default
      ,"If specified, a file with this basename will be written to with\n"
      "the final value of the parameters.  A different file for each\n"
      "process will be created.  Note that these files can be used for the\n"
      "initial guess for the parameters."
      );
  }
  return pl;
}

// Misc Access/Setup

void MoochoThyraSolver::setSolveMode( const ESolveMode solveMode )
{
  solveMode_ = solveMode;
}

MoochoThyraSolver::ESolveMode
MoochoThyraSolver::getSolveMode() const
{
  return solveMode_;
}

MoochoSolver& MoochoThyraSolver::getSolver()
{
  return solver_;
}

const MoochoSolver& MoochoThyraSolver::getSolver() const
{
  return solver_;
}

// Model specification, setup, solve, and solution extraction.

void MoochoThyraSolver::setModel(
  const Teuchos::RCP<Thyra::ModelEvaluator<value_type> > &origModel,
  const int p_idx,
  const int g_idx
  )
{

  using Teuchos::rcp;
  using Teuchos::sublist;
  using NLPInterfacePack::NLP;
  using NLPInterfacePack::NLPDirectThyraModelEvaluator;
  using NLPInterfacePack::NLPFirstOrderThyraModelEvaluator;

  origModel_ = origModel;
  p_idx_ = p_idx;
  g_idx_ = g_idx;

  //const int procRank = Teuchos::GlobalMPISession::getRank();
  //const int numProcs = Teuchos::GlobalMPISession::getNProc();

  //
  // Wrap the orginal model in different decorators
  //

  outerModel_ = origModel_;

  // Inverse objective (and parameter regularization)?
  if (useInvObjFunc_) {
    Teuchos::RCP<Thyra::DefaultInverseModelEvaluator<Scalar> >
      inverseModel = Thyra::defaultInverseModelEvaluator<Scalar>(outerModel_);
    inverseModel->setVerbLevel(Teuchos::VERB_LOW);
    inverseModel->set_observationTargetIO(get_stateVectorIO());
    inverseModel->set_parameterBaseIO(get_parameterVectorIO());
    inverseModel->setParameterList(
      Teuchos::sublist(paramList_,InverseObjectiveFunctionSettings_name) );
    outerModel_ = inverseModel; 
    g_idx_ = inverseModel->Ng()-1;
  }
  
  // Evaluation loging
  Teuchos::RCP<std::ostream>
    modelEvalLogOut = Teuchos::fancyOStream(
      solver_.generate_output_file("ModelEvaluationLog")
      );
  Teuchos::RCP<Thyra::DefaultEvaluationLoggerModelEvaluator<Scalar> >
    loggerThyraModel
    = rcp(
      new Thyra::DefaultEvaluationLoggerModelEvaluator<Scalar>(
        outerModel_,modelEvalLogOut
        )
      );
  outerModel_ = loggerThyraModel; 
  
  // Manipulating the nominal values
  nominalModel_
    = rcp(
      new Thyra::DefaultNominalBoundsOverrideModelEvaluator<Scalar>(outerModel_,Teuchos::null)
      );
  outerModel_ = nominalModel_; 

  // Capturing the final point
  finalPointModel_
    = rcp(
      new Thyra::DefaultFinalPointCaptureModelEvaluator<value_type>(outerModel_)
      );
  outerModel_ = finalPointModel_;

  // Parameter lumping?
  if (useParameterLumping_) {
    Teuchos::RCP<Thyra::DefaultLumpedParameterModelEvaluator<Scalar> >
      lumpedParamModel = Thyra::defaultLumpedParameterModelEvaluator<Scalar>(outerModel_);
    lumpedParamModel->setVerbLevel(Teuchos::VERB_LOW);
    lumpedParamModel->setParameterList(
      sublist(paramList_,LumpedParametersSettings_name));
    outerModel_ = lumpedParamModel; 
  }
  // Note, above we put parameter lumping on the very top so that everything
  // like the final point capture and the nominal bounds overrider deal with
  // the original parameters, not the lumped parameters.

  //
  // Create the NLP
  //
    
  Teuchos::RCP<NLP> nlp;

  switch(solveMode_) {
    case SOLVE_MODE_FORWARD: {
      RCP<NLPFirstOrderThyraModelEvaluator>
        nlpFirstOrder = rcp(
          new NLPFirstOrderThyraModelEvaluator(outerModel_,-1,-1)
          );
      nlpFirstOrder->showModelEvaluatorTrace(showModelEvaluatorTrace_);
      nlp = nlpFirstOrder;
      break;
    }
    case SOLVE_MODE_OPTIMIZE: {
      // Setup finite difference object
      RCP<Thyra::DirectionalFiniteDiffCalculator<Scalar> > objDirecFiniteDiffCalculator;
      if(use_finite_diff_for_obj_) {
        objDirecFiniteDiffCalculator = rcp(new Thyra::DirectionalFiniteDiffCalculator<Scalar>());
        if(paramList_.get())
          objDirecFiniteDiffCalculator->setParameterList(
            Teuchos::sublist(paramList_,ObjectiveFiniteDifferenceSettings_name)
            );
      }
      RCP<Thyra::DirectionalFiniteDiffCalculator<Scalar> > conDirecFiniteDiffCalculator;
      if(use_finite_diff_for_con_) {
        conDirecFiniteDiffCalculator = rcp(new Thyra::DirectionalFiniteDiffCalculator<Scalar>());
        if(paramList_.get())
          conDirecFiniteDiffCalculator->setParameterList(
            Teuchos::sublist(paramList_,ConstraintsFiniteDifferenceSettings_name)
            );
      }
      if( nonlinearlyElimiateStates_ ) {
        // Create a Thyra::NonlinearSolverBase object to solve and eliminate the
        // state variables and the state equations
        Teuchos::RCP<Thyra::DampenedNewtonNonlinearSolver<Scalar> >
          stateSolver = rcp(new Thyra::DampenedNewtonNonlinearSolver<Scalar>()); // ToDo: Replace with MOOCHO!
        stateSolver->defaultTol(fwd_newton_tol_);
        stateSolver->defaultMaxNewtonIterations(fwd_newton_max_iters_);
        stateSolver->useDampenedLineSearch(fwd_newton_dampening_);
        stateSolver->maxLineSearchIterations(fwd_newton_max_ls_iters_);
        // Create the reduced Thyra::ModelEvaluator object for p -> g_hat(p)
        Teuchos::RCP<Thyra::DefaultStateEliminationModelEvaluator<Scalar> >
          reducedThyraModel = rcp(new Thyra::DefaultStateEliminationModelEvaluator<Scalar>(outerModel_,stateSolver));
        Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >
          finalReducedThyraModel;
        if(use_finite_diff_for_obj_) {
          // Create the finite-difference wrapped Thyra::ModelEvaluator object
          Teuchos::RCP<Thyra::DefaultFiniteDifferenceModelEvaluator<Scalar> >
            fdReducedThyraModel = Thyra::defaultFiniteDifferenceModelEvaluator<Scalar>(
              reducedThyraModel, objDirecFiniteDiffCalculator
              );
          finalReducedThyraModel = fdReducedThyraModel;
        }
        else {
          finalReducedThyraModel = reducedThyraModel;
        }
        // Wrap the reduced NAND Thyra::ModelEvaluator object in an NLP object
        RCP<NLPFirstOrderThyraModelEvaluator>
          nlpFirstOrder = rcp(
            new NLPFirstOrderThyraModelEvaluator(finalReducedThyraModel,p_idx_,g_idx_)
            );
        nlpFirstOrder->showModelEvaluatorTrace(showModelEvaluatorTrace_);
        nlp = nlpFirstOrder;
      }
      else {
        switch(nlpType_) {
          case NLP_TYPE_DIRECT: {
            Teuchos::RCP<NLPDirectThyraModelEvaluator>
              nlpDirect = rcp(
                new NLPDirectThyraModelEvaluator(
                  outerModel_,p_idx_,g_idx_
                  ,objDirecFiniteDiffCalculator
                  ,conDirecFiniteDiffCalculator
                  )
                );
            nlpDirect->showModelEvaluatorTrace(showModelEvaluatorTrace_);
            nlp = nlpDirect;
            break;
          }
          case NLP_TYPE_FIRST_ORDER: {
            RCP<NLPFirstOrderThyraModelEvaluator>
              nlpFirstOrder = rcp(
                new NLPFirstOrderThyraModelEvaluator(outerModel_,p_idx_,g_idx_)
                );
            nlpFirstOrder->showModelEvaluatorTrace(showModelEvaluatorTrace_);
            nlp = nlpFirstOrder;
            break;
          }
          default:
            TEUCHOS_TEST_FOR_EXCEPT(true);
        }
      }
      break;
    }
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);
  }
    
  // Set the NLP
  solver_.set_nlp(nlp);

}

const Teuchos::RCP<Thyra::ModelEvaluator<value_type> >
MoochoThyraSolver::getOrigModel() const
{
  return origModel_;
}

const Teuchos::RCP<Thyra::ModelEvaluator<value_type> >
MoochoThyraSolver::getOuterModel() const
{
  return outerModel_;
}

void MoochoThyraSolver::readInitialGuess(
  std::ostream *out_arg
  )
{
  using Teuchos::OSTab;
  using Teuchos::RCP;
  using Thyra::clone;
  typedef Thyra::ModelEvaluatorBase MEB;

  RCP<Teuchos::FancyOStream>
    out = Teuchos::getFancyOStream(Teuchos::rcp(out_arg,false));
  RCP<MEB::InArgs<value_type> >
    initialGuess = clone(origModel_->getNominalValues()),
    lowerBounds = clone(origModel_->getLowerBounds()),
    upperBounds = clone(origModel_->getUpperBounds());

  RCP<const Thyra::VectorSpaceBase<value_type> >
    x_space = origModel_->get_x_space();
  if( 0 != x_space.get() ) {
    x_reader_.set_vecSpc(origModel_->get_x_space());
    if(this->get_stateVectorIO().get())
      x_reader_.set_fileIO(this->get_stateVectorIO());
    Teuchos::VerboseObjectTempState<Thyra::ParameterDrivenMultiVectorInput<value_type> >
      vots_x_reader(RCP<Thyra::ParameterDrivenMultiVectorInput<value_type> >(&x_reader_,false),out,Teuchos::VERB_LOW);
    initialGuess->set_x(
      readVectorOverride(
        x_reader_,"initial guess for the state \'x\'",initialGuess->get_x()
        )
      );
  }
  if( origModel_->Np() > 0 ) {
    p_reader_.set_vecSpc(origModel_->get_p_space(p_idx_));
    p_l_reader_.set_vecSpc(p_reader_.get_vecSpc());
    p_u_reader_.set_vecSpc(p_reader_.get_vecSpc());
    if(this->get_parameterVectorIO().get()) {
      p_reader_.set_fileIO(this->get_parameterVectorIO());
      p_l_reader_.set_fileIO(p_reader_.get_fileIO());
      p_u_reader_.set_fileIO(p_reader_.get_fileIO());
    }
    Teuchos::VerboseObjectTempState<Thyra::ParameterDrivenMultiVectorInput<value_type> >
      vots_p_reader(RCP<Thyra::ParameterDrivenMultiVectorInput<value_type> >(&p_reader_,false),out,Teuchos::VERB_LOW);
    initialGuess->set_p(
      p_idx_,
      readVectorOverride(
        p_reader_,"initial guess for the parameters \'p\'",
        initialGuess->get_p(p_idx_)
        )
      );
    lowerBounds->set_p(
      p_idx_,
      readVectorOverride(
        p_l_reader_,"lower bounds for the parameters \'p\'",lowerBounds->get_p(p_idx_)
        )
      );
    upperBounds->set_p(
      p_idx_,
      readVectorOverride(
        p_u_reader_,"upper bounds for the parameters \'p\'",upperBounds->get_p(p_idx_)
        )
      );
  }
  nominalModel_->setNominalValues(initialGuess);
  nominalModel_->setLowerBounds(lowerBounds);
  nominalModel_->setUpperBounds(upperBounds);
}

void MoochoThyraSolver::setInitialGuess(
  const Teuchos::RCP<const Thyra::ModelEvaluatorBase::InArgs<value_type> > &initialGuess
  )
{
  nominalModel_->setNominalValues(initialGuess);
}

void MoochoThyraSolver::setInitialGuess(
  const Thyra::ModelEvaluatorBase::InArgs<value_type> &initialGuess
  )
{
  nominalModel_->setNominalValues(
    Teuchos::rcp(new Thyra::ModelEvaluatorBase::InArgs<value_type>(initialGuess))
    );
}
  
MoochoSolver::ESolutionStatus MoochoThyraSolver::solve()
{
  using Teuchos::RCP; using Teuchos::null; using Teuchos::describe;
  solver_.update_solver();
  std::ostringstream os;
  os
    << "\n**********************************"
    << "\n*** MoochoThyraSolver::solve() ***"
    << "\n**********************************\n";
  const RCP<const Thyra::VectorSpaceBase<value_type> >
    x_space = outerModel_->get_x_space(),
    p_space = (
      ( p_idx_ >= 0 && outerModel_->Np() > 0 )
      ? outerModel_->get_p_space(p_idx_)
      : null
      );
  if( !is_null(x_space) )
    os << "\nx_space: " << describe(*x_space,Teuchos::VERB_MEDIUM);
  if( !is_null(p_space) )
    os << "\np_space: " << describe(*p_space,Teuchos::VERB_MEDIUM);
  *solver_.get_console_out() << os.str();
  *solver_.get_summary_out() << os.str();
  *solver_.get_journal_out() << os.str();
  outerModel_->setOStream(Teuchos::getFancyOStream(solver_.get_journal_out()));
  return solver_.solve_nlp();
}

const Thyra::ModelEvaluatorBase::InArgs<value_type>&
MoochoThyraSolver::getFinalPoint() const
{
  return finalPointModel_->getFinalPoint();
}

void MoochoThyraSolver::writeFinalSolution(
  std::ostream *out_arg
  ) const
{
  using Teuchos::OSTab;
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::getFancyOStream(Teuchos::rcp(out_arg,false));
  if( stateSoluFileBase_ != "" && finalPointModel_->getFinalPoint().get_x().get() ) {
    if(out.get())
      *out << "\nWriting the state solution \'x\' to the file(s) with base name \""<<stateSoluFileBase_<<"\" ...\n";
    stateVectorIO().writeMultiVectorToFile(
      *finalPointModel_->getFinalPoint().get_x(),stateSoluFileBase_
      );
  }
  if(
    ( "" != paramSoluFileBase_ )
    && ( origModel_->Np() > 0 )
    && ( 0 != finalPointModel_->getFinalPoint().get_p(p_idx_).get() )
    )
  {
    if(out.get())
      *out << "\nWriting the parameter solution \'p\' to the file(s) with base name \""<<paramSoluFileBase_<<"\" ...\n";
    parameterVectorIO().writeMultiVectorToFile(
      *finalPointModel_->getFinalPoint().get_p(p_idx_),paramSoluFileBase_
      );
  }
}

void MoochoThyraSolver::writeParamsFile(
  const std::string &outputXmlFileName
  ) const
{
  std::string xmlOutputFile
    = ( outputXmlFileName.length() ? outputXmlFileName : paramsUsedXmlOutFileName() );
  if( paramList_.get() && xmlOutputFile.length() ) {
    Teuchos::writeParameterListToXmlFile(*paramList_,xmlOutputFile);
  }
}

} // namespace MoochoPack
