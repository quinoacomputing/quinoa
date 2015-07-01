/* @HEADER@ */
// ************************************************************************
//
//                             Sundance
//                 Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Kevin Long (kevin.long@ttu.edu)
//

/* @HEADER@ */




#include "Sundance.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "PlayaMPIComm.hpp"
#include "PlayaMPISession.hpp"
#include "SundanceAssembler.hpp"
#include "SundanceQuadratureIntegral.hpp"
#include "SundanceQuadratureEvalMediator.hpp"
#include "SundanceDiscreteFuncEvaluator.hpp"
#include "SundanceRefIntegral.hpp"
#include "SundancePathUtils.hpp"
#include "SundanceEvaluator.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceCellJacobianBatch.hpp"
#include "SundanceUnaryFunctor.hpp"
#include "SundanceGrouperBase.hpp"
#include "SundanceSparsitySuperset.hpp"
#include "SundanceDefaultPath.hpp"
#include "SundanceDOFMapBuilder.hpp"
#include "SundanceVersionString.hpp"
#include "SundanceProductTransformation.hpp"
#include <unistd.h>
#ifndef _MSC_VER
#include <sys/unistd.h>
#else
#include "winprocess.h"
#endif
#ifdef HAVE_MPI
#include "mpi.h"
#endif

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#endif

namespace Sundance
{
  using Playa::MPIDataType;
  using Playa::MPIOp;
  using std::ostringstream;


void handleException(std::exception& e)
{SundanceGlobal::handleException(e);}


bool passFailTest(bool pass)
{return SundanceGlobal::passFailTest(pass);}

bool passFailTest(double error, double tol)
{return SundanceGlobal::passFailTest(error, tol);}


bool passFailTest(const std::string& statusMsg,
  bool status, double error, double tol)
{return SundanceGlobal::passFailTest(statusMsg, status, error, tol);}

int& testStatus()
{return SundanceGlobal::testStatus();}


CommandLineProcessor& clp()
{return SundanceGlobal::clp();}


int init(int* argc, char*** argv)
{return SundanceGlobal::init(argc, argv);}


int finalize() {return SundanceGlobal::finalize();}



void setOption(const std::string& optionName,
  int& value,
  const std::string& helpMsg)
{
  SundanceGlobal::setOption(optionName, value, helpMsg);
}

void setOption(const std::string& optionName,
  double& value,
  const std::string& helpMsg)
{
  SundanceGlobal::setOption(optionName, value, helpMsg);
}

void setOption(const std::string& optionName,
  std::string& value,
  const std::string& helpMsg)
{
  SundanceGlobal::setOption(optionName, value, helpMsg);
}

void setOption(const std::string& optionTrueName,
  const std::string& optionFalseName,
  bool& value,
  const std::string& helpMsg)
{
  SundanceGlobal::setOption(optionTrueName,
    optionFalseName,
    value,
    helpMsg);
}
} // namespace Sundance

static Time& totalTimer()
{
  static RCP<Time> rtn
    = TimeMonitor::getNewTimer("total Sundance time");
  return *rtn;
}

int SundanceGlobal::init(int* argc, char*** argv)
{

  try
  {
    /* start up MPI. In a serial run, this will be a no-op */
    //      MPISession::init(argc, argv);
    MPISession::init(argc, (char***) argv);

    /* Start a stopwatch. It will be stopped upon a call to finalize() */
    totalTimer().start();

    Tabs tab;

    /* read standard command line flags */
    std::string configFilename = "";

    bool defaultFpCheck = false;
    bool debugWait = false;
    bool showVersion = false;
    bool showBanner = true;
    bool showTimings = false;
    bool cmdFpCheck = defaultFpCheck;
    int defaultWorkSetSize = 400;
    int cmdWorkSetSize = defaultWorkSetSize;

    Assembler::workSetSize() = defaultWorkSetSize;

    clp().setOption("config", &configFilename, "Configuration file");
    clp().setOption("fpcheck", "nofpcheck", &cmdFpCheck,
      "Check results of math lib calls in expr evals");
    clp().setOption("version", "noversion", &showVersion,
      "Show Sundance version number and exit");
    clp().setOption("banner", "nobanner", &showBanner,
      "Show Sundance banner on root processor at start of run");
    clp().setOption("timings", "notimings", &showTimings,
      "Show timings at end of run");

    clp().setOption("workset", &cmdWorkSetSize,
      "Work set size");


    clp().setOption("debug", "nodebug", &debugWait, "Whether to attach a debugger to this process, holding until 'wait' is set to 0");


    clp().throwExceptions(false);
    clp().recogniseAllOptions(false);

    CommandLineProcessor::EParseCommandLineReturn rtn
      = clp().parse(*argc, (char**) *argv);

    TEUCHOS_TEST_FOR_EXCEPTION(rtn != CommandLineProcessor::PARSE_SUCCESSFUL,
      std::runtime_error,
      "Command-line parsing failed");


    if (showVersion)
    {
      if (MPIComm::world().getRank()==0)
      {
        cout << "Simulation built using Sundance version "
             << VersionString::number()
             << " (" << VersionString::date() << ")" << std::endl;
        exit(0);
      }
    }
    if (showBanner && MPIComm::world().getRank()==0)
    {
      ostringstream oss;
      oss << "Simulation built using Sundance version "
           << VersionString::number()
           << " (" << VersionString::date() << ")" << std::endl;

      oss << "Sundance is copyright"
           << std::endl << " (C) 2005-2012 Sandia National Laboratories "
           << std::endl
           << " (C) 2007-2012 Texas Tech University"
           << std::endl;
      oss << "and is licensed under the BSD License" << std::endl;
      oss << std::endl;
      cout << oss.str() << std::flush;
    }

    MPIComm::world().synchronize();
    if (!showTimings) skipTimingOutput() = true;

    //      debugWait = true;
    if (debugWait)
    {
      int wait=1;
      int pid = getpid();
      std::string myCommandName=((char**)(*argv))[0];
      std::string debugCmd = "ddd --gdb -x ~/.gdbinit " + myCommandName
        + " " + Teuchos::toString(pid) + " &";
      cout << "launching " << debugCmd << std::endl;
      TEUCHOS_ASSERT_EQUALITY(0, system(debugCmd.c_str()));
      while (wait) {;}
    }


    bool worksetSetOnCmdLine = cmdWorkSetSize != defaultWorkSetSize;
    if (worksetSetOnCmdLine)
    {
      Assembler::workSetSize() = (int) cmdWorkSetSize;
    }
  }
  catch(std::exception& e)
  {
    handleException(e);
    return 1;
  }
  return 0;
}


bool& SundanceGlobal::showStartupMessage()
{
  return MPISession::showStartupMessage();
}



void SundanceGlobal::handleException(std::exception& e)
{
  cout << "Sundance detected exception: " << std::endl;
  cout << e.what() << std::endl;
  cout << "test FAILED" << std::endl;
  testStatus() = -1;
}


int SundanceGlobal::finalize()
{
  totalTimer().stop();

  try
  {
    Tabs tab;
    /* we may need to skip timing summaries because of a Trilinos 6.0.x bug */
    if (!skipTimingOutput()) TimeMonitor::summarize();
    //  MPISession::finalize();
  }
  catch(std::exception& e)
  {
    handleException(e);
    return 1;
  }
  return 0;
}






bool SundanceGlobal::checkTest(double error, double tol)
{
  int myFail = error > tol;
  int anyFail = myFail;
  MPIComm::world().allReduce((void*) &myFail, (void*) &anyFail, 1, MPIDataType::intType(),
    MPIOp::sumOp());
  return (anyFail == 0);
}

bool SundanceGlobal:: passFailTest(double error, double tol)
{
  MPIComm::world().synchronize();
  bool pass;
  if (MPIComm::world().getRank()==0)
  {
    cout << "error norm = " << error << std::endl;
    cout << "tolerance = " << tol << std::endl;
  }
  pass = checkTest(error, tol);
  if (MPIComm::world().getRank()==0)
  {
    if (pass)
    {
      cout << "test PASSED" << std::endl;
    }
    else
    {
      cout << "test FAILED" << std::endl;
    }
  }
  testStatus() = pass!=true;
  return pass;
}


bool SundanceGlobal:: passFailTest(const std::string& statusMsg,
  bool status, double error, double tol)
{
  MPIComm::world().synchronize();
  bool pass;
  if (MPIComm::world().getRank()==0)
  {

    cout << statusMsg << ": ";
    if (status) cout << "true" << std::endl;
    else cout << "false" << std::endl;
    cout << "error norm = " << error << std::endl;
    cout << "tolerance = " << tol << std::endl;
  }
  pass = checkTest(error, tol);
  if (MPIComm::world().getRank()==0)
  {
    if (status && pass)
    {
      cout << "test PASSED" << std::endl;
    }
    else
    {
      cout << "test FAILED" << std::endl;
    }
  }
  testStatus() = pass!=true;
  return pass;
}

bool SundanceGlobal:: passFailTest(bool pass)
{
  MPIComm::world().synchronize();
  if (MPIComm::world().getRank()==0)
  {
    if (pass)
    {
      cout << "test PASSED" << std::endl;
    }
    else
    {
      cout << "test FAILED" << std::endl;
    }
  }
  testStatus() = pass!=true;
  return pass;
}


void SundanceGlobal::setOption(const std::string& optionName,
  int& value,
  const std::string& helpMsg)
{
  clp().setOption(optionName.c_str(), &value, helpMsg.c_str());
}

void SundanceGlobal::setOption(const std::string& optionName,
  double& value,
  const std::string& helpMsg)
{
  clp().setOption(optionName.c_str(), &value, helpMsg.c_str());
}

void SundanceGlobal::setOption(const std::string& optionName,
  std::string& value,
  const std::string& helpMsg)
{
  clp().setOption(optionName.c_str(), &value, helpMsg.c_str());
}

void SundanceGlobal::setOption(const std::string& optionTrueName,
  const std::string& optionFalseName,
  bool& value,
  const std::string& helpMsg)
{
  clp().setOption(optionTrueName.c_str(),
    optionFalseName.c_str(),
    &value,
    helpMsg.c_str());
}


