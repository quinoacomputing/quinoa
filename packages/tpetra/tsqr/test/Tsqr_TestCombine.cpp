//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

#include <Tsqr_ConfigDefs.hpp>
#include "Teuchos_ConfigDefs.hpp" // HAVE_MPI
#include "Teuchos_Tuple.hpp"
#ifdef HAVE_MPI
#  include "Teuchos_GlobalMPISession.hpp"
#  include "Teuchos_oblackholestream.hpp"
#endif // HAVE_MPI
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Tsqr_CombineBenchmark.hpp"
#include "Tsqr_CombineTest.hpp"

#ifdef HAVE_KOKKOSTSQR_COMPLEX
#  include <complex>
#endif // HAVE_KOKKOSTSQR_COMPLEX

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>


namespace {
  using Teuchos::RCP;

  //
  // Short description of this program, printed if the --help option
  // is given at the command line.
  //
  static const char docString[] =
    "This program tests accuracy and performance of TSQR::Combine.";

  //
  // The TestParameters struct encapsulates values of command-line
  // parameters.
  //
  struct TestParameters {
    TestParameters () :
      verify (false),
      benchmark (false),
      numRows (100),
      numCols (5),
      numTrials (3),
      calibrate (false),
      averageTimings (true),
      testReal (true),
#ifdef HAVE_KOKKOSTSQR_COMPLEX
      testComplex (true),
#endif // HAVE_KOKKOSTSQR_COMPLEX
      printFieldNames (true),
      printTrilinosTestStuff (true),
      strictPerfTests (false),
      allowance (1.2),
      verbose (true),
      debug (false)
      {}

    // Whether to run the accuracy test.
    bool verify;
    // Whether to run the performance test.
    bool benchmark;
    // Number of rows in the test matrix.
    int numRows;
    // Number of columns in the test matrix.
    int numCols;
    // Number of trials (benchmark only).
    int numTrials;
    // Whether to pick the number of trials automatically, using an
    // iterative calibration process (benchmark only).
    bool calibrate;
    // Whether to print averaged timings over all trials (true), or the cumulative timing over all trials (false).
    bool averageTimings;
    // Whether to test real-arithmetic routines.
    bool testReal;
#ifdef HAVE_KOKKOSTSQR_COMPLEX
    // Whether to test complex-arithmetic routines.  We don't let this
    // option exist unless TSQR was built with complex arithmetic
    // support.
    bool testComplex;
#endif // HAVE_KOKKOSTSQR_COMPLEX
    // Whether to print column (field) names.
    bool printFieldNames;
    // Whether to print output that the Trilinos test framework
    // expects, in order to judge a test as passed or failed.
    bool printTrilinosTestStuff;
    // Whether the benchmark should fail if performance of
    // TSQR::CombineNative (and TSQR::CombineFortran, if applicable)
    // relative to that of TSQR::CombineDefault is not good enough.
    bool strictPerfTests;
    // If strictPerfTests is true: how much slower CombineNative (and
    // CombineFortran, if applicable) is allowed to be, relative to
    // CombineDefault.
    double allowance;
    // Whether to print verbose status output.
    bool verbose;
    // Whether to print debugging output to stderr.
    bool debug;
    std::string additionalFieldNames, additionalData;
  };

  // Benchmark TSQR::Combine.
  //
  // out [out] output stream for benchmark results.
  //   It will only be used on rank 0.
  //
  // params [in] test parameter struct.  This method reads
  //   the following fields: numRows, numCols, numTrials,
  //   testReal, testComplex.
  //
  // Warning: Call only on (MPI) rank 0.  Otherwise, you'll run the
  //   test routine on every MPI rank simultaneously, but only report
  //   results on rank 0.
  void
    benchmark (std::ostream& out,
        const TestParameters& params)
    {
      std::vector<int> seed(4);
      const bool useSeedValues = false; // Fill in seed with defaults.

      using TSQR::Test::benchmarkCombine;
      typedef Teuchos::Time timer_type;

      TSQR::Test::CombineBenchmarkParameters testParams;
      testParams.numRows = params.numRows;
      testParams.numCols = params.numCols;
      testParams.testReal = params.testReal;
#ifdef HAVE_KOKKOSTSQR_COMPLEX
      testParams.testComplex = params.testComplex;
#else
      testParams.testComplex = false;
#endif // HAVE_KOKKOSTSQR_COMPLEX
      testParams.numTrials = params.numTrials;
      testParams.calibrate = params.calibrate;
      testParams.averageTimings = params.averageTimings;
      testParams.strictPerfTests = params.strictPerfTests;
      testParams.allowance = params.allowance;
      testParams.seed = seed;
      testParams.useSeedValues = useSeedValues;
      testParams.additionalFieldNames = params.additionalFieldNames;
      testParams.additionalData = params.additionalData;
      testParams.printFieldNames = params.printFieldNames;
      testParams.debug = params.debug;

      benchmarkCombine<timer_type> (out, testParams);
    }

  // Test accuracy of TSQR::Combine.
  //
  // out [out] output stream for benchmark results.
  //   It will only be used on rank 0.
  //
  // params [in] test parameter struct.  This method reads
  //   the following fields: numRows, numCols, numTrials,
  //   testReal, testComplex.
  //
  // Warning: Call only on (MPI) rank 0.  Otherwise, you'll run
  //   the test routine on every MPI rank simultaneously, but
  //   only report results on rank 0.
  void
    verify (std::ostream& out,
        const TestParameters& params)
    {
      typedef int ordinal_type;

      const ordinal_type numRows = params.numRows;
      const ordinal_type numCols = params.numCols;
#ifdef HAVE_KOKKOSTSQR_COMPLEX
      const bool testComplex = params.testComplex;
#else
      const bool testComplex = false;
#endif // HAVE_KOKKOSTSQR_COMPLEX
      const bool printFieldNames = params.printFieldNames;
      const bool simulateSequentialTsqr = false;
      const bool debug = false;

      using TSQR::Test::verifyCombine;
      verifyCombine (numRows, numCols, params.testReal, testComplex,
          printFieldNames, simulateSequentialTsqr, debug);
    }

  // \brief Parse command-line options for this test
  //
  // argc [in] As usual in C(++).
  // argv [in] As usual in C(++).
  //
  // allowedToPrint [in] Whether this (MPI) process is allowed
  //   to print to stdout/stderr.  Different per (MPI) process.
  //
  // printedHelp [out] Whether this (MPI) process printed the
  //   "help" display (summary of command-line options)
  //
  // Return: Encapsulation of command-line options.
  TestParameters
    parseOptions (int argc,
        char* argv[],
        const bool allowedToPrint,
        bool& printedHelp)
    {
      using std::cerr;
      using std::endl;

      printedHelp = false;

      // Command-line parameters, set to their default values.
      TestParameters params;
      try {
        using Teuchos::CommandLineProcessor;

        CommandLineProcessor cmdLineProc (/* throwExceptions=*/ true,
            /* recognizeAllOptions=*/ true);
        cmdLineProc.setDocString (docString);
        cmdLineProc.setOption ("verify",
            "noverify",
            &params.verify,
            "Test accuracy of TSQR::Combine implementations.");
        cmdLineProc.setOption ("benchmark",
            "nobenchmark",
            &params.benchmark,
            "Test performance of TSQR::Combine implementations.");
        cmdLineProc.setOption ("debug",
            "nodebug",
            &params.debug,
            "Print copious debugging information to stderr.");
        cmdLineProc.setOption ("numRows",
            &params.numRows,
            "Number of rows in the cache block test.");
        cmdLineProc.setOption ("numCols",
            &params.numCols,
            "Number of columns in the cache block test, and "
            "number of rows and columns in each upper triangular "
            "matrix in the pair test.");
        cmdLineProc.setOption ("numTrials",
            &params.numTrials,
            "For benchmarks: Number of trials.  "
            "Ignored if --calibrate option is set.");
        cmdLineProc.setOption ("calibrate",
            "noCalibrate",
            &params.calibrate,
            "For benchmarks: ignore numTrials, and calibrate "
            "the number of trials based on computed timer "
            "resolution and problem size (numRows and "
            "numCols).");
        cmdLineProc.setOption ("meanTimings",
            "sumTimings",
            &params.averageTimings,
            "For benchmarks: whether timings should be "
            "computed as an arithmetic mean (true) or as a "
            "sum (false) over all trials.");
        cmdLineProc.setOption ("testReal",
            "noTestReal",
            &params.testReal,
            "Test real-arithmetic routines.");
#ifdef HAVE_KOKKOSTSQR_COMPLEX
        cmdLineProc.setOption ("testComplex",
            "noTestComplex",
            &params.testComplex,
            "Test complex-arithmetic routines.  This option "
            "may only be set if Trilinos was built with "
            "complex arithmetic support.");
#endif // HAVE_KOKKOSTSQR_COMPLEX
        cmdLineProc.setOption ("strictPerfTests",
            "noStrictPerfTests",
            &params.strictPerfTests,
            "For benchmarks: whether the test should fail if "
            "run time of TSQR::CombineNative / run time of "
            "TSQR::CombineDefault (both for the cache block "
            "benchmark) is greater than the given slowdown "
            "allowance.  Ditto for TSQR::CombineFortran, if "
            "TSQR was built with Fortran support.");
        cmdLineProc.setOption ("allowance",
            &params.allowance,
            "For benchmarks: if strictPerfTests is true: "
            "allowed slowdown factor.  If exceeded, the test "
            "fails.");
        cmdLineProc.setOption ("additionalFieldNames",
            &params.additionalFieldNames,
            "Any additional field name(s) (comma-delimited "
            "string) to add to the benchmark output.  Empty "
            "by default.  Good for things known when invoking "
            "the benchmark executable, but not (easily) known "
            "inside the benchmark -- e.g., environment "
            "variables.");
        cmdLineProc.setOption ("additionalData",
            &params.additionalData,
            "Any additional data to add to the output, "
            "corresponding to the above field name(s). "
            "Empty by default.");
        cmdLineProc.setOption ("printFieldNames",
            "noPrintFieldNames",
            &params.printFieldNames,
            "Print field names for benchmark output (including "
            "any arguments to --fieldNames).");
        cmdLineProc.setOption ("printTrilinosTestStuff",
            "noPrintTrilinosTestStuff",
            &params.printTrilinosTestStuff,
            "Print output that makes the Trilinos test "
            "framework happy (but makes benchmark results "
            "parsing scripts unhappy)");
        cmdLineProc.parse (argc, argv);
      }
      catch (Teuchos::CommandLineProcessor::UnrecognizedOption& e) {
        if (allowedToPrint)
          cerr << "Unrecognized command-line option: " << e.what() << endl;
        throw e;
      }
      catch (Teuchos::CommandLineProcessor::HelpPrinted& e) {
        printedHelp = true;
        return params; // Don't verify parameters in this case
      }

      // Validate.  TODO (mfh 08 Jul 2010) Figure out how to do this with
      // ParameterList validators.
      if (params.numRows <= 0)
        throw std::invalid_argument ("Number of rows must be positive");
      else if (params.numCols <= 0)
        throw std::invalid_argument ("Number of columns must be positive");
      else if (params.numRows < params.numCols)
        throw std::invalid_argument ("Number of rows must be >= number of columns");
      else if (params.benchmark && params.numTrials < 1)
        throw std::invalid_argument ("Benchmark requires numTrials >= 1");

      return params;
    }
} // namespace (anonymous)


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

  int
main (int argc, char *argv[])
{
  using Teuchos::RCP;

#ifdef HAVE_MPI
  typedef RCP< const Teuchos::Comm<int> > comm_ptr;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackhole);
  comm_ptr comm = Teuchos::DefaultComm<int>::getComm();
  const int myRank = comm->getRank();
  // Only Rank 0 gets to write to stdout.  The other MPI process ranks
  // send their output to something that looks like /dev/null (and
  // likely is, on Unix-y operating systems).
  std::ostream& out = (myRank == 0) ? std::cout : blackhole;
  // Only Rank 0 performs the tests.
  const bool performingTests = (myRank == 0);
  const bool allowedToPrint = (myRank == 0);

#else // Don't HAVE_MPI: single-node test

  const bool performingTests = true;
  const bool allowedToPrint = true;
  std::ostream& out = std::cout;
#endif // HAVE_MPI

  // Fetch command-line parameters.
  bool printedHelp = false;
  TestParameters params =
    parseOptions (argc, argv, allowedToPrint, printedHelp);
  if (printedHelp)
    return 0;

  bool success = false;
  bool verbose = false;
  try {
    if (performingTests)
    {
      using std::endl;

      if (params.benchmark)
        benchmark (out, params);

      // We allow the same run to do both benchmark and verify.
      if (params.verify)
        verify (out, params);

      success = true;

      if (params.printTrilinosTestStuff)
        // The Trilinos test framework expects a message like this.
        out << "\nEnd Result: TEST PASSED" << endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);
  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
