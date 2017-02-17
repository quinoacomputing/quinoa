// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include <Teuchos_MatrixMarket_Raw_Checker.hpp>
#include <Teuchos_MatrixMarket_Raw_Reader.hpp>
#include <Teuchos_MatrixMarket_SetScientific.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_DefaultSerialComm.hpp>
#include <algorithm>

using std::endl;

namespace {
  // Sample Matrix Market sparse matrix file.  We include this so we
  // can test without needing to read in a file.  Notice that all the
  // decimal floating-point values in this example can be represented
  // exactly in binary floating point.  This example has correct
  // syntax, so you won't need to use tolerant mode to parse it.
  const char sampleMatrixMarketFile[] =
    "%%MatrixMarket matrix coordinate real general\n"
    "5 5 10\n"
    "5 5 55.0\n"
    "4 4 44.0\n"
    "3 3 33.0\n"
    "2 2 22.0\n"
    "1 1 11.0\n"
    "4 5 45.0\n"
    "3 4 34.0\n"
    "2 3 23.0\n"
    "1 2 12.0\n"
    "1 5 15.0\n";

  // Given the three arrays of a CSR data structure, along with the
  // numbers of rows and columns, print the result to the given output
  // stream as a MatrixMarket file.
  template<class OrdinalType, class ScalarType>
  void
  csrToMatrixMarket (std::ostream& out,
                     Teuchos::ArrayView<OrdinalType> ptr,
                     Teuchos::ArrayView<OrdinalType> ind,
                     Teuchos::ArrayView<ScalarType> val,
                     const OrdinalType numRows,
                     const OrdinalType numCols)
  {
    using Teuchos::ArrayView;
    using std::endl;
    typedef typename ArrayView<OrdinalType>::size_type size_type;
    typedef Teuchos::ScalarTraits<ScalarType> STS;

    // Make the output stream write floating-point numbers in
    // scientific notation.  It will politely put the output
    // stream back to its state on input, when this scope
    // terminates.
    Teuchos::MatrixMarket::details::SetScientific<ScalarType> sci (out);

    out << "%%MatrixMarket matrix coordinate ";
    if (STS::isComplex) {
      out << "complex ";
    }
    else {
      out << "real ";
    }
    out << "general" << endl;
    out << numRows << " " << numCols << " " << ptr[numRows] << endl;
    OrdinalType k;
    for (OrdinalType i = 0; i < numRows; ++i) {
      for (k = ptr[i]; k < ptr[i+1]; ++k) {
        // Matrix Market files use 1-based row and column indices.
        out << (i+1) << " " << (ind[k]+1) << " ";
        if (STS::isComplex) {
          out << STS::real (val[k]) << " " << STS::imag (val[k]);
        }
        else {
          out << val[k];
        }
        out << endl;
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION(k != ptr[numRows], std::logic_error,
      "csrToMatrixMarket: Failed to print all the matrix entries!  The last k "
      "index value is " << k << ", but the number of entries is " << ptr[numRows]
      << ".");
  }
} // namespace (anonymous)

// Benchmark driver
int
main (int argc, char *argv[])
{
  using Teuchos::MatrixMarket::Raw::Checker;
  using Teuchos::MatrixMarket::Raw::Reader;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::Comm;
  using Teuchos::CommandLineProcessor;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Teuchos::SerialComm;
  using std::cout;
  using std::cerr;
  typedef double scalar_type;
  typedef int ordinal_type;

  // Name of the Matrix Market sparse matrix file to read.  If empty,
  // use the Matrix Market example embedded as a string in this file.
  std::string filename;
  // If true, just check the sparse matrix file.  Otherwise,
  // do a full conversion to CSR (compressed sparse row) format.
  bool checkOnly = false;
  // Whether to echo the sparse matrix to stdout after reading it
  // successfully.
  bool echo = false;
  // Whether to parse the Matrix Market file tolerantly.
  bool tolerant = false;
  // Verbosity of output
  bool verbose = false;
  // Whether to print debugging-level output
  bool debug = false;

  CommandLineProcessor cmdp (false, true);
  cmdp.setOption ("filename", &filename,
                  "Name of the Matrix Market sparse matrix file to read.");
  cmdp.setOption ("checkOnly", "fullTest", &checkOnly,
                  "If true, just check the syntax of the input file.  "
                  "Otherwise, do a full test.");
  cmdp.setOption ("echo", "noecho", &echo,
                  "Whether to echo the sparse matrix contents to stdout "
                  "after reading it successfully.");
  cmdp.setOption ("tolerant", "strict", &tolerant,
                  "Whether to tolerate syntax errors in the Matrix Market file.");
  cmdp.setOption ("verbose", "quiet", &verbose,
                  "Print status output to stdout.");
  cmdp.setOption ("debug", "nodebug", &debug,
                  "Print possibly copious debugging output to stderr.");
  // Parse the command-line arguments.
  {
    const CommandLineProcessor::EParseCommandLineReturn parseResult =
      cmdp.parse (argc,argv);
    // If the caller asks us to print the documentation, or does not
    // explicitly say to run the benchmark, we let this "test" pass
    // trivially.
    if (parseResult == CommandLineProcessor::PARSE_HELP_PRINTED) {
      std::cout << "End Result: TEST PASSED" << endl;
      return EXIT_SUCCESS;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
       parseResult != CommandLineProcessor::PARSE_SUCCESSFUL,
       std::invalid_argument, "Failed to parse command-line arguments.");
  }

  // Test reading in the sparse matrix.  If no filename or an empty
  // filename is specified, the test passes trivially.
  bool success = true;
  if (checkOnly) {
    typedef Checker<scalar_type, ordinal_type> checker_type;
    checker_type checker (echo, tolerant, debug);

    RCP<const Comm<int> > comm = rcp (new SerialComm<int>);
    if (filename != "") {
      if (verbose) {
        cout << "Checking syntax of the Matrix Market file \"" << filename
             << "\"" << endl;
      }
      success = success && checker.readFile (*comm, filename);
      if (verbose) {
        if (success) {
          cout << "The given file is a valid Matrix Market file." << endl;
        }
        else {
          cout << "The given file has syntax errors." << endl;
        }
      }
    }
    else {
      if (verbose) {
        cout << "Checking syntax of the Matrix Market string example" << endl
             << std::flush;// for debug output next
      }
      if (debug) {
        cerr << "Matrix Market string example: " << endl
             << sampleMatrixMarketFile << endl;
      }
      std::istringstream in (sampleMatrixMarketFile);
      RCP<std::istream> inStream = rcpFromRef (in);
      success = success && checker.read (*comm, inStream);
      if (verbose) {
        if (success) {
          cout << "The example has valid Matrix Market syntax." << endl;
        }
        else {
          cout << "The example has syntax errors." << endl;
        }
      }
    }
  }
  else {
    typedef Reader<scalar_type, ordinal_type> reader_type;
    reader_type reader (tolerant, debug);
    ArrayRCP<ordinal_type> ptr, ind;
    ArrayRCP<scalar_type> val;
    ordinal_type numRows, numCols;
    if (filename != "") {
      if (verbose) {
        cout << "Reading the Matrix Market file \"" << filename << "\"" << endl;
      }
      success = success && reader.readFile (ptr, ind, val,
                                            numRows, numCols, filename);
    }
    else {
      if (verbose) {
        cout << "Reading the Matrix Market string example" << endl;
      }
      if (debug) {
        cerr << "Matrix Market string example:" << endl
             << sampleMatrixMarketFile << endl;
      }
      std::istringstream inStr (sampleMatrixMarketFile);
      success = success && reader.read (ptr, ind, val, numRows, numCols, inStr);
    }
    TEUCHOS_TEST_FOR_EXCEPTION(! success, std::runtime_error, "Matrix Market "
      "reader failed to read the given file or input stream.");
    if (success && verbose) {
      cout << "Successfully read the Matrix Market string example" << endl
           << std::flush; // for following debug output
    }
    if (debug) {
      cerr << "CSR output info:" << endl
           << "  ptr.size() = " << ptr.size()
           << ", ind.size() = " << ind.size()
           << ", val.size() = " << val.size()
           << ", numRows = " << numRows
           << ", numCols = " << numCols
           << endl;
    }

    // Here's the fun part.  Output the CSR data to an output stream.
    // Then read in the output stream.  The resulting matrix should be
    // exactly the same (unless the original file had elements at the
    // same location that were added together with rounding error).
    std::ostringstream outStr;
    if (success && verbose) {
      cout << "Printing the CSR arrays to a Matrix Market output stream"
           << endl << std::flush;
    }
    csrToMatrixMarket<ordinal_type, scalar_type> (outStr, ptr (), ind (), val (),
                                                  numRows, numCols);
    if (debug && echo) {
      cerr << "CSR data:" << endl
           << "- ptr = [";
      for (ordinal_type i = 0; i < ptr.size(); ++i) {
        cerr << ptr[i];
        if (i+1 != ptr.size()) { // don't subtract from zero if unsigned
          cerr << ", ";
        }
      }
      cerr << "]" << endl
           << "- ind = [";
      for (ordinal_type i = 0; i < ind.size(); ++i) {
        cerr << ind[i];
        if (i+1 != ind.size()) { // don't subtract from zero if unsigned
          cerr << ", ";
        }
      }
      cerr << "]" << endl
           << "- val = [";
      for (ordinal_type i = 0; i < val.size(); ++i) {
        cerr << val[i];
        if (i+1 != val.size()) { // don't subtract from zero if unsigned
          cerr << ", ";
        }
      }
      cerr << "]" << endl;

      cerr << "CSR data, converted back to Matrix Market format" << endl;
      csrToMatrixMarket<ordinal_type, scalar_type> (cerr, ptr (), ind (), val (),
                                                    numRows, numCols);
      cerr << endl;
    }

    std::istringstream inStr (outStr.str ());
    ArrayRCP<ordinal_type> newptr, newind;
    ArrayRCP<scalar_type> newval;
    ordinal_type newNumRows, newNumCols;
    if (success && verbose) {
      cout << "Reading the Matrix Market output back into CSR arrays" << endl;
    }
    success = success && reader.read (newptr, newind, newval,
                                      newNumRows, newNumCols, inStr);
    TEUCHOS_TEST_FOR_EXCEPTION(! success, std::logic_error, "Matrix Market "
      "reader failed to read the output back into CSR arrays.");
    if (success && verbose) {
      cout << "Successfully read the Matrix Market output back into CSR arrays"
           << endl << std::flush;
    }
    if (debug) {
      cerr << "CSR output info:" << endl
           << "  newptr.size() = " << newptr.size()
           << ", newind.size() = " << newind.size()
           << ", newval.size() = " << newval.size()
           << ", newNumRows = " << newNumRows
           << ", newNumCols = " << newNumCols
           << endl;
    }

    // The old arrays should equal the new arrays.
    TEUCHOS_TEST_FOR_EXCEPTION(ptr.size () != newptr.size (), std::logic_error,
      "New ptr array has a different length than old ptr array");
    TEUCHOS_TEST_FOR_EXCEPTION(ind.size () != newind.size (), std::logic_error,
      "New ind array has a different length than old ind array");
    TEUCHOS_TEST_FOR_EXCEPTION(val.size () != newval.size (), std::logic_error,
      "New val array has a different length than old val array");
    TEUCHOS_TEST_FOR_EXCEPTION(newNumRows != numRows || newNumCols != numCols,
      std::logic_error, "New dimensions differ from old dimensions");
    TEUCHOS_TEST_FOR_EXCEPTION(ptr.size () != numRows+1, std::logic_error,
      "ptr.size() != numRows+1");
    TEUCHOS_TEST_FOR_EXCEPTION(newptr.size () != newNumRows+1, std::logic_error,
      "newptr.size() != newNumRows+1");

    for (ordinal_type rowIndex = 0; rowIndex < numRows; ++rowIndex) {
      TEUCHOS_TEST_FOR_EXCEPTION(ptr[rowIndex] != newptr[rowIndex],
        std::logic_error, "At row index " << rowIndex << ", ptr[rowIndex] = "
        << ptr[rowIndex] << " != newptr[rowIndex] = " << newptr[rowIndex]
        << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(ptr[rowIndex+1] != newptr[rowIndex+1],
        std::logic_error, "At row index " << rowIndex << ", ptr[rowIndex+1] = "
        << ptr[rowIndex+1] << " != newptr[rowIndex+1] = " << newptr[rowIndex+1]
        << ".");
      for (ordinal_type k = ptr[rowIndex]; k < ptr[rowIndex+1]; ++k) {
        TEUCHOS_TEST_FOR_EXCEPTION(ind[k] != newind[k], std::logic_error,
          "At row index " << rowIndex << ", ind[k=" << k << "] = "
          << ind[k] << " != newind[k] = " << newind[k] << ".");
        // You may want to relax this inequality if the original
        // Matrix Market file had multiple entries at the same
        // location and if adding them together resulted in rounding
        // error.
        TEUCHOS_TEST_FOR_EXCEPTION(val[k] != newval[k], std::logic_error,
          "At row index " << rowIndex << ", val[k=" << k << "] = "
          << val[k] << " != newval[k] = " << newval[k] << ".");
      }
    }
  }

  if (success) {
    std::cout << "End Result: TEST PASSED" << endl;
    return EXIT_SUCCESS;
  }
  else {
    std::cout << "End Result: TEST FAILED" << endl;
    return EXIT_FAILURE;
  }
}



