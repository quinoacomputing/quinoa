/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER
*/

#include <Tpetra_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Map.hpp>

#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION
#include "Tpetra_Map_def.hpp"
#include "Tpetra_Directory_def.hpp"
#endif

namespace {

  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::Array;
  using Tpetra::global_size_t;
  using Teuchos::outArg;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
  }

  RCP<const Teuchos::Comm<int> > getDefaultComm()
  {
    return Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  }

  RCP<Tpetra::DefaultPlatform::DefaultPlatformType::NodeType> getDefaultNode()
  {
    return Tpetra::DefaultPlatform::getDefaultPlatform().getNode();
  }

  //
  // UNIT TESTS
  //

  ////
  TEUCHOS_UNIT_TEST( Map, Bug5401_NegativeBaseIndex )
  {
    using std::endl;
    out << "Bug 5401 (negative index base) test" << endl;
    Teuchos::OSTab tab0 (out);

    // failure reading 1x4 matrix under MPI
    typedef int                          LO;
    typedef int                          GO;
    typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType node_type;
    typedef Tpetra::Map<LO,GO,node_type> map_type;
    typedef Teuchos::Comm<int>           comm_type;
    typedef Teuchos::Array<GO>::size_type size_type;

    // create a comm
    RCP<const comm_type> comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myRank = comm->getRank ();

    out << "Comm has " << numImages << " process"
        << (numImages != 1 ? "es" : "") << endl;
    TEUCHOS_TEST_FOR_EXCEPTION(
      numImages != 2,
      std::logic_error,
      "This test is appropriate only for MPI runs of rank 2.")
    RCP<node_type> node = getDefaultNode ();

    const GO numElements = 78;
    const GO baseIndexIsNegOne = -1;
    const global_size_t GINV   = Teuchos::OrdinalTraits<global_size_t>::invalid();
    Array<int> elements (numElements);

    out << "Create array of global indices.  All processes have the same "
        << "global index.  The first global index on all processes is "
        << baseIndexIsNegOne << "." << endl;

    // first global element is -1
    for (size_type i = 0; i < elements.size (); ++i) {
      elements[i] = i - 1;
    }

    //int localMapCtorSuccess = 0;
    RCP<map_type> map;

    map = Teuchos::rcp (new map_type (GINV, elements(), baseIndexIsNegOne, comm));

    // try {
    //   map = Teuchos::rcp (new map_type (GINV, elements(), baseIndexIsNegOne, comm));
    //   int localMapCtorSuccess = 1;
    // } catch (std::exception& e) {
    //   out << "Process " << myRank << ": Noncontiguous Map constructor failed "
    //       << "with the following exception message: " << e.what () << endl;
    // }
    // int globalMapCtorSuccess = 0;
    // Teuchos::reduceAll<int, int> (*comm, Teuchos::REDUCE_MIN,
    //                               localMapCtorSuccess,
    //                               outArg (globalMapCtorSuccess));
    // TEST_EQUALITY_CONST( globalMapCtorSuccess, 1 );

    out << "Process " << myRank << ":";
    {
      Teuchos::OSTab tab1 (out);
      out << "My number of global indices: " << map->getNodeNumElements () << endl
          << "Global number of global indices: " << map->getGlobalNumElements () << endl
          << "Index base: " << map->getIndexBase () << endl
          << "My min global index: " << map->getMinGlobalIndex () << endl
          << "Global min global index: " << map->getMinAllGlobalIndex () << endl;
    }

    TEST_EQUALITY( map->getNodeNumElements(),   as<size_t> (numElements) );
    TEST_EQUALITY( map->getGlobalNumElements(), as<global_size_t> (numElements*numImages) );
    TEST_EQUALITY( map->getIndexBase(),         as<GO> (-1) );
    TEST_EQUALITY( map->getMinGlobalIndex(),    as<GO> (-1) );
    TEST_EQUALITY( map->getMinAllGlobalIndex(), as<GO> (-1) );

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

}


