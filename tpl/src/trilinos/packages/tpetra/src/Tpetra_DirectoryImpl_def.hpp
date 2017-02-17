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

#ifndef __Tpetra_DirectoryImpl_def_hpp
#define __Tpetra_DirectoryImpl_def_hpp

#ifdef DOXYGEN_USE_ONLY
#  include <Tpetra_DirectoryImpl_decl.hpp>
#endif
#include <Tpetra_Distributor.hpp>

namespace Tpetra {
  namespace Details {
    template<class LO, class GO, class NT>
    Directory<LO, GO, NT>::
    Directory (const Teuchos::RCP<const typename Directory<LO, GO, NT>::map_type>& map) :
      map_ (map)
    {}

    template<class LO, class GO, class NT>
    LookupStatus
    Directory<LO, GO, NT>::
    getEntries (const Teuchos::ArrayView<const GO> &globalIDs,
                const Teuchos::ArrayView<int> &nodeIDs,
                const Teuchos::ArrayView<LO> &localIDs,
                const bool computeLIDs) const
    {
      // Ensure that globalIDs, nodeIDs, and localIDs (if applicable)
      // all have the same size, before modifying any output arguments.
      TEUCHOS_TEST_FOR_EXCEPTION(nodeIDs.size() != globalIDs.size(),
        std::invalid_argument, Teuchos::typeName(*this) << "::getEntries(): "
        "Output arrays do not have the right sizes.  nodeIDs.size() = "
        << nodeIDs.size() << " != globalIDs.size() = " << globalIDs.size()
        << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(
        computeLIDs && localIDs.size() != globalIDs.size(),
        std::invalid_argument, Teuchos::typeName(*this) << "::getEntries(): "
        "Output array do not have the right sizes.  localIDs.size() = "
        << localIDs.size() << " != globalIDs.size() = " << globalIDs.size()
        << ".");

      // Initially, fill nodeIDs and localIDs (if applicable) with
      // invalid values.  The "invalid" process ID is -1 (this means
      // the same thing as MPI_ANY_SOURCE to Teuchos, so it's an
      // "invalid" process ID); the invalid local ID comes from
      // OrdinalTraits.
      std::fill (nodeIDs.begin(), nodeIDs.end(), -1);
      if (computeLIDs) {
        std::fill (localIDs.begin(), localIDs.end(),
                   Teuchos::OrdinalTraits<LO>::invalid ());
      }
      // Actually do the work.
      return this->getEntriesImpl (globalIDs, nodeIDs, localIDs, computeLIDs);
    }

    template<class LO, class GO, class NT>
    ReplicatedDirectory<LO, GO, NT>::
    ReplicatedDirectory (const Teuchos::RCP<const typename ReplicatedDirectory<LO, GO, NT>::map_type>& map) :
      Directory<LO, GO, NT> (map)
    {}

    template<class LO, class GO, class NT>
    std::string
    ReplicatedDirectory<LO, GO, NT>::description () const
    {
      std::ostringstream os;
      os << "ReplicatedDirectory"
         << "<" << Teuchos::TypeNameTraits<LO>::name ()
         << ", " << Teuchos::TypeNameTraits<GO>::name ()
         << ", " << Teuchos::TypeNameTraits<NT>::name () << ">";
      return os.str ();
    }

    template<class LO, class GO, class NT>
    DistributedContiguousDirectory<LO, GO, NT>::
    DistributedContiguousDirectory (const Teuchos::RCP<const typename DistributedContiguousDirectory<LO, GO, NT>::map_type>& map) :
      Directory<LO, GO, NT> (map)
    {
      using Teuchos::gatherAll;
      using Teuchos::RCP;

      RCP<const Teuchos::Comm<int> > comm = map->getComm ();

      TEUCHOS_TEST_FOR_EXCEPTION(! map->isDistributed (), std::invalid_argument,
        Teuchos::typeName (*this) << " constructor: Map is not distributed.");
      TEUCHOS_TEST_FOR_EXCEPTION(! map->isContiguous (), std::invalid_argument,
        Teuchos::typeName (*this) << " constructor: Map is not contiguous.");

      // Make room for the min global ID on each proc, plus one
      // entry at the end for the max cap.
      allMinGIDs_.resize (comm->getSize () + 1);
      // Get my process' min global ID.
      GO minMyGID = map->getMinGlobalIndex ();
      // Gather all of the min global IDs into the first getSize()
      // entries of allMinGIDs_.
      Teuchos::gatherAll<int, GO> (*comm, 1, &minMyGID, comm->getSize (),
                                   allMinGIDs_.getRawPtr ());
      // Put the max cap at the end.  Adding one lets us write loops
      // over the global IDs with the usual strict less-than bound.
      allMinGIDs_.back() = map->getMaxAllGlobalIndex() +
        Teuchos::OrdinalTraits<GO>::one();
    }

    template<class LO, class GO, class NT>
    std::string
    DistributedContiguousDirectory<LO, GO, NT>::description () const
    {
      std::ostringstream os;
      os << "DistributedContiguousDirectory"
         << "<" << Teuchos::TypeNameTraits<LO>::name ()
         << ", " << Teuchos::TypeNameTraits<GO>::name ()
         << ", " << Teuchos::TypeNameTraits<NT>::name () << ">";
      return os.str ();
    }

    template<class LO, class GO, class NT>
    LookupStatus
    ReplicatedDirectory<LO, GO, NT>::
    getEntriesImpl (const Teuchos::ArrayView<const GO> &globalIDs,
                    const Teuchos::ArrayView<int> &nodeIDs,
                    const Teuchos::ArrayView<LO> &localIDs,
                    const bool computeLIDs) const
    {
      using Teuchos::Array;
      using Teuchos::ArrayRCP;
      using Teuchos::ArrayView;
      using Teuchos::as;
      using Teuchos::Comm;
      using Teuchos::RCP;

      LookupStatus res = AllIDsPresent;
      RCP<const map_type> map = this->getMap ();
      RCP<const Teuchos::Comm<int> > comm = map->getComm ();
      const int myRank = comm->getRank ();

      // Map is on one process or is locally replicated.
      typename ArrayView<int>::iterator procIter = nodeIDs.begin();
      typename ArrayView<LO>::iterator lidIter = localIDs.begin();
      typename ArrayView<const GO>::iterator gidIter;
      for (gidIter = globalIDs.begin(); gidIter != globalIDs.end(); ++gidIter) {
        if (map->isNodeGlobalElement (*gidIter)) {
          *procIter++ = myRank;
          if (computeLIDs) {
            *lidIter++ = map->getLocalElement (*gidIter);
          }
        }
        else {
          // Advance the pointers, leaving these values set to invalid
          procIter++;
          if (computeLIDs) {
            lidIter++;
          }
          res = IDNotPresent;
        }
      }
      return res;
    }

    template<class LO, class GO, class NT>
    LookupStatus
    DistributedContiguousDirectory<LO, GO, NT>::
    getEntriesImpl (const Teuchos::ArrayView<const GO> &globalIDs,
                    const Teuchos::ArrayView<int> &nodeIDs,
                    const Teuchos::ArrayView<LO> &localIDs,
                    const bool computeLIDs) const
    {
      using Teuchos::Array;
      using Teuchos::ArrayRCP;
      using Teuchos::ArrayView;
      using Teuchos::as;
      using Teuchos::Comm;
      using Teuchos::RCP;

      const LO LINVALID = Teuchos::OrdinalTraits<LO>::invalid();
      LookupStatus res = AllIDsPresent;

      RCP<const map_type> map  = this->getMap ();
      RCP<const Teuchos::Comm<int> > comm = map->getComm ();
      const int numProcs = comm->getSize ();
      const global_size_t nOverP = map->getGlobalNumElements () / numProcs;

      // Map is distributed but contiguous.
      typename ArrayView<int>::iterator procIter = nodeIDs.begin();
      typename ArrayView<LO>::iterator lidIter = localIDs.begin();
      typename ArrayView<const GO>::iterator gidIter;
      for (gidIter = globalIDs.begin(); gidIter != globalIDs.end(); ++gidIter) {
        LO LID = LINVALID; // Assume not found until proven otherwise
        int image = -1;
        GO GID = *gidIter;
        // Guess uniform distribution and start a little above it
        // TODO: replace by a binary search
        int curRank;
        { // We go through all this trouble to avoid overflow and
          // signed / unsigned casting mistakes (that were made in
          // previous versions of this code).
          const GO one = Teuchos::OrdinalTraits<GO>::one();
          const GO two = one + one;
          const GO nOverP_GID = static_cast<GO> (nOverP);
          const GO lowerBound = GID / std::max(nOverP_GID, one) + two;
          // It's probably not OK to cast this to int in general.  It
          // works as long as |GID| <= the global number of entries
          // and nOverP is appropriately sized for int.  Trouble may
          // ensue if the index base has an exotic value.
          const int lowerBound_int = static_cast<int> (lowerBound);
          curRank = std::min (lowerBound_int, numProcs - 1);
        }
        bool found = false;
        while (curRank >= 0 && curRank < numProcs) {
          if (allMinGIDs_[curRank] <= GID) {
            if (GID < allMinGIDs_[curRank + 1]) {
              found = true;
              break;
            }
            else {
              curRank++;
            }
          }
          else {
            curRank--;
          }
        }
        if (found) {
          image = curRank;
          LID = as<LO> (GID - allMinGIDs_[image]);
        }
        else {
          res = IDNotPresent;
        }
        *procIter++ = image;
        if (computeLIDs) {
          *lidIter++ = LID;
        }
      }
      return res;
    }

    template<class LO, class GO, class NT>
    DistributedNoncontiguousDirectory<LO, GO, NT>::
    DistributedNoncontiguousDirectory (const Teuchos::RCP<const typename DistributedNoncontiguousDirectory<LO, GO, NT>::map_type>& map) :
      Directory<LO, GO, NT> (map)
    {
      using Teuchos::Array;
      using Teuchos::ArrayView;
      using Teuchos::as;
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::typeName;
      using Teuchos::TypeNameTraits;

      // This class' implementation of getEntriesImpl() currently
      // encodes the following assumptions:
      //
      // 1. global_size_t >= GO
      // 2. global_size_t >= int
      // 3. global_size_t >= LO
      //
      // We check these assumptions here.
      TEUCHOS_TEST_FOR_EXCEPTION(sizeof(global_size_t) < sizeof(GO),
        std::logic_error, typeName (*this) << ": sizeof(Tpetra::"
        "global_size_t) = " << sizeof(global_size_t) << " < sizeof(Global"
        "Ordinal = " << TypeNameTraits<LO>::name () << ") = " << sizeof(GO)
        << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(sizeof(global_size_t) < sizeof(int),
        std::logic_error, typeName (*this) << ": sizeof(Tpetra::"
        "global_size_t) = " << sizeof(global_size_t) << " < sizeof(int) = "
        << sizeof(int) << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(sizeof(global_size_t) < sizeof(LO),
        std::logic_error, typeName (*this) << ": sizeof(Tpetra::"
        "global_size_t) = " << sizeof(global_size_t) << " < sizeof(Local"
        "Ordinal = " << TypeNameTraits<LO>::name () << ") = " << sizeof(LO)
        << ".");

      RCP<const Teuchos::Comm<int> > comm = map->getComm ();
      const LO LINVALID = Teuchos::OrdinalTraits<LO>::invalid ();
      const GO minAllGID = map->getMinAllGlobalIndex ();
      const GO maxAllGID = map->getMaxAllGlobalIndex ();

      // The "Directory Map" (see below) will have a range of elements
      // from the minimum to the maximum GID of the user Map, and a
      // minimum GID of minAllGID from the user Map.
      const global_size_t numGlobalEntries = maxAllGID - minAllGID + 1;

      // We can't afford to replicate the whole directory on each
      // process, so create the "Directory Map", a uniform contiguous
      // Map that describes how we will distribute the directory over
      // processes.
      //
      // FIXME (mfh 08 May 2012) Here we're setting minAllGID to be
      // the index base.  The index base should be separate from the
      // minimum GID.
      directoryMap_ = rcp (new map_type (numGlobalEntries, minAllGID, comm, GloballyDistributed, map->getNode ()));

      // The number of Directory elements that my process owns.
      const size_t dir_numMyEntries = directoryMap_->getNodeNumElements ();

      // NOTE (mfh 21 Mar 2012): I rephrased the comment below so that
      // it made a little more sense, but I don't fully understand what
      // it means yet.
      //
      // Allocate process ID List and LID list.  Initialize to
      // invalid values, in case the user global element list does
      // fill all IDs from minAllGID to maxAllGID (e.g., allows
      // global indices to be all even integers).
      nodeIDs_.resize (dir_numMyEntries, -1);
      LIDs_.resize (dir_numMyEntries, LINVALID);

      // Get list of process IDs that own the directory entries for the
      // Map GIDs.  These will be the targets of the sends that the
      // Distributor will do.
      const int myRank = comm->getRank ();
      const size_t numMyEntries = map->getNodeNumElements ();
      Array<int> sendImageIDs (numMyEntries);
      ArrayView<const GO> myGlobalEntries = map->getNodeElementList ();
      // An ID not present in this lookup indicates that it lies outside
      // of the range [minAllGID,maxAllGID] (from map_).  this means
      // something is wrong with map_, our fault.
      TEUCHOS_TEST_FOR_EXCEPTION(
        directoryMap_->getRemoteIndexList (myGlobalEntries, sendImageIDs) == IDNotPresent,
        std::logic_error, Teuchos::typeName(*this) << " constructor: the "
        "Directory Map could not find out where one or more of my Map's "
        "elements should go.  This probably means there is a bug in Map or "
        "Directory.  Please report this bug to the Tpetra developers.");

      // Initialize the distributor using the list of process IDs to
      // which to send.  We'll use the distributor to send out triples
      // of (GID, process ID, LID).  We're sending the entries to the
      // processes that the Directory Map says should own them, which is
      // why we called directoryMap_->getRemoteIndexList() above.
      Distributor distor (comm);
      const size_t numReceives = distor.createFromSends (sendImageIDs);

      // NOTE (mfh 21 Mar 2012) The following code assumes that
      // sizeof(GO) >= sizeof(int) and sizeof(GO) >= sizeof(LO).
      //
      // Create and fill buffer of (GID, process ID, LID) triples to
      // send out.  We pack the (GID, process ID, LID) triples into a
      // single Array of GO, casting the process ID from int to GO and
      // the LID from LO to GO as we do so.
      const int packetSize = 3; // We're sending triples, so packet size is 3.
      Array<GO> exportEntries (packetSize * numMyEntries); // data to send out
      {
        typename Array<GO>::iterator iter = exportEntries.begin();
        for (size_t i=0; i < numMyEntries; ++i) {
          *iter++ = myGlobalEntries[i];
          *iter++ = as<GO> (myRank);
          *iter++ = as<GO> (i);
        }
      }
      // Buffer of data to receive.  The Distributor figured out for
      // us how many packets we're receiving, when we called its
      // createFromSends() method to set up the distribution plan.
      Array<GO> importElements (packetSize * distor.getTotalReceiveLength ());

      // Distribute the triples of (GID, process ID, LID).
      distor.doPostsAndWaits (exportEntries ().getConst (), packetSize, importElements ());

      // Unpack the redistributed data.
      {
        typename Array<GO>::iterator iter = importElements.begin();
        for (size_t i = 0; i < numReceives; ++i) {
          // Each "packet" (contiguous chunk of importElements) contains
          // a triple: (GID, process ID, LID).
          //
          // Convert incoming GID to Directory LID.
          const GO currGID = *iter++;
          const LO currLID = directoryMap_->getLocalElement (currGID);
          TEUCHOS_TEST_FOR_EXCEPTION(currLID == LINVALID, std::logic_error,
            Teuchos::typeName(*this) << " constructor: Incoming global ID "
            << currGID << " does not have a corresponding local ID in the "
            "Directory Map.  Please report this bug to the Tpetra developers.");
          nodeIDs_[currLID] = *iter++;
          LIDs_[currLID]    = *iter++;
        }
      }
    }

    template<class LO, class GO, class NT>
    std::string
    DistributedNoncontiguousDirectory<LO, GO, NT>::description () const
    {
      std::ostringstream os;
      os << "DistributedNoncontiguousDirectory"
         << "<" << Teuchos::TypeNameTraits<LO>::name ()
         << ", " << Teuchos::TypeNameTraits<GO>::name ()
         << ", " << Teuchos::TypeNameTraits<NT>::name () << ">";
      return os.str ();
    }

    template<class LO, class GO, class NT>
    LookupStatus
    DistributedNoncontiguousDirectory<LO, GO, NT>::
    getEntriesImpl (const Teuchos::ArrayView<const GO> &globalIDs,
                    const Teuchos::ArrayView<int> &nodeIDs,
                    const Teuchos::ArrayView<LO> &localIDs,
                    const bool computeLIDs) const
    {
      using Teuchos::Array;
      using Teuchos::ArrayRCP;
      using Teuchos::ArrayView;
      using Teuchos::as;
      using Teuchos::RCP;

      const LO LINVALID = Teuchos::OrdinalTraits<LO>::invalid();
      LookupStatus res = AllIDsPresent;

      RCP<const map_type> map = this->getMap ();
      RCP<const Teuchos::Comm<int> > comm = map->getComm ();
      const size_t numEntries = globalIDs.size ();

      //
      // Set up directory structure.
      //

      // If we're computing LIDs, we also have to include them in each
      // packet, along with the GID and process ID.
      const int packetSize = computeLIDs ? 3 : 2;

      // For data distribution, we use: Surprise!  A Distributor!
      Distributor distor (comm);

      // Get directory locations for the requested list of entries.
      Array<int> dirImages (numEntries);
      res = directoryMap_->getRemoteIndexList (globalIDs, dirImages ());
      // Check for unfound globalIDs and set corresponding nodeIDs to -1
      size_t numMissing = 0;
      if (res == IDNotPresent) {
        for (size_t i=0; i < numEntries; ++i) {
          if (dirImages[i] == -1) {
            nodeIDs[i] = -1;
            if (computeLIDs) {
              localIDs[i] = LINVALID;
            }
            numMissing++;
          }
        }
      }

      ArrayRCP<GO> sendGIDs;
      ArrayRCP<int> sendImages;
      distor.createFromRecvs (globalIDs, dirImages (), sendGIDs, sendImages);
      const size_t numSends = sendGIDs.size ();

      //    global_size_t >= GO
      //    global_size_t >= size_t >= int
      //    global_size_t >= size_t >= LO
      // Therefore, we can safely store all of these in a global_size_t
      Array<global_size_t> exports (packetSize * numSends);
      {
        LO curLID;
        typename Array<global_size_t>::iterator exportsIter = exports.begin();
        typename ArrayRCP<GO>::const_iterator gidIter = sendGIDs.begin();
        for ( ; gidIter != sendGIDs.end(); ++gidIter) {
          *exportsIter++ = as<global_size_t> (*gidIter);
          curLID = directoryMap_->getLocalElement (*gidIter);
          TEUCHOS_TEST_FOR_EXCEPTION(curLID == LINVALID, std::logic_error,
            Teuchos::typeName (*this) << "::getEntriesImpl(): The Directory "
            "Map's global ID " << *gidIter << " does not have a corresponding "
            "local ID.  Please report this bug to the Tpetra developers.");
          *exportsIter++ = as<global_size_t> (nodeIDs_[curLID]);
          if (computeLIDs) {
            *exportsIter++ = as<global_size_t> (LIDs_[curLID]);
          }
        }
      }

      Array<global_size_t> imports (packetSize * distor.getTotalReceiveLength ());
      distor.doPostsAndWaits (exports ().getConst (), packetSize, imports ());

      typename Array<global_size_t>::iterator ptr = imports.begin();
      const size_t numRecv = numEntries - numMissing;

      Array<GO> sortedIDs (globalIDs);
      ArrayRCP<GO> offset = arcp<GO> (numEntries);
      GO ii=0;
      for (typename ArrayRCP<GO>::iterator oo = offset.begin(); oo != offset.end(); ++oo, ++ii) {
        *oo = ii;
      }
      sort2 (sortedIDs.begin(), sortedIDs.begin() + numEntries, offset.begin());

      typedef typename Array<GO>::iterator IT;
      // we know these conversions are in range, because we loaded this data
      for (size_t i = 0; i < numRecv; ++i) {
        GO curGID = as<GO> (*ptr++);
        std::pair<IT, IT> p1 = std::equal_range (sortedIDs.begin(), sortedIDs.end(), curGID);
        if (p1.first != p1.second) {
          //found it
          size_t j = p1.first - sortedIDs.begin();
          nodeIDs[offset[j]] = as<int>(*ptr++);
          if (computeLIDs) {
            localIDs[offset[j]] = as<LO>(*ptr++);
          }
          if (nodeIDs[offset[j]] == -1) {
            res = IDNotPresent;
          }
        }
      }
      return res;
    }
  } // namespace Details
} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra::Details namespace!
//
#define TPETRA_DIRECTORY_IMPL_INSTANT(LO,GO,NODE)                     \
  template<> class Directory< LO , GO , NODE >;                       \
  template<> class ReplicatedDirectory< LO , GO , NODE >;             \
  template<> class DistributedContiguousDirectory< LO , GO , NODE >;  \
  template<> class DistributedNoncontiguousDirectory< LO , GO , NODE >;

#endif // __Tpetra_DirectoryImpl_def_hpp
