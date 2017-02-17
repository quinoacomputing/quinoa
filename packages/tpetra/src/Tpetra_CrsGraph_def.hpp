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

// FINISH: need to check that fill is active before performing a number of the methods here; adding to the tests currently

#ifndef TPETRA_CRSGRAPH_DEF_HPP
#define TPETRA_CRSGRAPH_DEF_HPP

#include <Kokkos_NodeTrace.hpp>
#include <Teuchos_Assert.hpp>
#include <Teuchos_NullIteratorTraits.hpp>
#include <Teuchos_as.hpp>
#include <algorithm>
#include <string>
#include <utility>
#include <Teuchos_SerialDenseMatrix.hpp>

#ifdef DOXYGEN_USE_ONLY
  #include "Tpetra_CrsGraph_decl.hpp"
#endif

namespace Tpetra {

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  CrsGraph (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap,
            size_t maxNumEntriesPerRow,
            ProfileType pftype,
            const RCP<ParameterList>& params)
  : DistObject<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Node>(rowMap)
  , rowMap_(rowMap)
  , nodeNumEntries_(0)
  , nodeNumAllocated_(OrdinalTraits<size_t>::invalid())
  , pftype_(pftype)
  , numAllocForAllRows_(maxNumEntriesPerRow)
  , indicesAreAllocated_(false)
  , indicesAreLocal_(false)
  , indicesAreGlobal_(false)
  {
    typedef Teuchos::OrdinalTraits<size_t> OTST;
    staticAssertions();
    TEUCHOS_TEST_FOR_EXCEPTION(maxNumEntriesPerRow == OTST::invalid(),
      std::invalid_argument, "The allocation hint must be a valid size_t value, "
      "which in this case means it must not be Teuchos::OrdinalTraits<size_t>::"
      "invalid().");
    resumeFill(params);
    checkInternalState();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  CrsGraph (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap,
            const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &colMap,
            size_t maxNumEntriesPerRow,
            ProfileType pftype,
            const RCP<ParameterList>& params)
  : DistObject<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Node>(rowMap)
  , rowMap_(rowMap)
  , colMap_(colMap)
  , nodeNumEntries_(0)
  , nodeNumAllocated_(OrdinalTraits<size_t>::invalid())
  , pftype_(pftype)
  , numAllocForAllRows_(maxNumEntriesPerRow)
  , indicesAreAllocated_(false)
  , indicesAreLocal_(false)
  , indicesAreGlobal_(false)
  {
    typedef Teuchos::OrdinalTraits<size_t> OTST;
    staticAssertions();
    TEUCHOS_TEST_FOR_EXCEPTION(maxNumEntriesPerRow == OTST::invalid(),
      std::invalid_argument, "The allocation hint must be a valid size_t value, "
      "which in this case means it must not be Teuchos::OrdinalTraits<size_t>::"
      "invalid().");
    resumeFill(params);
    checkInternalState();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  CrsGraph (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap,
            const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc,
            ProfileType pftype,
            const RCP<ParameterList>& params)
  : DistObject<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Node>(rowMap)
  , rowMap_(rowMap)
  , nodeNumEntries_(0)
  , nodeNumAllocated_(OrdinalTraits<size_t>::invalid())
  , pftype_(pftype)
  , numAllocPerRow_(NumEntriesPerRowToAlloc)
  , numAllocForAllRows_(0)
  , indicesAreAllocated_(false)
  , indicesAreLocal_(false)
  , indicesAreGlobal_(false)
  {
    typedef Teuchos::OrdinalTraits<size_t> OTST;
    std::string tfecfFuncName("CrsGraph(rowMap,NumEntriesPerRowToAlloc)");
    staticAssertions();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC((size_t)NumEntriesPerRowToAlloc.size() != getNodeNumRows(), std::invalid_argument,
        ": NumEntriesPerRowToAlloc must have as many entries as specified by rowMap for this node.");
    //size_t numMin = OTST::max();
    //size_t numMax = OTST::zero();
    for (size_t r=0; r < getNodeNumRows(); ++r) {
      const size_t curRowCount = NumEntriesPerRowToAlloc[r];
      TEUCHOS_TEST_FOR_EXCEPTION(curRowCount == OTST::invalid(),
        std::invalid_argument, "NumEntriesPerRowToAlloc[" << r << "] specifies "
        "an invalid number of entries (Teuchos::OrdinalTraits<size_t>::"
        "invalid()).");
      //numMin = std::min<size_t> (numMin, curRowCount);
      //numMax = std::max<size_t> (numMax, curRowCount);
    }
    resumeFill(params);
    checkInternalState();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  CrsGraph (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap,
            const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &colMap,
            const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc,
            ProfileType pftype,
            const RCP<ParameterList>& params)
  : DistObject<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Node>(rowMap)
  , rowMap_(rowMap)
  , colMap_(colMap)
  , nodeNumEntries_(0)
  , nodeNumAllocated_(OrdinalTraits<size_t>::invalid())
  , pftype_(pftype)
  , numAllocPerRow_(NumEntriesPerRowToAlloc)
  , numAllocForAllRows_(0)
  , indicesAreAllocated_(false)
  , indicesAreLocal_(false)
  , indicesAreGlobal_(false)
  {
    typedef Teuchos::OrdinalTraits<size_t> OTST;
    std::string tfecfFuncName("CrsGraph(rowMap,colMap,NumEntriesPerRowToAlloc)");
    staticAssertions();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC((size_t)NumEntriesPerRowToAlloc.size() != getNodeNumRows(), std::invalid_argument,
        ": NumEntriesPerRowToAlloc must have as many entries as specified by rowMap for this node.");
    //size_t numMin = OTST::max();
    //size_t numMax = OTST::zero();
    for (size_t r=0; r < getNodeNumRows(); ++r) {
      const size_t curRowCount = NumEntriesPerRowToAlloc[r];
      TEUCHOS_TEST_FOR_EXCEPTION(curRowCount == OTST::invalid(),
        std::invalid_argument, "NumEntriesPerRowToAlloc[" << r << "] specifies "
        "an invalid number of entries (Teuchos::OrdinalTraits<size_t>::"
        "invalid()).");
      //numMin = std::min<size_t> (numMin, curRowCount);
      //numMax = std::max<size_t> (numMax, curRowCount);
    }
    resumeFill(params);
    checkInternalState();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::~CrsGraph()
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const ParameterList>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  getValidParameters () const
  {
    RCP<ParameterList> params = parameterList ("Tpetra::CrsGraph");

    // Make a sublist for the Import.
    RCP<ParameterList> importSublist = parameterList ("Import");

    // FIXME (mfh 02 Apr 2012) We should really have the Import and
    // Export objects fill in these lists.  However, we don't want to
    // create an Import or Export unless we need them.  For now, we
    // know that the Import and Export just pass the list directly to
    // their Distributor, so we can create a Distributor here
    // (Distributor's constructor is a lightweight operation) and have
    // it fill in the list.

    // Fill in Distributor default parameters by creating a
    // Distributor and asking it to do the work.
    Distributor distributor (rowMap_->getComm(), importSublist);
    params->set ("Import", *importSublist, "How the Import performs communication.");

    // Make a sublist for the Export.  For now, it's a clone of the
    // Import sublist.  It's not a shallow copy, though, since we
    // might like the Import to do communication differently than the
    // Export.
    params->set ("Export", *importSublist, "How the Export performs communication.");

    return params;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  setParameterList (const RCP<ParameterList>& params)
  {
    RCP<const ParameterList> validParams = getValidParameters ();
    params->validateParametersAndSetDefaults (*validParams);
    this->setMyParamList (params);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  global_size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getGlobalNumRows() const
  {
    return rowMap_->getGlobalNumElements();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  global_size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getGlobalNumCols() const
  {
    std::string tfecfFuncName("getGlobalNumCols()");
    TEUCHOS_TEST_FOR_EXCEPTION(isFillComplete() == false, std::runtime_error, ": requires domain map (which requires fillComplete()).");
    return getDomainMap()->getGlobalNumElements();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getNodeNumRows() const
  {
    return rowMap_->getNodeNumElements();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getNodeNumCols() const
  {
    std::string tfecfFuncName("getNodeNumCols()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(hasColMap() != true, std::runtime_error, ": requires column map.");
    return colMap_->getNodeNumElements();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getNodeNumDiags() const
  {
    return nodeNumDiags_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  global_size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getGlobalNumDiags() const
  {
    return globalNumDiags_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Node>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getNode() const
  {
    return rowMap_->getNode();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getRowMap() const
  {
    return rowMap_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getColMap() const
  {
    return colMap_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getDomainMap() const
  {
    return domainMap_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getRangeMap() const
  {
    return rangeMap_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> >
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getImporter() const
  {
    return importer_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> >
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getExporter() const
  {
    return exporter_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::hasColMap() const
  {
    return (colMap_ != null);
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::isStorageOptimized() const
  {
    bool isOpt = indicesAreAllocated_ && (numRowEntries_ == null) && (getNodeNumRows() > 0);
#ifdef HAVE_TPETRA_DEBUG
    std::string tfecfFuncName("isStorageOptimized()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( (isOpt == true) && (getProfileType() == DynamicProfile), std::logic_error,
        ": Violated stated post-conditions. Please contact Tpetra team.");
#endif
    return isOpt;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  ProfileType CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getProfileType() const
  {
    return pftype_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  global_size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getGlobalNumEntries() const
  {
    return globalNumEntries_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getNodeNumEntries() const
  {
    return nodeNumEntries_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  global_size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getGlobalMaxNumRowEntries() const
  {
    return globalMaxNumRowEntries_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getNodeMaxNumRowEntries() const
  {
    return nodeMaxNumRowEntries_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::isFillComplete() const
  {
    return fillComplete_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::isFillActive() const
  {
    return !fillComplete_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::isUpperTriangular() const
  {
    return upperTriangular_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::isLowerTriangular() const
  {
    return lowerTriangular_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::isLocallyIndexed() const
  {
    return indicesAreLocal_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::isGloballyIndexed() const
  {
    return indicesAreGlobal_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getNodeAllocationSize() const
  {
    return nodeNumAllocated_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  const RCP<const Comm<int> > &
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getComm() const
  {
    return rowMap_->getComm();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  GlobalOrdinal CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getIndexBase() const
  {
    return rowMap_->getIndexBase();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::indicesAreAllocated() const
  {
    return indicesAreAllocated_;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                    Internal utility methods                             //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::isSorted() const
  {
    return indicesAreSorted_;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::isMerged() const
  {
    return noRedundancies_;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setSorted(bool sorted)
  {
    indicesAreSorted_ = sorted;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setMerged(bool merged)
  {
    noRedundancies_ = merged;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::allocateIndices(ELocalGlobal lg)
  {
    // this is a protected function, only callable by us. if it was called incorrectly, it is our fault.
    std::string tfecfFuncName("allocateIndices()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( isLocallyIndexed()  && lg==GlobalIndices, std::logic_error, ": Internal logic error. Please contact Tpetra team.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( isGloballyIndexed() && lg==LocalIndices,  std::logic_error, ": Internal logic error. Please contact Tpetra team.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( indicesAreAllocated() == true,            std::logic_error, ": Internal logic error. Please contact Tpetra team.");
    const size_t numRows = getNodeNumRows();
    indicesAreLocal_  = (lg == LocalIndices);
    indicesAreGlobal_ = (lg == GlobalIndices);
    nodeNumAllocated_ = 0;
    if (numAllocPerRow_ == null && getNodeNumRows() > 0) {
      // this wastes memory, temporarily, but it simplifies the code and interfaces to follow
      // TODO: it is a candidate for change
      ArrayRCP<size_t> tmpnumallocperrow = arcp<size_t>(numRows);
      std::fill(tmpnumallocperrow.begin(), tmpnumallocperrow.end(), numAllocForAllRows_);
      numAllocPerRow_ = tmpnumallocperrow;
    }
    if (getProfileType() == StaticProfile) {
      //
      //  STATIC ALLOCATION PROFILE
      //
      // have the graph do this for us, with thread appropriate allocation if possible
      // not important for globals, because we never multiply out of them, but we might as well use the code that we have
      rowPtrs_ = LocalMatOps::allocRowPtrs( getRowMap()->getNode(), numAllocPerRow_() );
      if (lg == LocalIndices) lclInds1D_ = LocalMatOps::template allocStorage< LocalOrdinal>( getRowMap()->getNode(), rowPtrs_() );
      else                    gblInds1D_ = LocalMatOps::template allocStorage<GlobalOrdinal>( getRowMap()->getNode(), rowPtrs_() );
      nodeNumAllocated_ = rowPtrs_[numRows];
    }
    else {
      //
      //  DYNAMIC ALLOCATION PROFILE
      //
      ArrayRCP<const size_t>::iterator numalloc = numAllocPerRow_.begin();
      size_t howmany = numAllocForAllRows_;
      if (lg == LocalIndices) {
        lclInds2D_ = arcp< ArrayRCP<LocalOrdinal> >(numRows);
        nodeNumAllocated_ = 0;
        for (size_t i=0; i < numRows; ++i) {
          if (numAllocPerRow_ != null) howmany = *numalloc++;
          nodeNumAllocated_ += howmany;
          if (howmany > 0) lclInds2D_[i] = arcp<LocalOrdinal>(howmany);
        }
      }
      else { // allocate global indices
        gblInds2D_ = arcp< ArrayRCP<GlobalOrdinal> >(numRows);
        nodeNumAllocated_ = 0;
        for (size_t i=0; i < numRows; ++i) {
          if (numAllocPerRow_ != null) howmany = *numalloc++;
          nodeNumAllocated_ += howmany;
          if (howmany > 0) gblInds2D_[i] = arcp<GlobalOrdinal>(howmany);
        }
      }
    }
    if (numRows > 0) {
      numRowEntries_ = arcp<size_t>(numRows);
      std::fill( numRowEntries_.begin(), numRowEntries_.end(), (size_t)0 );
    }
    // done with these
    numAllocForAllRows_ = 0;
    numAllocPerRow_     = null;
    indicesAreAllocated_ = true;
    checkInternalState();
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  template <class T>
  ArrayRCP<T> CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::allocateValues1D() const
  {
    std::string tfecfFuncName("allocateValues()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( 
        indicesAreAllocated() == false, 
        std::runtime_error, ": graph indices must be allocated before values."
    );
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( 
        getProfileType() != StaticProfile, 
        std::runtime_error, ": graph indices must be allocated in a static profile."
    );
    ArrayRCP<T> values1D;
    values1D = LocalMatOps::template allocStorage<T>( getRowMap()->getNode(), rowPtrs_() );
    return values1D;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  template <class T>
  ArrayRCP<ArrayRCP<T> > CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::allocateValues2D() const
  {
    std::string tfecfFuncName("allocateValues2D()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( 
        indicesAreAllocated() == false, 
        std::runtime_error, ": graph indices must be allocated before values."
    );
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( 
        getProfileType() != DynamicProfile, 
        std::runtime_error, ": graph indices must be allocated in a dynamic profile."
    );
    ArrayRCP<ArrayRCP<T> > values2D;
    values2D = arcp<ArrayRCP<T> >(getNodeNumRows());
    if (lclInds2D_ != null) {
      for (size_t r=0; r < (size_t)lclInds2D_.size(); ++r) {
        if (lclInds2D_[r] != null) {
          values2D[r] = arcp<T>(lclInds2D_[r].size());
        }
      }
    }
    else if (gblInds2D_ != null) {
      for (size_t r=0; r < (size_t)gblInds2D_.size(); ++r) {
        if (gblInds2D_[r] != null) {
          values2D[r] = arcp<T>(gblInds2D_[r].size());
        }
      }
    }
    return values2D;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  template <ELocalGlobal lg>
  RowInfo CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::updateAlloc(RowInfo rowinfo, size_t newAllocSize)
  {
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPT( rowMap_->isNodeLocalElement(rowinfo.localRow) == false );
    TEUCHOS_TEST_FOR_EXCEPT( newAllocSize < rowinfo.allocSize );
    TEUCHOS_TEST_FOR_EXCEPT( (lg == LocalIndices && isLocallyIndexed() == false) || (lg == GlobalIndices && isGloballyIndexed() == false) );
    TEUCHOS_TEST_FOR_EXCEPT( newAllocSize == 0 );
    TEUCHOS_TEST_FOR_EXCEPT( indicesAreAllocated() == false );
#endif
    // allocate a larger space for row "lrow"
    // copy any existing data from previous allocation to new allocation
    // update sizes
    if (lg == LocalIndices) {
      ArrayRCP<LocalOrdinal> old_alloc = lclInds2D_[rowinfo.localRow];
      lclInds2D_[rowinfo.localRow] = arcp<LocalOrdinal>(newAllocSize);
      std::copy(old_alloc.begin(), old_alloc.begin() + rowinfo.numEntries, lclInds2D_[rowinfo.localRow].begin());
    }
    else /* if lg == GlobalIndices */ {
      ArrayRCP<GlobalOrdinal> old_alloc = gblInds2D_[rowinfo.localRow];
      gblInds2D_[rowinfo.localRow] = arcp<GlobalOrdinal>(newAllocSize);
      std::copy(old_alloc.begin(), old_alloc.begin() + rowinfo.numEntries, gblInds2D_[rowinfo.localRow].begin());
    }
    //
    nodeNumAllocated_ += (newAllocSize - rowinfo.allocSize);
    rowinfo.allocSize = newAllocSize;
    return rowinfo;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  template <ELocalGlobal lg, class T>
  RowInfo CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::updateAllocAndValues(RowInfo rowinfo, size_t newAllocSize, ArrayRCP<T> &rowVals)
  {
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPT( ! rowMap_->isNodeLocalElement(rowinfo.localRow) );
    TEUCHOS_TEST_FOR_EXCEPT( newAllocSize < rowinfo.allocSize );
    TEUCHOS_TEST_FOR_EXCEPT( (lg == LocalIndices && ! isLocallyIndexed()) || (lg == GlobalIndices && ! isGloballyIndexed()) );
    TEUCHOS_TEST_FOR_EXCEPT( newAllocSize == 0 );
    TEUCHOS_TEST_FOR_EXCEPT( ! indicesAreAllocated() );
#endif
    // allocate a larger space for row "lrow"
    // copy any existing data from previous allocation to new allocation
    // update sizes
    if (lg == LocalIndices) {
      ArrayRCP<LocalOrdinal> old_alloc = lclInds2D_[rowinfo.localRow];
      lclInds2D_[rowinfo.localRow] = arcp<LocalOrdinal>(newAllocSize);
      std::copy(old_alloc.begin(), old_alloc.begin() + rowinfo.numEntries, lclInds2D_[rowinfo.localRow].begin());
    }
    else /* if lg == GlobalIndices */ {
      ArrayRCP<GlobalOrdinal> old_alloc = gblInds2D_[rowinfo.localRow];
      gblInds2D_[rowinfo.localRow] = arcp<GlobalOrdinal>(newAllocSize);
      std::copy(old_alloc.begin(), old_alloc.begin() + rowinfo.numEntries, gblInds2D_[rowinfo.localRow].begin());
    }
    ArrayRCP<const T> oldVals = rowVals;
    rowVals = arcp<T>(newAllocSize);
    std::copy(oldVals.begin(), oldVals.begin() + rowinfo.numEntries, rowVals.begin());
    //
    nodeNumAllocated_ += (newAllocSize - rowinfo.allocSize);
    rowinfo.allocSize = newAllocSize;
    return rowinfo;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  ArrayView<const LocalOrdinal> CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getLocalView(RowInfo rowinfo) const
  {
    ArrayView<const LocalOrdinal> view;
    if (rowinfo.allocSize > 0) {
      if (lclInds1D_ != null) {
        view = lclInds1D_(rowinfo.offset1D,rowinfo.allocSize);
      }
      else if (lclInds2D_[rowinfo.localRow] != null) {
        view = lclInds2D_[rowinfo.localRow]();
      }
    }
    return view;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  ArrayView<LocalOrdinal> CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getLocalViewNonConst(RowInfo rowinfo)
  {
    ArrayView<LocalOrdinal> view;
    if (rowinfo.allocSize > 0) {
      if (lclInds1D_ != null) {
        view = lclInds1D_(rowinfo.offset1D,rowinfo.allocSize);
      }
      else if (lclInds2D_[rowinfo.localRow] != null) {
        view = lclInds2D_[rowinfo.localRow]();
      }
    }
    return view;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  ArrayView<const GlobalOrdinal> CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getGlobalView(RowInfo rowinfo) const
  {
    ArrayView<const GlobalOrdinal> view;
    if (rowinfo.allocSize > 0) {
      if (gblInds1D_ != null) {
        view = gblInds1D_(rowinfo.offset1D,rowinfo.allocSize);
      }
      else if (gblInds2D_[rowinfo.localRow] != null) {
        view = gblInds2D_[rowinfo.localRow]();
      }
    }
    return view;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  ArrayView<GlobalOrdinal> CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getGlobalViewNonConst(RowInfo rowinfo)
  {
    ArrayView<GlobalOrdinal> view;
    if (rowinfo.allocSize > 0) {
      if (gblInds1D_ != null) {
        view = gblInds1D_(rowinfo.offset1D,rowinfo.allocSize);
      }
      else if (gblInds2D_[rowinfo.localRow] != null) {
        view = gblInds2D_[rowinfo.localRow]();
      }
    }
    return view;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RowInfo CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getRowInfo(size_t myRow) const
  {
#ifdef HAVE_TPETRA_DEBUG
    std::string tfecfFuncName("getRowInfo()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(rowMap_->isNodeLocalElement(myRow) == false, std::logic_error, ": Internal logic error. Please contact Tpetra team.");
#endif
    const size_t STINV = OrdinalTraits<size_t>::invalid();
    RowInfo ret;
    ret.localRow = myRow;
    if (nodeNumAllocated_ != 0 && nodeNumAllocated_ != STINV) {
      // graph data structures have the info that we need
      //
      // if static graph, offsets tell us the allocation size
      if (getProfileType() == StaticProfile) {
        ret.offset1D   = rowPtrs_[myRow];
        ret.allocSize  = rowPtrs_[myRow+1] - rowPtrs_[myRow];
        if (numRowEntries_ == null) ret.numEntries = ret.allocSize;
        else                        ret.numEntries = numRowEntries_[myRow];
      }
      else {
        ret.offset1D = STINV;
        if (isLocallyIndexed()) {
          ret.allocSize = lclInds2D_[myRow].size();
        }
        else {
          ret.allocSize = gblInds2D_[myRow].size();
        }
        ret.numEntries = numRowEntries_[myRow];
      }
    }
    else if (nodeNumAllocated_ == 0) {
      // have performed allocation, but the graph has no allocation or entries
      ret.allocSize = 0;
      ret.numEntries = 0;
      ret.offset1D = STINV;
    }
    else if (indicesAreAllocated() == false) {
      // haven't performed allocation yet; probably won't hit this code
      if (numAllocPerRow_ == null) {
        ret.allocSize = numAllocForAllRows_;
      }
      else {
        ret.allocSize = numAllocPerRow_[myRow];
      }
      ret.numEntries = 0;
      ret.offset1D = STINV;
    }
    else {
      // don't know how we ended up here...
      TEUCHOS_TEST_FOR_EXCEPT(true);
    }
    return ret;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::staticAssertions() const
  {
    // Assumption: sizeof(GlobalOrdinal) >= sizeof(LocalOrdinal)
    //    This is so that we can store LocalOrdinals in the memory formerly occupied by GlobalOrdinals
    // Assumption: max(GlobalOrdinal) >= max(LocalOrdinal)  and  max(size_t) >= max(LocalOrdinal)
    //    This is so that we can represent any LocalOrdinal as a size_t, and any LocalOrdinal as a GlobalOrdinal
    Teuchos::CompileTimeAssert<sizeof(GlobalOrdinal) < sizeof(LocalOrdinal)> cta_size1; (void)cta_size1;
    Teuchos::CompileTimeAssert<sizeof(global_size_t) < sizeof(size_t)      > cta_size2; (void)cta_size2;
    // can't call max() with CompileTimeAssert, because it isn't a constant expression; will need to make this a runtime check
    std::string msg = typeName(*this) + ": Object cannot be allocated with stated template arguments: size assumptions are not valid.";
    TEUCHOS_TEST_FOR_EXCEPTION( (size_t)OrdinalTraits<LocalOrdinal>::max() > OrdinalTraits<size_t>::max(),          std::runtime_error, msg);
    TEUCHOS_TEST_FOR_EXCEPTION( (global_size_t)OrdinalTraits<LocalOrdinal>::max() > (global_size_t)OrdinalTraits<GlobalOrdinal>::max(),           std::runtime_error, msg);
    TEUCHOS_TEST_FOR_EXCEPTION( (size_t)OrdinalTraits<GlobalOrdinal>::max() > OrdinalTraits<global_size_t>::max(),  std::runtime_error, msg);
    TEUCHOS_TEST_FOR_EXCEPTION( OrdinalTraits<size_t>::max() > OrdinalTraits<global_size_t>::max(),                 std::runtime_error, msg);
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  template <ELocalGlobal lg>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::filterIndices(const SLocalGlobalNCViews &inds) const
  {
    const Map<LocalOrdinal,GlobalOrdinal,Node> &cmap = *colMap_;
    Teuchos::CompileTimeAssert<lg != GlobalIndices && lg != LocalIndices> cta_lg;
    (void)cta_lg;
    size_t numFiltered = 0;
#ifdef HAVE_TPETRA_DEBUG
    size_t numFiltered_debug = 0;
#endif
    if (lg == GlobalIndices) {
      ArrayView<GlobalOrdinal> ginds = inds.ginds;
      typename ArrayView<GlobalOrdinal>::iterator fend = ginds.begin(),
                                                  cptr = ginds.begin();
      while (cptr != ginds.end()) {
        if (cmap.isNodeGlobalElement(*cptr)) {
          *fend++ = *cptr;
#ifdef HAVE_TPETRA_DEBUG
          ++numFiltered_debug;
#endif
        }
        ++cptr;
      }
      numFiltered = fend - ginds.begin();
    }
    else if (lg == LocalIndices) {
      ArrayView<LocalOrdinal> linds = inds.linds;
      typename ArrayView<LocalOrdinal>::iterator fend = linds.begin(),
                                                 cptr = linds.begin();
      while (cptr != linds.end()) {
        if (cmap.isNodeLocalElement(*cptr)) {
          *fend++ = *cptr;
#ifdef HAVE_TPETRA_DEBUG
          ++numFiltered_debug;
#endif
        }
        ++cptr;
      }
      numFiltered = fend - linds.begin();
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPT( numFiltered != numFiltered_debug );
#endif
    return numFiltered;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  template <ELocalGlobal lg, class T>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::filterIndicesAndValues(const SLocalGlobalNCViews &inds, const ArrayView<T> &vals) const
  {
    const Map<LocalOrdinal,GlobalOrdinal,Node> &cmap = *colMap_;
    Teuchos::CompileTimeAssert<lg != GlobalIndices && lg != LocalIndices> cta_lg;
    (void)cta_lg;
    size_t numFiltered = 0;
    typename ArrayView<T>::iterator fvalsend = vals.begin(),
                                    valscptr = vals.begin();
#ifdef HAVE_TPETRA_DEBUG
    size_t numFiltered_debug = 0;
#endif
    if (lg == GlobalIndices) {
      ArrayView<GlobalOrdinal> ginds = inds.ginds;
      typename ArrayView<GlobalOrdinal>::iterator fend = ginds.begin(),
                                                  cptr = ginds.begin();
      while (cptr != ginds.end()) {
        if (cmap.isNodeGlobalElement(*cptr)) {
          *fend++ = *cptr;
          *fvalsend++ = *valscptr;
#ifdef HAVE_TPETRA_DEBUG
          ++numFiltered_debug;
#endif
        }
        ++cptr;
        ++valscptr;
      }
      numFiltered = fend - ginds.begin();
    }
    else if (lg == LocalIndices) {
      ArrayView<LocalOrdinal> linds = inds.linds;
      typename ArrayView<LocalOrdinal>::iterator fend = linds.begin(),
                                                 cptr = linds.begin();
      while (cptr != linds.end()) {
        if (cmap.isNodeLocalElement(*cptr)) {
          *fend++ = *cptr;
          *fvalsend++ = *valscptr;
#ifdef HAVE_TPETRA_DEBUG
          ++numFiltered_debug;
#endif
        }
        ++cptr;
        ++valscptr;
      }
      numFiltered = fend - linds.begin();
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPT( numFiltered != numFiltered_debug );
    TEUCHOS_TEST_FOR_EXCEPT( valscptr != vals.end() );
    TEUCHOS_TEST_FOR_EXCEPT( numFiltered != (size_t)(fvalsend - vals.begin()) );
#endif
    return numFiltered;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  template <ELocalGlobal lg, ELocalGlobal I>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::insertIndices(RowInfo rowinfo, const SLocalGlobalViews &newInds)
  {
    Teuchos::CompileTimeAssert<lg != GlobalIndices && lg != LocalIndices> cta_lg; (void)cta_lg;
    size_t numNewInds = 0;
    if (lg == GlobalIndices) {
      ArrayView<const GlobalOrdinal> new_ginds = newInds.ginds;
      numNewInds = new_ginds.size();
      if (I == GlobalIndices) {
        ArrayView<GlobalOrdinal> gind_view = getGlobalViewNonConst(rowinfo);
        std::copy(new_ginds.begin(), new_ginds.end(), gind_view.begin()+rowinfo.numEntries);
      }
      else if (I == LocalIndices) {
        ArrayView<LocalOrdinal> lind_view = getLocalViewNonConst(rowinfo);
        typename ArrayView<const GlobalOrdinal>::iterator         in = new_ginds.begin();
        const typename ArrayView<const GlobalOrdinal>::iterator stop = new_ginds.end();
        typename ArrayView<LocalOrdinal>::iterator out = lind_view.begin()+rowinfo.numEntries;
        while (in != stop) {
          *out++ = colMap_->getLocalElement (*in++);
        }
      }
    }
    else if (lg == LocalIndices) {
      ArrayView<const LocalOrdinal> new_linds = newInds.linds;
      numNewInds = new_linds.size();
      if (I == LocalIndices) {
        ArrayView<LocalOrdinal> lind_view = getLocalViewNonConst(rowinfo);
        std::copy(new_linds.begin(), new_linds.end(), lind_view.begin()+rowinfo.numEntries);
      }
      else if (I == GlobalIndices) {
        // not needed yet
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Tpetra::CrsGraph::"
          "insertIndices: the case where the input indices are local and the "
          "indices to write are global (lg=LocalIndices, I=GlobalIndices) has "
          "not yet been implemented.");
      }
    }
    numRowEntries_[rowinfo.localRow] += numNewInds;
    nodeNumEntries_ += numNewInds;
    setSorted(false);
    noRedundancies_ = false;
    return numNewInds;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  template <ELocalGlobal lg, ELocalGlobal I, class IterO, class IterN>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::insertIndicesAndValues(RowInfo rowinfo, const SLocalGlobalViews &newInds, IterO rowVals, IterN newVals)
  {
    size_t numNewInds = insertIndices<lg,I>(rowinfo,newInds);
    std::copy( newVals, newVals + numNewInds, rowVals + rowinfo.numEntries );
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  template <ELocalGlobal lg, class IterO, class IterN, class BinaryFunction>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::transformValues(RowInfo rowinfo, const SLocalGlobalViews &inds, IterO rowVals, IterN newVals, BinaryFunction f) const
  {
    Teuchos::CompileTimeAssert<lg != GlobalIndices && lg != LocalIndices> cta_lg; (void)cta_lg;
    const size_t STINV = OrdinalTraits<size_t>::invalid();
    if (lg == GlobalIndices) {
      ArrayView<const GlobalOrdinal> search_ginds = inds.ginds;
      for (size_t j=0; j < (size_t)search_ginds.size(); ++j) {
        const size_t k = findGlobalIndex(rowinfo, search_ginds[j]);
        if (k != STINV) {
          rowVals[k] = f( rowVals[k], newVals[j] );
        }
      }
    }
    else if (lg == LocalIndices) {
      ArrayView<const LocalOrdinal> search_linds = inds.linds;
      for (size_t j=0; j < (size_t)search_linds.size(); ++j) {
        const size_t k = findLocalIndex(rowinfo, search_linds[j]);
        if (k != STINV) {
          rowVals[k] = f( rowVals[k], newVals[j] );
        }
      }
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::sortRowIndices(RowInfo rowinfo)
  {
    if (rowinfo.numEntries > 0) {
      ArrayView<LocalOrdinal> inds_view = getLocalViewNonConst(rowinfo);
      std::sort(inds_view.begin(), inds_view.begin() + rowinfo.numEntries);
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  // in the future, perhaps this could use std::sort with a boost::zip_iterator
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  template <class Scalar>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::sortRowIndicesAndValues(RowInfo rowinfo, ArrayView<Scalar> values)
  {
    if (rowinfo.numEntries > 0) {
      ArrayView<LocalOrdinal> inds_view = getLocalViewNonConst(rowinfo);
      sort2(inds_view.begin(), inds_view.begin()+rowinfo.numEntries, values.begin());
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::mergeRowIndices(RowInfo rowinfo)
  {
    ArrayView<LocalOrdinal> inds_view = getLocalViewNonConst(rowinfo);
    typename ArrayView<LocalOrdinal>::iterator beg, end, newend;
    beg = inds_view.begin();
    end = inds_view.begin() + rowinfo.numEntries;
    newend = std::unique(beg,end);
    const size_t mergedEntries = newend - beg;
#ifdef HAVE_TPETRA_DEBUG
    // merge should not have eliminated any entries; if so, the assignment below will destory the packed structure
    TEUCHOS_TEST_FOR_EXCEPT( isStorageOptimized() && mergedEntries != rowinfo.numEntries );
#endif
    numRowEntries_[rowinfo.localRow] = mergedEntries;
    nodeNumEntries_ -= (rowinfo.numEntries - mergedEntries);
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  // in the future, this could use std::unique with a boost::zip_iterator
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  template <class Iter, class BinaryFunction>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::mergeRowIndicesAndValues(RowInfo rowinfo, Iter rowValueIter, BinaryFunction f)
  {
    ArrayView<LocalOrdinal> inds_view = getLocalViewNonConst(rowinfo);
    typename ArrayView<LocalOrdinal>::iterator beg, end, newend;
    beg = inds_view.begin();
    end = inds_view.begin() + rowinfo.numEntries;
    newend = beg;
    if (beg != end) {
      typename ArrayView<LocalOrdinal>::iterator cur = beg + 1;
      Iter vcur = rowValueIter + 1,
           vend = rowValueIter;
      cur = beg+1;
      while (cur != end) {
        if (*cur != *newend) {
          // new entry; save it
          ++newend;
          ++vend;
          (*newend) = (*cur);
          (*vend) = (*vcur);
        }
        else {
          // old entry; merge it
          (*vend) = f(*vend,*vcur);
        }
        ++cur;
        ++vcur;
      }
      ++newend; // point one past the last entry, per typical [beg,end) semantics
    }
    const size_t mergedEntries = newend - beg;
#ifdef HAVE_TPETRA_DEBUG
    // merge should not have eliminated any entries; if so, the assignment below will destory the packed structure
    TEUCHOS_TEST_FOR_EXCEPT( isStorageOptimized() && mergedEntries != rowinfo.numEntries );
#endif
    numRowEntries_[rowinfo.localRow] = mergedEntries;
    nodeNumEntries_ -= (rowinfo.numEntries - mergedEntries);
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setDomainRangeMaps(
                                const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap,
                                const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap)
  {
    // simple pointer comparison for equality
    if (domainMap_ != domainMap) {
      domainMap_ = domainMap;
      importer_ = null;
    }
    if (rangeMap_ != rangeMap) {
      rangeMap_  = rangeMap;
      exporter_ = null;
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::findLocalIndex(RowInfo rowinfo, LocalOrdinal ind) const
  {
    typedef typename ArrayView<const LocalOrdinal>::iterator IT;
    bool found = true;
    // get a view of the row, if it wasn't passed by the caller
    ArrayView<const LocalOrdinal> rowinds = getLocalView(rowinfo);
    IT rptr, locptr = Teuchos::NullIteratorTraits<IT>::getNull();
    rptr = rowinds.begin();
    if (isSorted()) {
      // binary search
      std::pair<IT,IT> p = std::equal_range(rptr,rptr+rowinfo.numEntries,ind);
      if (p.first == p.second) found = false;
      else locptr = p.first;
    }
    else {
      // direct search
      locptr = std::find(rptr,rptr+rowinfo.numEntries,ind);
      if (locptr == rptr+rowinfo.numEntries) found = false;
    }
    size_t ret = OrdinalTraits<size_t>::invalid();
    if (found) {
      ret = (locptr - rptr);
    }
    return ret;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::findGlobalIndex(RowInfo rowinfo, GlobalOrdinal ind) const
  {
    typedef typename ArrayView<const GlobalOrdinal>::iterator IT;
    bool found = true;
    // get a view of the row, if it wasn't passed by the caller
    ArrayView<const GlobalOrdinal> rowinds = getGlobalView(rowinfo);
    IT rptr, locptr = Teuchos::NullIteratorTraits<IT>::getNull();
    rptr = rowinds.begin();
    if (isSorted()) {
      // binary search
      std::pair<IT,IT> p = std::equal_range(rptr,rptr+rowinfo.numEntries,ind);
      if (p.first == p.second) found = false;
      else locptr = p.first;
    }
    else {
      // direct search
      locptr = std::find(rptr,rptr+rowinfo.numEntries,ind);
      if (locptr == rptr+rowinfo.numEntries) found = false;
    }
    size_t ret = OrdinalTraits<size_t>::invalid();
    if (found) {
      ret = (locptr - rptr);
    }
    return ret;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::clearGlobalConstants()
  {
    globalNumEntries_       = OrdinalTraits<global_size_t>::invalid();
    globalNumDiags_         = OrdinalTraits<global_size_t>::invalid();
    globalMaxNumRowEntries_ = OrdinalTraits<global_size_t>::invalid();
    nodeNumDiags_           = OrdinalTraits<       size_t>::invalid();
    nodeMaxNumRowEntries_   = OrdinalTraits<       size_t>::invalid();
    haveGlobalConstants_    = false;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::checkInternalState() const
  {
#ifdef HAVE_TPETRA_DEBUG
    const global_size_t GSTI = OrdinalTraits<global_size_t>::invalid();
    const size_t         STI = OrdinalTraits<size_t>::invalid();
    std::string err = typeName(*this) + "::checkInternalState(): Likely internal logic error. Please contact Tpetra team.";
    // check the internal state of this data structure
    // this is called by numerous state-changing methods, in a debug build, to ensure that the object
    // always remains in a valid state
    // the graph should have been allocated with a row map
    TEUCHOS_TEST_FOR_EXCEPTION( rowMap_ == null,                     std::logic_error, err );
    // am either complete or active
    TEUCHOS_TEST_FOR_EXCEPTION( isFillActive() == isFillComplete(),  std::logic_error, err );
    // if active, i have no local graph
    TEUCHOS_TEST_FOR_EXCEPTION( isFillActive() && lclGraph_ != null, std::logic_error, err );
    // if the graph has been fill completed, then all maps should be present
    TEUCHOS_TEST_FOR_EXCEPTION( isFillComplete() == true && (colMap_ == null || rangeMap_ == null || domainMap_ == null), std::logic_error, err );
    // if storage has been optimized, then indices should have been allocated (even if trivially so)
    TEUCHOS_TEST_FOR_EXCEPTION( isStorageOptimized() == true && indicesAreAllocated() == false, std::logic_error, err );
    // if storage has been optimized, then number of allocated is now the number of entries
    TEUCHOS_TEST_FOR_EXCEPTION( isStorageOptimized() == true && getNodeAllocationSize() != getNodeNumEntries(), std::logic_error, err );
    // if graph doesn't have the global constants, then they should all be marked as invalid
    TEUCHOS_TEST_FOR_EXCEPTION( haveGlobalConstants_ == false && ( globalNumEntries_ != GSTI || globalNumDiags_ != GSTI || globalMaxNumRowEntries_ != GSTI ), std::logic_error, err );
    // if the graph has global cosntants, then they should be valid.
    TEUCHOS_TEST_FOR_EXCEPTION( haveGlobalConstants_ == true && ( globalNumEntries_ == GSTI || globalNumDiags_ == GSTI || globalMaxNumRowEntries_ == GSTI ), std::logic_error, err );
    TEUCHOS_TEST_FOR_EXCEPTION( haveGlobalConstants_ == true && ( globalNumEntries_ < nodeNumEntries_ || globalNumDiags_ < nodeNumDiags_ || globalMaxNumRowEntries_ < nodeMaxNumRowEntries_ ),
                        std::logic_error, err );
    // if indices are allocated, then the allocation specifications should have been released
    TEUCHOS_TEST_FOR_EXCEPTION( indicesAreAllocated() == true  && (numAllocForAllRows_ != 0 || numAllocPerRow_ != null),                        std::logic_error, err );
    // if indices are not allocated, then information dictating allocation quantities should be present
    TEUCHOS_TEST_FOR_EXCEPTION( indicesAreAllocated() == false && (nodeNumAllocated_ != STI || nodeNumEntries_ != 0),                           std::logic_error, err );
    // if storage is optimized, then profile should be static
    TEUCHOS_TEST_FOR_EXCEPTION( isStorageOptimized() && pftype_ != StaticProfile,                                                               std::logic_error, err );
    // rowPtrs_ is required to have N+1 rows; rowPtrs_[N] == gblInds1D_.size()/lclInds1D_.size()
    TEUCHOS_TEST_FOR_EXCEPTION( isGloballyIndexed() && rowPtrs_ != null && ((size_t)rowPtrs_.size() != getNodeNumRows()+1 || rowPtrs_[getNodeNumRows()] != (size_t)gblInds1D_.size()), std::logic_error, err );
    TEUCHOS_TEST_FOR_EXCEPTION(  isLocallyIndexed() && rowPtrs_ != null && ((size_t)rowPtrs_.size() != getNodeNumRows()+1 || rowPtrs_[getNodeNumRows()] != (size_t)lclInds1D_.size()), std::logic_error, err );
    // if profile is dynamic and we have allocated, then 2D allocations should be present
    TEUCHOS_TEST_FOR_EXCEPTION( pftype_ == DynamicProfile && indicesAreAllocated() && getNodeNumRows() > 0 && lclInds2D_ == null && gblInds2D_ == null,
                                                                                                                                        std::logic_error, err );
    // if profile is dynamic, then numentries and 2D indices are needed and should be present
    TEUCHOS_TEST_FOR_EXCEPTION( pftype_ == DynamicProfile && indicesAreAllocated() && getNodeNumRows() > 0 && (numRowEntries_ == null || (lclInds2D_ == null && gblInds2D_ == null)),
                                                                                                                                        std::logic_error, err );
    // if profile is dynamic, then 1D allocations should not be present
    TEUCHOS_TEST_FOR_EXCEPTION( pftype_ == DynamicProfile && (lclInds1D_ != null || gblInds1D_ != null),                                std::logic_error, err );
    // if profile is dynamic, then row offsets should not be present
    TEUCHOS_TEST_FOR_EXCEPTION( pftype_ == DynamicProfile && rowPtrs_ != null,                                                          std::logic_error, err );
    // if profile is static and we have allocated non-trivially, then 1D allocations should be present
    TEUCHOS_TEST_FOR_EXCEPTION( pftype_ == StaticProfile && indicesAreAllocated() && getNodeAllocationSize() > 0 && lclInds1D_ == null && gblInds1D_ == null,
                                                                                                                                        std::logic_error, err );
    // if profile is static and we have a non-trivial application, then row offsets should be allocated
    TEUCHOS_TEST_FOR_EXCEPTION( pftype_ == StaticProfile && indicesAreAllocated() && getNodeNumRows() > 0 && (rowPtrs_ == null),
                                                                                                                                        std::logic_error, err );
    // if profile is static, then 2D allocations should not be present
    TEUCHOS_TEST_FOR_EXCEPTION( pftype_ == StaticProfile && (lclInds2D_ != null || gblInds2D_ != null),                                         std::logic_error, err );
    // if indices are not allocated, then none of the buffers should be.
    TEUCHOS_TEST_FOR_EXCEPTION( indicesAreAllocated() == false && (rowPtrs_ != null || numRowEntries_ != null ||
                                                           lclInds1D_ != null || lclInds2D_ != null ||
                                                           gblInds1D_ != null || gblInds2D_ != null),                                   std::logic_error, err );
    // indices may be local or global only if they are allocated (numAllocated is redundant; could simply be indicesAreLocal_ || indicesAreGlobal_)
    TEUCHOS_TEST_FOR_EXCEPTION( (indicesAreLocal_ == true || indicesAreGlobal_ == true) && indicesAreAllocated_ == false,         std::logic_error, err );
    // indices may be local or global, but not both
    TEUCHOS_TEST_FOR_EXCEPTION( indicesAreLocal_ == true && indicesAreGlobal_ == true,                                            std::logic_error, err );
    // if indices are local, then global allocations should not be present
    TEUCHOS_TEST_FOR_EXCEPTION( indicesAreLocal_ == true && (gblInds1D_ != null || gblInds2D_ != null),                           std::logic_error, err );
    // if indices are global, then local allocations should not be present
    TEUCHOS_TEST_FOR_EXCEPTION( indicesAreGlobal_ == true && (lclInds1D_ != null || lclInds2D_ != null),                          std::logic_error, err );
    // if indices are local, then local allocations should be present
    TEUCHOS_TEST_FOR_EXCEPTION( indicesAreLocal_ == true && getNodeAllocationSize() > 0 && lclInds1D_ == null && getNodeNumRows() > 0 && lclInds2D_ == null,
                                                                                                                          std::logic_error, err );
    // if indices are global, then global allocations should be present
    TEUCHOS_TEST_FOR_EXCEPTION( indicesAreGlobal_ == true && getNodeAllocationSize() > 0 && gblInds1D_ == null && getNodeNumRows() > 0 && gblInds2D_ == null,
                                                                                                                          std::logic_error, err );
    // if indices are allocated, then we should have recorded how many were allocated
    TEUCHOS_TEST_FOR_EXCEPTION( indicesAreAllocated() == true  && nodeNumAllocated_ == STI,                                       std::logic_error, err );
    // if indices are not allocated, then the allocation size should be marked invalid
    TEUCHOS_TEST_FOR_EXCEPTION( indicesAreAllocated() == false && nodeNumAllocated_ != STI,                                       std::logic_error, err );
    // check the actual allocations
    if (indicesAreAllocated()) {
      size_t actualNumAllocated = 0;
      if (pftype_ == DynamicProfile) {
        if (isGloballyIndexed() && gblInds2D_ != null) {
          for (size_t r = 0; r < getNodeNumRows(); ++r) {
            actualNumAllocated += gblInds2D_[r].size();
          }
        }
        else if (isLocallyIndexed() && lclInds2D_ != null) {
          for (size_t r = 0; r < getNodeNumRows(); ++r) {
            actualNumAllocated += lclInds2D_[r].size();
          }
        }
      }
      else { // pftype_ == StaticProfile)
        actualNumAllocated = rowPtrs_[getNodeNumRows()];
        TEUCHOS_TEST_FOR_EXCEPTION(  isLocallyIndexed() == true && (size_t)lclInds1D_.size() != actualNumAllocated, std::logic_error, err );
        TEUCHOS_TEST_FOR_EXCEPTION( isGloballyIndexed() == true && (size_t)gblInds1D_.size() != actualNumAllocated, std::logic_error, err );
      }
      TEUCHOS_TEST_FOR_EXCEPTION(actualNumAllocated != nodeNumAllocated_, std::logic_error, err );
    }
#endif
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                  User-visible class methods                             //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const
  {
    const LocalOrdinal lrow = rowMap_->getLocalElement(globalRow);
    size_t ret = OrdinalTraits<size_t>::invalid();
    if (lrow != OrdinalTraits<LocalOrdinal>::invalid())
    {
      RowInfo rowinfo = getRowInfo(lrow);
      ret = rowinfo.numEntries;
    }
    return ret;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getNumEntriesInLocalRow(LocalOrdinal localRow) const
  {
    size_t ret = OrdinalTraits<size_t>::invalid();
    if (rowMap_->isNodeLocalElement(localRow)) {
      RowInfo rowinfo = getRowInfo(localRow);
      ret = rowinfo.numEntries;
    }
    return ret;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getNumAllocatedEntriesInGlobalRow(GlobalOrdinal globalRow) const
  {
    const LocalOrdinal lrow = rowMap_->getLocalElement(globalRow);
    size_t ret = OrdinalTraits<size_t>::invalid();
    if (lrow != OrdinalTraits<LocalOrdinal>::invalid())
    {
      RowInfo rowinfo = getRowInfo(lrow);
      ret = rowinfo.allocSize;
    }
    return ret;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getNumAllocatedEntriesInLocalRow(LocalOrdinal localRow) const
  {
    size_t ret = OrdinalTraits<size_t>::invalid();
    if (rowMap_->isNodeLocalElement(localRow)) {
      RowInfo rowinfo = getRowInfo(localRow);
      ret = rowinfo.allocSize;
    }
    return ret;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  ArrayRCP<const size_t> CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getNodeRowPtrs() const
  {
    return rowPtrs_;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  ArrayRCP<const LocalOrdinal> CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getNodePackedIndices() const
  {
    return lclInds1D_;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getLocalRowCopy(LocalOrdinal localRow, const ArrayView<LocalOrdinal> &indices, size_t &NumIndices) const
  {
    // can only do this if
    // * we have local indices: isLocallyIndexed()
    // or
    // * we are capable of producing them: isGloballyIndexed() && hasColMap()
    // short circuit if we aren't allocated
    std::string tfecfFuncName("getLocalRowCopy(localRow,...)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(isGloballyIndexed() == true && hasColMap() == false, std::runtime_error, ": local indices cannot be produced.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(rowMap_->isNodeLocalElement(localRow) == false, std::runtime_error,
        ": localRow (== " << localRow << ") is not valid on this node.");
    const RowInfo rowinfo = getRowInfo(localRow);
    NumIndices = rowinfo.numEntries;
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC((size_t)indices.size() < NumIndices, std::runtime_error,
        ": specified storage (size==" << indices.size()
        << ") is not large enough to hold all entries for this row (NumIndices == " << NumIndices << ").");
    if (isLocallyIndexed()) {
      ArrayView<const LocalOrdinal> lview = getLocalView(rowinfo);
      std::copy( lview.begin(), lview.begin() + NumIndices, indices.begin());
    }
    else if (isGloballyIndexed()) {
      ArrayView<const GlobalOrdinal> gview = getGlobalView(rowinfo);
      for (size_t j=0; j < NumIndices; ++j) {
        indices[j] = colMap_->getLocalElement(gview[j]);
      }
    }
    else {
#ifdef HAVE_TPETRA_DEBUG
      // should have fallen in one of the above
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( indicesAreAllocated() == true, std::logic_error, ": Internal logic error. Please contact Tpetra team.");
#endif
      NumIndices = 0;
    }
    return;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getGlobalRowCopy(GlobalOrdinal globalRow, const ArrayView<GlobalOrdinal> &indices, size_t &NumIndices) const
  {
    // we either currently store global indices, or we have a column map with which to transcribe our local indices for the user
    const LocalOrdinal lrow = rowMap_->getLocalElement(globalRow);
    std::string tfecfFuncName("getGlobalRowCopy(globalRow,...)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(lrow == OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error,
        ": globalRow (== " << globalRow << ") does not belong to this node.");
    const RowInfo rowinfo = getRowInfo((size_t)lrow);
    NumIndices = rowinfo.numEntries;
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC((size_t)indices.size() < NumIndices, std::runtime_error,
        ": specified storage (size==" << indices.size()
        << ") is not large enough to hold all entries for this row (NumIndices == " << NumIndices << ").");
    if (isLocallyIndexed()) {
      ArrayView<const LocalOrdinal> lview = getLocalView(rowinfo);
      for (size_t j=0; j < NumIndices; ++j) {
        indices[j] = colMap_->getGlobalElement(lview[j]);
      }
    }
    else if (isGloballyIndexed()) {
      ArrayView<const GlobalOrdinal> gview = getGlobalView(rowinfo);
      std::copy(gview.begin(), gview.begin() + NumIndices, indices.begin());
    }
    return;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getLocalRowView(LocalOrdinal localRow, ArrayView<const LocalOrdinal> &indices) const
  {
    std::string tfecfFuncName("getLocalRowView()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(isGloballyIndexed() == true, std::runtime_error, ": local indices cannot be provided.");
    indices = null;
    if (rowMap_->isNodeLocalElement(localRow) == true) {
      const RowInfo rowinfo = getRowInfo(localRow);
      if (rowinfo.numEntries > 0) {
        indices = getLocalView(rowinfo);
        indices = indices(0,rowinfo.numEntries);
      }
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( (size_t)indices.size() != getNumEntriesInLocalRow(localRow), std::logic_error, ": Violated stated post-conditions. Please contact Tpetra team.");
#endif
    return;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getGlobalRowView(GlobalOrdinal globalRow, ArrayView<const GlobalOrdinal> &indices) const
  {
    std::string tfecfFuncName("getGlobalRowView()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(isLocallyIndexed() == true, std::runtime_error, ": global indices cannot be provided.");
    indices = null;
    if (rowMap_->isNodeGlobalElement(globalRow) == true) {
      const LocalOrdinal lrow = rowMap_->getLocalElement(globalRow);
      const RowInfo rowinfo = getRowInfo(lrow);
      if (rowinfo.numEntries > 0) {
        indices = getGlobalView(rowinfo);
        indices = indices(0,rowinfo.numEntries);
      }
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( (size_t)indices.size() != getNumEntriesInGlobalRow(globalRow), std::logic_error, ": Violated stated post-conditions. Please contact Tpetra team.");
#endif
    return;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::insertLocalIndices(
            LocalOrdinal localRow, const ArrayView<const LocalOrdinal> &indices)
  {
    std::string tfecfFuncName("insertLocalIndices()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( isFillActive() == false,                        std::runtime_error, ": requires that fill is active.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( isGloballyIndexed() == true,                    std::runtime_error, ": graph indices are global; use insertGlobalIndices().");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( hasColMap() == false,                           std::runtime_error, ": cannot insert local indices without a column map.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( rowMap_->isNodeLocalElement(localRow) == false, std::runtime_error, ": row does not belong to this node.");
    if (indicesAreAllocated() == false) {
      allocateIndices(LocalIndices);
    }
    // use column map to filter the entries
    Array<LocalOrdinal> f_inds(indices);
    SLocalGlobalNCViews inds_ncview;
    inds_ncview.linds = f_inds();
    const size_t numFilteredEntries = filterIndices<LocalIndices>(inds_ncview);
    if (numFilteredEntries > 0) {
      RowInfo rowInfo = getRowInfo(localRow);
      const size_t curNumEntries = rowInfo.numEntries;
      const size_t newNumEntries = curNumEntries + numFilteredEntries;
      if (newNumEntries > rowInfo.allocSize) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(getProfileType() == StaticProfile, std::runtime_error, ": new indices exceed statically allocated graph structure.");
        TPETRA_EFFICIENCY_WARNING(true,std::runtime_error,
            "::insertLocalIndices(): Pre-allocated space has been exceeded, requiring new allocation. To improve efficiency, suggest larger allocation.");
        // update allocation only as much as necessary
        rowInfo = updateAlloc<LocalIndices>(rowInfo, newNumEntries);
      }
      SLocalGlobalViews inds_view;
      inds_view.linds = f_inds(0,numFilteredEntries);
      insertIndices<LocalIndices,LocalIndices>(rowInfo, inds_view);
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(indicesAreAllocated() == false || isLocallyIndexed() == false, std::logic_error,
        ": Violated stated post-conditions. Please contact Tpetra team.");
#endif
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::insertGlobalIndices(
            GlobalOrdinal grow, const ArrayView<const GlobalOrdinal> &indices)
  {
    std::string tfecfFuncName("insertGlobalIndices()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(isLocallyIndexed() == true, std::runtime_error, ": graph indices are local; use insertLocalIndices().");
    // this can't really be satisfied for now, because if we are fillComplete(), then we are local
    // in the future, this may change. however, the rule that modification require active fill will not change.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( isFillActive() == false, std::runtime_error, ": requires that fill is active.");
    if (indicesAreAllocated() == false) {
      allocateIndices(GlobalIndices);
    }
    const LocalOrdinal myRow = rowMap_->getLocalElement(grow);
    if (myRow != OrdinalTraits<LocalOrdinal>::invalid())
    {
      // if we have a column map, use it to filter the entries.
      Array<GlobalOrdinal> filtered_indices;
      SLocalGlobalViews inds_view;
      if (hasColMap()) {
        SLocalGlobalNCViews inds_ncview;
        // filter indices and values through the column map
        filtered_indices.assign(indices.begin(), indices.end());
        inds_ncview.ginds = filtered_indices();
        const size_t numFilteredEntries = filterIndices<GlobalIndices>(inds_ncview);
        inds_view.ginds = filtered_indices(0,numFilteredEntries);
      }
      else {
        inds_view.ginds = indices;
      }
      const size_t numFilteredEntries = inds_view.ginds.size();
      // add the new indices and values
      if (numFilteredEntries > 0) {
        RowInfo rowInfo = getRowInfo(myRow);
        const size_t curNumEntries = rowInfo.numEntries;
        const size_t newNumEntries = curNumEntries + numFilteredEntries;
        if (newNumEntries > rowInfo.allocSize) {
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(getProfileType() == StaticProfile, std::runtime_error, ": new indices exceed statically allocated graph structure.");
          TPETRA_EFFICIENCY_WARNING(true,std::runtime_error,
              "::insertGlobalValues(): Pre-allocated space has been exceeded, requiring new allocation. To improve efficiency, suggest larger allocation.");
          // update allocation only as much as necessary
          rowInfo = updateAlloc<GlobalIndices>(rowInfo, newNumEntries);
        }
        insertIndices<GlobalIndices,GlobalIndices>(rowInfo, inds_view);
#ifdef HAVE_TPETRA_DEBUG
        {
          const size_t chkNewNumEntries = getNumEntriesInLocalRow(myRow);
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(chkNewNumEntries != newNumEntries, std::logic_error, ": Internal logic error. Please contact Tpetra team.");
        }
#endif
      }
    }
    else {
      // a nonlocal row
      for (typename ArrayView<const GlobalOrdinal>::iterator i=indices.begin(); i != indices.end(); ++i) {
        nonlocals_[grow].push_back(*i);
      }
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(indicesAreAllocated() == false || isGloballyIndexed() == false, std::logic_error,
        ": Violated stated post-conditions. Please contact Tpetra team.");
#endif
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::removeLocalIndices(LocalOrdinal lrow)
  {
    std::string tfecfFuncName("removeLocalIndices()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( isFillActive() == false,                    std::runtime_error, ": requires that fill is active.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( isStorageOptimized() == true,               std::runtime_error, ": cannot remove indices after optimizeStorage() has been called.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( isGloballyIndexed() == true,                std::runtime_error, ": graph indices are global.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( rowMap_->isNodeLocalElement(lrow) == false, std::runtime_error, ": row does not belong to this node.");
    if (indicesAreAllocated() == false) {
      allocateIndices(LocalIndices);
    }
    //
    clearGlobalConstants();
    //
    if (numRowEntries_ != null) {
      nodeNumEntries_ -= numRowEntries_[lrow];
      numRowEntries_[lrow] = 0;
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(getNumEntriesInLocalRow(lrow) != 0 || indicesAreAllocated() == false || isLocallyIndexed() == false, std::logic_error,
        ": Violated stated post-conditions. Please contact Tpetra team.");
#endif
  }


  // TODO: in the future, globalAssemble() should use import/export functionality
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::globalAssemble()
  {
    using Teuchos::as;
    using Teuchos::gatherAll;
    using Teuchos::ireceive;
    using Teuchos::isend;
    using Teuchos::outArg;
    using Teuchos::REDUCE_MAX;
    using Teuchos::reduceAll;
    using Teuchos::waitAll;
    using std::deque;
    using std::pair;
    using std::make_pair;
    typedef GlobalOrdinal GO;
    typedef typename std::map<GO, std::deque<GO> >::const_iterator NLITER;

    std::string tfecfFuncName("globalAssemble()"); // for exception macro
    RCP<const Comm<int> > comm = getComm();

    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
#ifdef HAVE_TPETRA_DEBUG
    Teuchos::barrier (*comm);
#endif // HAVE_TPETRA_DEBUG

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( ! isFillActive(), std::runtime_error,
      ": requires that fill is active.");
    // Determine if any nodes have global entries to share
    {
      size_t MyNonlocals = nonlocals_.size(), MaxGlobalNonlocals;
      reduceAll (*comm, REDUCE_MAX, MyNonlocals, outArg (MaxGlobalNonlocals));
      if (MaxGlobalNonlocals == 0) {
        return;  // no entries to share
      }
    }

    // compute a list of NLRs from nonlocals_ and use it to compute:
    //      IdsAndRows: a vector of (id,row) pairs
    //          NLR2Id: a map from NLR to the Id that owns it
    // globalNeighbors: a global graph of connectivity between images: globalNeighbors(i,j) indicates that j sends to i
    //         sendIDs: a list of all images I send to
    //         recvIDs: a list of all images I receive from (constructed later)
    Array<pair<int, GO> > IdsAndRows;
    std::map<GO, int> NLR2Id;
    Teuchos::SerialDenseMatrix<int, char> globalNeighbors;
    Array<int> sendIDs, recvIDs;
    {
      // nonlocals_ contains the entries we are holding for all non-local rows
      // we want a list of the rows for which we have data
      Array<GO> NLRs;
      std::set<GO> setOfRows;
      for (NLITER iter = nonlocals_.begin(); iter != nonlocals_.end(); ++iter) {
        setOfRows.insert(iter->first);
      }
      // copy the elements in the set into an Array
      NLRs.resize(setOfRows.size());
      std::copy(setOfRows.begin(), setOfRows.end(), NLRs.begin());

      // get a list of ImageIDs for the non-local rows (NLRs)
      Array<int> NLRIds(NLRs.size());
      {
        LookupStatus stat = rowMap_->getRemoteIndexList(NLRs(),NLRIds());
        char lclerror = ( stat == IDNotPresent ? 1 : 0 );
        char gblerror;
        reduceAll (*getComm(), REDUCE_MAX, lclerror, outArg (gblerror));
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(gblerror != 0, std::runtime_error,
          ": nonlocal entries correspond to invalid rows.");
      }

      // build up a list of neighbors, as well as a map between NLRs and Ids
      // localNeighbors[i] != 0 iff I have data to send to image i
      // put NLRs,Ids into an array of pairs
      IdsAndRows.reserve(NLRs.size());
      Array<char> localNeighbors(numImages, 0);
      typename Array<GO>::const_iterator nlr;
      typename Array<int>::const_iterator id;
      for (nlr = NLRs.begin(), id = NLRIds.begin();
           nlr != NLRs.end(); ++nlr, ++id) {
        NLR2Id[*nlr] = *id;
        localNeighbors[*id] = 1;
        // IdsAndRows.push_back(make_pair<int,GlobalOrdinal>(*id,*nlr));
        IdsAndRows.push_back(make_pair(*id,*nlr));
      }
      for (int j=0; j<numImages; ++j) {
        if (localNeighbors[j]) {
          sendIDs.push_back(j);
        }
      }
      // sort IdsAndRows, by Ids first, then rows
      std::sort(IdsAndRows.begin(),IdsAndRows.end());
      // gather from other nodes to form the full graph
      globalNeighbors.shapeUninitialized(numImages,numImages);
      gatherAll (*getComm(), numImages, localNeighbors.getRawPtr(),
                 numImages * numImages, globalNeighbors.values());
      // globalNeighbors at this point contains (on all images) the
      // connectivity between the images.
      // globalNeighbors(i,j) != 0 means that j sends to i/that i receives from j
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // FIGURE OUT WHO IS SENDING TO WHOM AND HOW MUCH
    // DO THIS IN THE PROCESS OF PACKING ALL OUTGOING DATA ACCORDING TO DESTINATION ID
    //////////////////////////////////////////////////////////////////////////////////////

    // loop over all columns to know from which images I can expect to receive something
    for (int j=0; j<numImages; ++j) {
      if (globalNeighbors(myImageID,j)) {
        recvIDs.push_back(j);
      }
    }
    const size_t numRecvs = recvIDs.size();

    // we know how many we're sending to already
    // form a contiguous list of all data to be sent
    // track the number of entries for each ID
    Array<pair<GO, GO> > IJSendBuffer;
    Array<size_t> sendSizes(sendIDs.size(), 0);
    size_t numSends = 0;
    for (typename Array<pair<int, GO> >::const_iterator IdAndRow = IdsAndRows.begin();
         IdAndRow != IdsAndRows.end(); ++IdAndRow) {
      int id = IdAndRow->first;
      GO row = IdAndRow->second;
      // have we advanced to a new send?
      if (sendIDs[numSends] != id) {
        numSends++;
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(sendIDs[numSends] != id,
          std::logic_error, ": internal logic error. Contact Tpetra team.");
      }
      // copy data for row into contiguous storage
      for (typename deque<GO>::const_iterator j = nonlocals_[row].begin(); j != nonlocals_[row].end(); ++j)
      {
        IJSendBuffer.push_back( pair<GlobalOrdinal,GlobalOrdinal>(row,*j) );
        sendSizes[numSends]++;
      }
    }
    if (IdsAndRows.size() > 0) {
      numSends++; // one last increment, to make it a count instead of an index
    }
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(as<typename Array<int>::size_type>(numSends) != sendIDs.size(), std::logic_error, ": internal logic error. Contact Tpetra team.");

    // don't need this data anymore
    nonlocals_.clear();

    //////////////////////////////////////////////////////////////////////////////////////
    // TRANSMIT SIZE INFO BETWEEN SENDERS AND RECEIVERS
    //////////////////////////////////////////////////////////////////////////////////////

    // Array of pending nonblocking communication requests.  It's OK
    // to mix nonblocking send and receive requests in the same
    // waitAll() call.
    Array<RCP<Teuchos::CommRequest> > requests;

    // perform non-blocking sends: send sizes to our recipients
    for (size_t s = 0; s < numSends ; ++s) {
      // We're using a nonowning RCP because all communication
      // will be local to this method and the scope of our data
      requests.push_back (isend<int, size_t> (*comm,
                                              rcp (&sendSizes[s], false),
                                              sendIDs[s]));
    }
    // perform non-blocking receives: receive sizes from our senders
    Array<size_t> recvSizes (numRecvs);
    for (size_t r = 0; r < numRecvs; ++r) {
      // We're using a nonowning RCP because all communication
      // will be local to this method and the scope of our data
      requests.push_back (ireceive (*comm, rcp (&recvSizes[r], false), recvIDs[r]));
    }
    // Wait on all the nonblocking sends and receives.
    if (! requests.empty()) {
      waitAll (*comm, requests());
    }
#ifdef HAVE_TPETRA_DEBUG
    Teuchos::barrier (*comm);
#endif // HAVE_TPETRA_DEBUG

    // This doesn't necessarily deallocate the array.
    requests.resize (0);

    ////////////////////////////////////////////////////////////////////////////////////
    // NOW SEND/RECEIVE ALL ROW DATA
    ////////////////////////////////////////////////////////////////////////////////////
    // from the size info, build the ArrayViews into IJSendBuffer
    Array<ArrayView<pair<GO,GO> > > sendBuffers(numSends,null);
    {
      size_t cur = 0;
      for (size_t s=0; s<numSends; ++s) {
        sendBuffers[s] = IJSendBuffer(cur,sendSizes[s]);
        cur += sendSizes[s];
      }
    }
    // perform non-blocking sends
    for (size_t s=0; s < numSends ; ++s)
    {
      // We're using a nonowning RCP because all communication
      // will be local to this method and the scope of our data
      ArrayRCP<pair<GO,GO> > tmpSendBuf =
        arcp (sendBuffers[s].getRawPtr(), 0, sendBuffers[s].size(), false);
      requests.push_back (isend<int, pair<GO,GO> > (*comm, tmpSendBuf, sendIDs[s]));
    }
    // calculate amount of storage needed for receives
    // setup pointers for the receives as well
    size_t totalRecvSize = std::accumulate (recvSizes.begin(), recvSizes.end(), 0);
    Array<pair<GO,GO> > IJRecvBuffer (totalRecvSize);
    // from the size info, build the ArrayViews into IJRecvBuffer
    Array<ArrayView<pair<GO,GO> > > recvBuffers (numRecvs, null);
    {
      size_t cur = 0;
      for (size_t r=0; r<numRecvs; ++r) {
        recvBuffers[r] = IJRecvBuffer(cur,recvSizes[r]);
        cur += recvSizes[r];
      }
    }
    // perform non-blocking recvs
    for (size_t r = 0; r < numRecvs; ++r) {
      // We're using a nonowning RCP because all communication
      // will be local to this method and the scope of our data
      ArrayRCP<pair<GO,GO> > tmpRecvBuf =
        arcp (recvBuffers[r].getRawPtr(), 0, recvBuffers[r].size(), false);
      requests.push_back (ireceive (*comm, tmpRecvBuf, recvIDs[r]));
    }
    // perform waits
    if (! requests.empty()) {
      waitAll (*comm, requests());
    }
#ifdef HAVE_TPETRA_DEBUG
    Teuchos::barrier (*comm);
#endif // HAVE_TPETRA_DEBUG

    ////////////////////////////////////////////////////////////////////////////////////
    // NOW PROCESS THE RECEIVED ROW DATA
    ////////////////////////////////////////////////////////////////////////////////////
    // TODO: instead of adding one entry at a time, add one row at a time.
    //       this requires resorting; they arrived sorted by sending node, so that entries could be non-contiguous if we received
    //       multiple entries for a particular row from different processors.
    //       it also requires restoring the data, which may make it not worth the trouble.
    for (typename Array<pair<GO,GO> >::const_iterator ij = IJRecvBuffer.begin();
         ij != IJRecvBuffer.end(); ++ij)
    {
      insertGlobalIndices(ij->first, tuple<GO> (ij->second));
    }
    checkInternalState();
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::resumeFill(const RCP<ParameterList> &params)
  {
#ifdef HAVE_TPETRA_DEBUG
    Teuchos::barrier( *rowMap_->getComm() );
#endif
    clearGlobalConstants();
    lclGraph_ = null;
    if (params != null) this->setParameterList (params);
    setSorted(true);
    lowerTriangular_  = false;
    upperTriangular_  = false;
    noRedundancies_   = false;
    fillComplete_ = false;
#ifdef HAVE_TPETRA_DEBUG
    std::string tfecfFuncName("resumeFill()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( isFillActive() == false || isFillComplete() == true, std::logic_error, ": Violated stated post-conditions. Please contact Tpetra team.");
#endif
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::fillComplete(const RCP<ParameterList> &params)
  {
    fillComplete(rowMap_,rowMap_,params);
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::fillComplete(
                                    const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap,
                                    const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap,
                                    const RCP<ParameterList> &params)
  {
#ifdef HAVE_TPETRA_DEBUG
    Teuchos::barrier( *rowMap_->getComm() );
#endif // HAVE_TPETRA_DEBUG
    std::string tfecfFuncName("fillComplete()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( ! isFillActive() || isFillComplete(),
      std::runtime_error, ": Graph fill state must be active (isFillActive() "
      "must be true) before calling fillComplete().");

    // allocate if unallocated
    if (! indicesAreAllocated()) {
      // allocate global, in case we do not have a column map
      allocateIndices( GlobalIndices );
    }
    // Global assemble, if we need to (we certainly don't need to if
    // there's only one process).  This call only costs a single
    // all-reduce if we don't need global assembly.
    if (getComm()->getSize() > 1) {
      globalAssemble();
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          nonlocals_.size() > 0, 
          std::runtime_error, ": cannot have non-local entries on a serial run. Invalid entries were submitted to the CrsMatrix.");
    }
    // set domain/range map: may clear the import/export objects
    setDomainRangeMaps(domainMap,rangeMap);
    // make column map
    if (! hasColMap()) {
      makeColMap();
    }
    if (isGloballyIndexed()) {
      makeIndicesLocal();
    }
    if (! isSorted()) {
      sortAllIndices();
    }
    if (! isMerged()) {
      mergeAllIndices();
    }
    makeImportExport(); // Make Import and Export objects
    computeGlobalConstants();
    // fill local objects
    fillLocalGraph(params);
    //
    fillComplete_ = true;
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( isFillActive() == true || isFillComplete() == false, std::logic_error, ": Violated stated post-conditions. Please contact Tpetra team.");
#endif
    //
    checkInternalState();
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::fillLocalGraph(const RCP<ParameterList> &params)
  {
    const size_t numRows = getNodeNumRows();
    ArrayRCP<LocalOrdinal> inds;
    ArrayRCP<size_t>       ptrs;
    bool requestOptimizedStorage = true;
    if (params != null && params->get("Optimize Storage",true) == false) requestOptimizedStorage = false;
    if (getProfileType() == DynamicProfile) {
      // 2d -> 1d packed
      ptrs = LocalMatOps::allocRowPtrs( getRowMap()->getNode(), numRowEntries_() );
      inds = LocalMatOps::template allocStorage<LocalOrdinal>( getRowMap()->getNode(), ptrs() );
      for (size_t row=0; row < numRows; ++row) {
        const size_t numentrs = numRowEntries_[row];
        std::copy( lclInds2D_[row].begin(), lclInds2D_[row].begin()+numentrs, inds+ptrs[row] );
      }
    }
    else if (getProfileType() == StaticProfile) {
      // 1d non-packed -> 1d packed
      if (nodeNumEntries_ != nodeNumAllocated_) {
        ptrs = LocalMatOps::allocRowPtrs( getRowMap()->getNode(), numRowEntries_() );
        inds = LocalMatOps::template allocStorage<LocalOrdinal>( getRowMap()->getNode(), ptrs() );
        for (size_t row=0; row < numRows; ++row) {
          const size_t numentrs = numRowEntries_[row];
          std::copy( lclInds1D_+rowPtrs_[row], lclInds1D_+rowPtrs_[row]+numentrs, inds+ptrs[row] );
        }
      }
      else {
        inds = lclInds1D_;
        ptrs = rowPtrs_;
      }
    }
    // can we ditch the old allocations for the packed one?
    if ( requestOptimizedStorage ) {
      lclInds2D_ = null;
      numRowEntries_ = null;
      // keep the new stuff
      lclInds1D_ = inds;
      rowPtrs_ = ptrs;
      nodeNumAllocated_ = nodeNumEntries_;
      pftype_ = StaticProfile;
    }
    // build the local graph, hand over the indices
    RCP<ParameterList> lclparams;
    if (params == null) lclparams = parameterList();
    else                lclparams = sublist(params,"Local Graph");
    lclGraph_ = rcp( new local_graph_type( getRowMap()->getNodeNumElements(), getRowMap()->getNode(), lclparams ) );
    lclGraph_->setStructure(ptrs,inds);
    ptrs = null;
    inds = null;
    // finalize local graph
    Teuchos::EDiag diag = ( getNodeNumDiags() < getNodeNumRows() ? Teuchos::UNIT_DIAG : Teuchos::NON_UNIT_DIAG );
    Teuchos::EUplo uplo = Teuchos::UNDEF_TRI;
    if      (isUpperTriangular()) uplo = Teuchos::UPPER_TRI;
    else if (isLowerTriangular()) uplo = Teuchos::LOWER_TRI;
    LocalMatOps::finalizeGraph(uplo,diag,*lclGraph_,params);
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  const RCP<const typename LocalMatOps::template graph<LocalOrdinal,Node>::graph_type> 
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getLocalGraph() const
  {
    return lclGraph_;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  const RCP<typename LocalMatOps::template graph<LocalOrdinal,Node>::graph_type> 
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getLocalGraphNonConst()
  {
    return lclGraph_;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::computeGlobalConstants()
  {
    // compute the local constants first
    const size_t nlrs = getNodeNumRows();
    // reset all local properties
    upperTriangular_ = true;
    lowerTriangular_ = true;
    nodeMaxNumRowEntries_ = 0;
    nodeNumDiags_         = 0;
    // indices are already sorted in each row
    const Map<LocalOrdinal,GlobalOrdinal,Node> &rowMap = *rowMap_;
    if (indicesAreAllocated() == true && nodeNumAllocated_ > 0) {
      for (size_t r=0; r < nlrs; ++r) {
        GlobalOrdinal rgid = rowMap.getGlobalElement(r);
        // determine the local column index for this row, used for delimiting the diagonal
        const LocalOrdinal rlcid = colMap_->getLocalElement(rgid);
        RowInfo rowinfo = getRowInfo(r);
        ArrayView<const LocalOrdinal> rview = getLocalView(rowinfo);
        typename ArrayRCP<const LocalOrdinal>::iterator beg, end, cur;
        beg = rview.begin();
        end = beg + rowinfo.numEntries;
        if (beg != end) {
          for (cur = beg; cur != end; ++cur) {
            // is this the diagonal?
            if (rlcid == (*cur)) ++nodeNumDiags_;
          }
          // because of sorting, smallest column index is (*beg); it indicates upper triangularity
          if (Teuchos::as<size_t>(beg[0]) < r) upperTriangular_ = false;
          // because of sorting, largest column index is (*newend); it indicates lower triangularity
          if (r < Teuchos::as<size_t>(end[-1])) lowerTriangular_ = false;
        }
        // compute num entries for this row, accumulate into nodeNumEntries_, update nodeMaxNumRowEntries_
        nodeMaxNumRowEntries_ = std::max( nodeMaxNumRowEntries_, rowinfo.numEntries );
      }
    }

    // compute global constants using computed local constants
    if (haveGlobalConstants_ == false) {
      global_size_t lcl[2], gbl[2];
      lcl[0] = nodeNumEntries_;
      lcl[1] = nodeNumDiags_;
      Teuchos::reduceAll<int,global_size_t>(*getComm(),Teuchos::REDUCE_SUM,2,lcl,gbl);
      globalNumEntries_ = gbl[0];
      globalNumDiags_   = gbl[1];
      Teuchos::reduceAll<int,global_size_t>(*getComm(),Teuchos::REDUCE_MAX,nodeMaxNumRowEntries_,
          outArg(globalMaxNumRowEntries_));
      haveGlobalConstants_ = true;
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::makeIndicesLocal()
  {
    // All nodes must be in the same index state.
    // Update index state by checking isLocallyIndexed/Global on all nodes
    computeIndexState();
    std::string tfecfFuncName("makeIndicesLocal()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(isLocallyIndexed() && isGloballyIndexed(), std::logic_error,
        ": inconsistent index state. Indices must be local on all nodes or global on all nodes.");
    // If user has not prescribed column map, create one from indices
    makeColMap();
    // Transform indices to local index space
    const size_t nlrs = getNodeNumRows();
    //
    if (isGloballyIndexed() && nlrs > 0) {
      // allocate data for local indices
      if (getProfileType() == StaticProfile) {
        // reinterpret the the compute buffer as LocalOrdinal if they are the same size
        // otherwise, just reallocate
        if (nodeNumAllocated_ && sizeof(LocalOrdinal) == sizeof(GlobalOrdinal) ) {
          lclInds1D_ = arcp_reinterpret_cast<LocalOrdinal>(gblInds1D_).persistingView(0,nodeNumAllocated_);
        }
        else {
          lclInds1D_ = LocalMatOps::template allocStorage<LocalOrdinal>( getRowMap()->getNode(), rowPtrs_() );
        }
        for (size_t r=0; r < getNodeNumRows(); ++r) {
          const size_t offset   = rowPtrs_[r],
                       numentry = numRowEntries_[r];
          for (size_t j=0; j<numentry; ++j) {
            GlobalOrdinal gid = gblInds1D_[offset + j];
            LocalOrdinal  lid = colMap_->getLocalElement(gid);
            lclInds1D_[offset + j] = lid;
#ifdef HAVE_TPETRA_DEBUG
            TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(lclInds1D_[offset + j] == OrdinalTraits<LocalOrdinal>::invalid(), std::logic_error,
                ": Internal error in fillComplete(). Please contact Tpetra team.");
#endif
          }
        }
        // done with this pointer (allocation will persist in lclInds1D_)
        gblInds1D_ = null;
      }
      else {  // getProfileType() == DynamicProfile
        lclInds2D_ = arcp< ArrayRCP<LocalOrdinal> >(nlrs);
        // if we have views, then make views
        for (size_t r=0; r < getNodeNumRows(); ++r) {
          if (gblInds2D_[r] != null) {
            ArrayRCP<GlobalOrdinal> ginds = gblInds2D_[r];
            const size_t rna = gblInds2D_[r].size();
            ArrayRCP< LocalOrdinal> linds = arcp_reinterpret_cast<LocalOrdinal>(ginds).persistingView(0,rna);
            // do the conversion in situ. this must be done from front to back.
            const size_t numentry = numRowEntries_[r];
            for (size_t j=0; j < numentry; ++j) {
              GlobalOrdinal gid = ginds[j];
              LocalOrdinal  lid = colMap_->getLocalElement(gid);
              linds[j] = lid;
#ifdef HAVE_TPETRA_DEBUG
              TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(linds[j] == OrdinalTraits<LocalOrdinal>::invalid(), std::logic_error,
                  ": Internal error in makeIndicesLocal(). Please contact Tpetra team.");
#endif
            }
            // reinterpret the the compute buffer as LocalOrdinal; keep the original size
            lclInds2D_[r] = linds;
            gblInds2D_[r] = null;
          }
        }
        gblInds2D_ = null;
      }
    }
    indicesAreLocal_  = true;
    indicesAreGlobal_ = false;
    checkInternalState();
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::computeIndexState()
  {
    char myIndices[2] = {0,0};
    if (indicesAreLocal_)  myIndices[0] = 1;
    if (indicesAreGlobal_) myIndices[1] = 1;
    char allIndices[2];
    Teuchos::reduceAll(*getComm(),Teuchos::REDUCE_MAX,2,myIndices,allIndices);
    indicesAreLocal_  = (allIndices[0]==1);  // If indices are local on one PE, should be local on all
    indicesAreGlobal_ = (allIndices[1]==1);  // If indices are global on one PE should be local on all
  }



  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::sortAllIndices()
  {
    TEUCHOS_TEST_FOR_EXCEPT(isGloballyIndexed()==true);   // this should be called only after makeIndicesLocal()
    if (isSorted() == false) {
      for (size_t row=0; row < getNodeNumRows(); ++row) {
        RowInfo rowInfo = getRowInfo(row);
        sortRowIndices(rowInfo);
      }
    }
    // we just sorted every row
    setSorted(true);
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::makeColMap()
  {
    std::string tfecfFuncName("makeColMap()");
    typedef OrdinalTraits<GlobalOrdinal> GOT;
    //
    if (hasColMap()) return;
    const size_t nlrs = getNodeNumRows();
    //
    computeIndexState();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(isLocallyIndexed() == true, std::runtime_error, ": indices must still be global when making the column map.");
    // ultimate goal: list of indices for column map
    Array<GlobalOrdinal> myColumns;
    // if isGloballyIndexed() == false and isLocallyIndexed() == false, then we have nothing to do
    if (isGloballyIndexed() == true) {
      // Construct two lists of columns for which we have a non-zero
      // Local GIDs: Column GIDs which are locally present on the Domain map
      // Remote GIDs: Column GIDs which are not locally present present on the Domain map
      //
      // instead of list of local GIDs, we will use a set<LocalOrdinal> LocalGID.
      //
      const LocalOrdinal LINV = OrdinalTraits<LocalOrdinal>::invalid();
      size_t numLocalColGIDs = 0, numRemoteColGIDs = 0;
      //
      // initial: partitioning into local and remote
      Array<char> GIDisLocal(domainMap_->getNodeNumElements(),0);
      std::set<GlobalOrdinal> RemoteGIDSet;
      for (size_t r=0; r < nlrs; ++r) {
        RowInfo rowinfo = getRowInfo(r);
        if (rowinfo.numEntries > 0) {
          ArrayView<const GlobalOrdinal> rowgids = getGlobalView(rowinfo);
          rowgids = rowgids(0,rowinfo.numEntries);
          for (typename ArrayView<const GlobalOrdinal>::iterator cind = rowgids.begin(); cind != rowgids.end(); ++cind)
          {
            GlobalOrdinal gid = (*cind);
            LocalOrdinal lid = domainMap_->getLocalElement(gid);
            if (lid != LINV) {
              char alreadyFound = GIDisLocal[lid];
              if (alreadyFound == 0) {
                GIDisLocal[lid] = 1;
                ++numLocalColGIDs;
              }
            }
            else {
              std::pair<typename std::set<GlobalOrdinal>::iterator, bool> ip;
              ip = RemoteGIDSet.insert(gid);
              if (ip.second == true) { // gid did not exist in the set and was actually inserted
                ++numRemoteColGIDs;
              }
            }
          }
        }
      }

      // Possible short-circuit for serial scenario
      // If the all domain GIDs are present as column indices, then set ColMap=DomainMap
      // By construction, LocalGIDs \subset DomainGIDs
      // If we have
      //   * Number of remote GIDs is 0, so that ColGIDs == LocalGIDs,
      // and
      //   * Number of local GIDs is number of domain GIDs
      // then
      //   * LocalGIDs \subset DomainGIDs && size(LocalGIDs) == size(DomainGIDs) => DomainGIDs == LocalGIDs == ColGIDs
      // on this node.
      // We will concern ourselves only with the special case of a serial DomainMap, obviating the need for a communication.
      // If
      //   * DomainMap has a serial comm
      // then we can set Column map as Domain map and return. Benefit: this graph won't need an Import object later.
      //
      // Note, for a serial domain map, there can be no RemoteGIDs, because there are no remote nodes.
      // Likely explanations for this are:
      //  * user submitted erroneous column indices
      //  * user submitted erroneous domain map
      if (Teuchos::size(*domainMap_->getComm()) == 1)
      {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(numRemoteColGIDs != 0, std::runtime_error,
            ": Some column IDs are not in the domain map." << std::endl
            << "Either these column IDs are invalid or the domain map is invalid." << std::endl
            << "Remember, for a rectangular matrix, the domain map must be passed to fillComplete().");
        if (numLocalColGIDs == domainMap_->getNodeNumElements()) {
          colMap_ = domainMap_;
          checkInternalState();
          return;
        }
      }

      // Now, populate myColumns() with a list of all column GIDs.
      // Put local GIDs at the front: they correspond to "same" and "permuted" entries between the column map and the domain map
      // Put remote GIDs at the back
      myColumns.resize(numLocalColGIDs + numRemoteColGIDs);
      // get pointers into myColumns for each part
      ArrayView<GlobalOrdinal> LocalColGIDs, RemoteColGIDs;
      LocalColGIDs  = myColumns(0,numLocalColGIDs);
      RemoteColGIDs = myColumns(numLocalColGIDs,numRemoteColGIDs);

      // Copy the Remote GIDs into myColumns
      std::copy(RemoteGIDSet.begin(), RemoteGIDSet.end(), RemoteColGIDs.begin());
      // We will make a list of corresponding node IDs
      Array<int> RemoteImageIDs(numRemoteColGIDs);
      // Lookup the Remote Nodes IDs in the Domain map
      {
        LookupStatus stat = domainMap_->getRemoteIndexList(RemoteColGIDs, RemoteImageIDs());
        // This error check is crucial: this tells us that one of the remote indices was not present in the domain map.
        // This means that the Import object cannot be constructed, because of incongruity between the column map and domain map.
        //   * The user has made a mistake in the column indices
        //   * The user has made a mistake w.r.t. the domain map
        // Same error message as above for serial case.
        char missingID_lcl = (stat == IDNotPresent ? 1 : 0);
        char missingID_gbl;
        Teuchos::reduceAll<int,char>(*getComm(),Teuchos::REDUCE_MAX,missingID_lcl,
          outArg(missingID_gbl));
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(missingID_gbl == 1, std::runtime_error,
            ": Some column IDs are not in the domain map." << std::endl
            << "Either these column IDs are invalid or the domain map is invalid." << std::endl
            << "Common cause: for a rectangular matrix, the domain map must be passed to fillComplete().");
      }
      // Sort External column indices so that all columns coming from a given remote processor are contiguous
      // This obviates the need for the Distributor associated with the Import from having to reorder data.
      sort2(RemoteImageIDs.begin(), RemoteImageIDs.end(), RemoteColGIDs.begin());

      // Copy the Local GIDs into myColumns. Two cases:
      // (1) If the number of Local column GIDs is the same as the number of Local domain GIDs, we
      //     can simply read the domain GIDs into the front part of ColIndices (see logic above from the serial short circuit case)
      // (2) We step through the GIDs of the DomainMap, checking to see if each domain GID is a column GID.
      //     we want to do this to maintain a consistent ordering of GIDs between the columns and the domain.
      ArrayView<const GlobalOrdinal> mge = domainMap_->getNodeElementList();
      if (numLocalColGIDs == domainMap_->getNodeNumElements()) {
        std::copy(mge.begin(), mge.end(), LocalColGIDs.begin());
      }
      else {
        size_t numlocalagain = 0;
        for (size_t i=0; i < domainMap_->getNodeNumElements(); ++i) {
          if (GIDisLocal[i]) {
            LocalColGIDs[numlocalagain++] = mge[i];
          }
        }
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(numlocalagain != numLocalColGIDs, std::logic_error, ": Internal logic error. Please contact Tpetra team.");
      }
    }
    colMap_ = rcp(new Map<LocalOrdinal,GlobalOrdinal,Node>(GOT::invalid(), myColumns, domainMap_->getIndexBase(), domainMap_->getComm(), domainMap_->getNode()) );
    checkInternalState();
    return;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::mergeAllIndices()
  {
    TEUCHOS_TEST_FOR_EXCEPT( isGloballyIndexed() ); // call only after makeIndicesLocal()
    TEUCHOS_TEST_FOR_EXCEPT( ! isSorted() ); // call only after sortIndices()
    if (! isMerged ()) {
      for (size_t row=0; row < getNodeNumRows(); ++row) {
        RowInfo rowInfo = getRowInfo(row);
        mergeRowIndices(rowInfo);
      }
      // we just merged every row
      setMerged(true);
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::makeImportExport()
  {
    typedef Import<LocalOrdinal,GlobalOrdinal,Node> import_type;
    typedef Export<LocalOrdinal,GlobalOrdinal,Node> export_type;

    TEUCHOS_TEST_FOR_EXCEPTION(! hasColMap (), std::logic_error, "Tpetra::"
      "CrsGraph: It's not allowed to call makeImportExport() unless the graph "
      "has a column Map.");
    RCP<ParameterList> params = this->getNonconstParameterList (); // could be null

    // Create the Import instance if necessary.
    if (domainMap_ != colMap_ && (! domainMap_->isSameAs (*colMap_))) {
      if (params.is_null () || ! params->isSublist ("Import")) {
        importer_ = rcp (new import_type (domainMap_, colMap_));
      }
      else {
        RCP<ParameterList> importSublist = sublist (params, "Import", true);
        importer_ = rcp (new import_type (domainMap_, colMap_, importSublist));
      }
    }
    else {
      importer_ = null;
    }

    // Create the Export instance if necessary.
    if (rangeMap_ != rowMap_ && (!rangeMap_->isSameAs(*rowMap_))) {
      if (params.is_null () || ! params->isSublist ("Export")) {
        exporter_ = rcp (new export_type (rowMap_, rangeMap_));
      }
      else {
        RCP<ParameterList> exportSublist = sublist (params, "Export", true);
        exporter_ = rcp (new export_type (rowMap_, rangeMap_, exportSublist));
      }
    }
    else {
      exporter_ = null;
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::description() const
  {
    std::ostringstream oss;
    oss << DistObject<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Node>::description();
    if (isFillComplete()) {
      oss << "{status = fill complete"
          << ", global rows = " << getGlobalNumRows()
          << ", global cols = " << getGlobalNumCols()
          << ", global num entries = " << getGlobalNumEntries()
          << "}";
    }
    else {
      oss << "{status = fill not complete"
          << ", global rows = " << getGlobalNumRows()
          << "}";
    }
    return oss.str();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  describe (Teuchos::FancyOStream &out,
            const Teuchos::EVerbosityLevel verbLevel) const
  {
    using std::endl;
    using std::setw;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;
    Teuchos::EVerbosityLevel vl = verbLevel;
    if (vl == VERB_DEFAULT) vl = VERB_LOW;
    RCP<const Comm<int> > comm = this->getComm();
    const int myImageID = comm->getRank(),
              numImages = comm->getSize();
    size_t width = 1;
    for (size_t dec=10; dec<getGlobalNumRows(); dec *= 10) {
      ++width;
    }
    width = std::max<size_t>(width,Teuchos::as<size_t>(11)) + 2;
    Teuchos::OSTab tab(out);
    //    none: print nothing
    //     low: print O(1) info from node 0
    //  medium: print O(P) info, num entries per node
    //    high: print O(N) info, num entries per row
    // extreme: print O(NNZ) info: print graph indices
    //
    // for medium and higher, print constituent objects at specified verbLevel
    if (vl != VERB_NONE) {
      if (myImageID == 0) out << this->description() << std::endl;
      // O(1) globals, minus what was already printed by description()
      if (isFillComplete() && myImageID == 0) {
        out << "Global number of diagonals = " << globalNumDiags_ << std::endl;
        out << "Global max number of row entries = " << globalMaxNumRowEntries_ << std::endl;
      }
      // constituent objects
      if (vl == VERB_MEDIUM || vl == VERB_HIGH || vl == VERB_EXTREME) {
        if (myImageID == 0) out << "\nRow map: " << std::endl;
        rowMap_->describe(out,vl);
        if (colMap_ != null) {
          if (myImageID == 0) out << "\nColumn map: " << std::endl;
          colMap_->describe(out,vl);
        }
        if (domainMap_ != null) {
          if (myImageID == 0) out << "\nDomain map: " << std::endl;
          domainMap_->describe(out,vl);
        }
        if (rangeMap_ != null) {
          if (myImageID == 0) out << "\nRange map: " << std::endl;
          rangeMap_->describe(out,vl);
        }
      }
      // O(P) data
      if (vl == VERB_MEDIUM || vl == VERB_HIGH || vl == VERB_EXTREME) {
        for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
          if (myImageID == imageCtr) {
            out << "Node ID = " << imageCtr << std::endl
                << "Node number of entries = " << nodeNumEntries_ << std::endl
                << "Node number of diagonals = " << nodeNumDiags_ << std::endl
                << "Node max number of entries = " << nodeMaxNumRowEntries_ << std::endl;
            if (indicesAreAllocated()) {
              out << "Node number of allocated entries = " << nodeNumAllocated_ << std::endl;
            }
            else {
              out << "Indices are not allocated." << std::endl;
            }
          }
          comm->barrier();
          comm->barrier();
          comm->barrier();
        }
      }
      // O(N) and O(NNZ) data
      if (vl == VERB_HIGH || vl == VERB_EXTREME) {
        for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
          if (myImageID == imageCtr) {
            out << std::setw(width) << "Node ID"
                << std::setw(width) << "Global Row"
                << std::setw(width) << "Num Entries";
            if (vl == VERB_EXTREME) {
              out << "Entries";
            }
            out << std::endl;
            for (size_t r=0; r < getNodeNumRows(); ++r) {
              RowInfo rowinfo = getRowInfo(r);
              GlobalOrdinal gid = rowMap_->getGlobalElement(r);
              out << std::setw(width) << myImageID
                  << std::setw(width) << gid
                  << std::setw(width) << rowinfo.numEntries;
              if (vl == VERB_EXTREME) {
                if (isGloballyIndexed()) {
                  ArrayView<const GlobalOrdinal> rowview = getGlobalView(rowinfo);
                  for (size_t j=0; j < rowinfo.numEntries; ++j) out << rowview[j] << " ";
                }
                else if (isLocallyIndexed()) {
                  ArrayView<const LocalOrdinal> rowview = getLocalView(rowinfo);
                  for (size_t j=0; j < rowinfo.numEntries; ++j) out << colMap_->getGlobalElement(rowview[j]) << " ";
                }
              }
              out << std::endl;
            }
          }
          comm->barrier();
          comm->barrier();
          comm->barrier();
        }
      }
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::checkSizes(const DistObject<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Node>& source)
  {
    // It's not clear what kind of compatibility checks on sizes can be performed here.
    // Epetra_CrsGraph doesn't check any sizes for compatibility.
    return true;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::copyAndPermute(
                          const DistObject<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Node> & source,
                          size_t numSameIDs,
                          const ArrayView<const LocalOrdinal> &permuteToLIDs,
                          const ArrayView<const LocalOrdinal> &permuteFromLIDs)
  {
    const CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>& src_graph = dynamic_cast<const CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>&>(source);
    std::string tfecfFuncName("copyAndPermute()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(permuteToLIDs.size() != permuteFromLIDs.size(), std::runtime_error, ": permuteToLIDs and permuteFromLIDs must have the same size.");
    bool src_filled = src_graph.isFillComplete();

    Array<GlobalOrdinal> row_copy;
    LocalOrdinal myid = 0;
    for (size_t i=0; i<numSameIDs; ++i, ++myid) {
      GlobalOrdinal gid = src_graph.getMap()->getGlobalElement(myid);
      if (src_filled) {
        size_t row_length = src_graph.getNumEntriesInGlobalRow(gid);
        row_copy.resize(row_length);
        size_t check_row_length = 0;
        src_graph.getGlobalRowCopy(gid, row_copy(), check_row_length);
        insertGlobalIndices( gid, row_copy() );
      }
      else {
        ArrayView<const GlobalOrdinal> row;
        src_graph.getGlobalRowView( gid,row );
        insertGlobalIndices( gid, row );
      }
    }

    for (LocalOrdinal i=0; i<permuteToLIDs.size(); ++i) {
      GlobalOrdinal mygid = this->getMap()->getGlobalElement(permuteToLIDs[i]);
      GlobalOrdinal srcgid= src_graph.getMap()->getGlobalElement(permuteFromLIDs[i]);
      if (src_filled) {
        size_t row_length = src_graph.getNumEntriesInGlobalRow(srcgid);
        row_copy.resize(row_length);
        size_t check_row_length = 0;
        src_graph.getGlobalRowCopy(srcgid, row_copy(), check_row_length);
        insertGlobalIndices( mygid, row_copy() );
      }
      else {
        ArrayView<const GlobalOrdinal> row;
        src_graph.getGlobalRowView( srcgid, row );
        insertGlobalIndices( mygid, row );
      }
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  packAndPrepare (const DistObject<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Node> & source,
                  const ArrayView<const LocalOrdinal> &exportLIDs,
                  Array<GlobalOrdinal> &exports,
                  const ArrayView<size_t> & numPacketsPerLID,
                  size_t& constantNumPackets,
                  Distributor &distor)
  {
    std::string tfecfFuncName("packAndPrepare()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(exportLIDs.size() != numPacketsPerLID.size(), std::runtime_error,
        ": exportLIDs and numPacketsPerLID must have the same size.");
    const CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>& src_graph = dynamic_cast<const CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>&>(source);
    // We don't check whether src_graph has had fillComplete called, because it doesn't matter whether the
    // *source* graph has been fillComplete'd. The target graph can not be fillComplete'd yet.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(this->isFillComplete() == true, std::runtime_error,
         ": import/export operations are not allowed on destination CrsGraph after fillComplete has been called.");
    constantNumPackets = 0;
    // first set the contents of numPacketsPerLID, and accumulate a total-num-packets:
    size_t totalNumPackets = 0;
    Array<GlobalOrdinal> row;
    for (LocalOrdinal i=0; i<exportLIDs.size(); ++i) {
      GlobalOrdinal GID = src_graph.getMap()->getGlobalElement(exportLIDs[i]);
      size_t row_length = src_graph.getNumEntriesInGlobalRow(GID);
      numPacketsPerLID[i] = row_length;
      totalNumPackets += row_length;
    }

    exports.resize(totalNumPackets);

    // now loop again and pack rows of indices into exports:
    size_t exportsOffset = 0;
    for (LocalOrdinal i=0; i<exportLIDs.size(); ++i) {
      GlobalOrdinal GID = src_graph.getMap()->getGlobalElement(exportLIDs[i]);
      size_t row_length = src_graph.getNumEntriesInGlobalRow(GID);
      row.resize(row_length);
      size_t check_row_length = 0;
      src_graph.getGlobalRowCopy(GID, row(), check_row_length);
      typename Array<GlobalOrdinal>::const_iterator
        row_iter = row.begin(), row_end = row.end();
      size_t j = 0;
      for (; row_iter != row_end; ++row_iter, ++j) {
        exports[exportsOffset+j] = *row_iter;
      }
      exportsOffset += row.size();
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::unpackAndCombine(
                            const ArrayView<const LocalOrdinal> &importLIDs,
                            const ArrayView<const GlobalOrdinal> &imports,
                            const ArrayView<size_t> &numPacketsPerLID,
                            size_t constantNumPackets,
                            Distributor & /* distor */,
                            CombineMode /* CM */)
  {
    // FIXME (mfh 02 Apr 2012) REPLACE combine mode has a perfectly
    // reasonable meaning, whether or not the matrix is fill complete.
    // It's just more work to implement.

    // We are not checking the value of the CombineMode input-argument.
    // For CrsGraph, we only support import/export operations if fillComplete has not yet been called.
    // Any incoming column-indices are inserted into the target graph. In this context, CombineMode values
    // of ADD vs INSERT are equivalent. What is the meaning of REPLACE for CrsGraph? If a duplicate column-index
    // is inserted, it will be compressed out when fillComplete is called.
    //
    // Note: I think REPLACE means that an existing row is replaced by the imported row, i.e., the existing indices are cleared. CGB, 6/17/2010

    std::string tfecfFuncName("unpackAndCombine()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(importLIDs.size() != numPacketsPerLID.size(), std::runtime_error,
        ": importLIDs and numPacketsPerLID must have the same size.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(this->isFillComplete() == true, std::runtime_error,
        ": import/export operations are not allowed on destination CrsGraph after fillComplete has been called.");
    size_t importsOffset = 0;
    typename ArrayView<const LocalOrdinal>::iterator
      impLIDiter = importLIDs.begin(), impLIDend = importLIDs.end();
    size_t i = 0;
    for (; impLIDiter != impLIDend; ++impLIDiter, ++i) {
      LocalOrdinal row_length = numPacketsPerLID[i];
      const ArrayView<const GlobalOrdinal> row(&imports[importsOffset], row_length);
      insertGlobalIndices(this->getMap()->getGlobalElement(*impLIDiter), row);
      importsOffset += row_length;
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                         Deprecated methods                              //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////////////////////
  // DEPRECATED
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::fillComplete(OptimizeOption os)
  {
    fillComplete(rowMap_,rowMap_,os);
  }


  /////////////////////////////////////////////////////////////////////////////
  // DEPRECATED
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::fillComplete(
                                    const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap,
                                    const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap,
                                    OptimizeOption os)
  {
    RCP<ParameterList> params = parameterList();
    if (os == DoOptimizeStorage) params->set<bool>("Optimize Storage",true);
    else                         params->set<bool>("Optimize Storage",false);
    fillComplete(domainMap,rangeMap,params);
  }

  /////////////////////////////////////////////////////////////////////////////
  // DEPRECATED
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  ArrayRCP<const LocalOrdinal> CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getLocalRowView(LocalOrdinal localRow) const
  {
    std::string tfecfFuncName("getLocalRowView(localRow)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(isGloballyIndexed() == true, std::runtime_error,
        ": local indices do not exist.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(rowMap_->isNodeLocalElement(localRow) == false, std::runtime_error,
        ": localRow (== " << localRow << ") is not valid on this node.");
    ArrayRCP<const LocalOrdinal> ret;
    RowInfo sizeInfo = getRowInfo(localRow);
    if (getProfileType() == StaticProfile && sizeInfo.numEntries > 0) {
      ret = lclInds1D_.persistingView(sizeInfo.offset1D,sizeInfo.numEntries);
    }
    else if (getProfileType() == DynamicProfile && sizeInfo.numEntries > 0) {
      ret = lclInds2D_[localRow].persistingView(0,sizeInfo.numEntries);
    }
    return ret;
  }


  /////////////////////////////////////////////////////////////////////////////
  // DEPRECATED
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  ArrayRCP<const GlobalOrdinal> CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getGlobalRowView(GlobalOrdinal globalRow) const
  {
    std::string tfecfFuncName("getGlobalRowView(globalRow)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(isLocallyIndexed() == true, std::runtime_error,
        ": global indices do not exist.");
    const LocalOrdinal lrow = rowMap_->getLocalElement(globalRow);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(lrow == OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error,
        ": globalRow (== " << globalRow << ") does not belong to this node.");
    ArrayRCP<const GlobalOrdinal> ret;
    RowInfo sizeInfo = getRowInfo(lrow);
    if (getProfileType() == StaticProfile && sizeInfo.numEntries > 0) {
      ret = gblInds1D_.persistingView(sizeInfo.offset1D,sizeInfo.numEntries);
    }
    else if (getProfileType() == DynamicProfile && sizeInfo.numEntries > 0) {
      ret = gblInds2D_[lrow].persistingView(0,sizeInfo.numEntries);
    }
    return ret;
  }


  /////////////////////////////////////////////////////////////////////////////
  // DEPRECATED
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::optimizeStorage ()
  {
    // provided only for backwards compatibility
    // previous semantics required that fillComplete() had been called.
    std::string tfecfFuncName("optimizeStorage()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(! isFillComplete(), std::runtime_error,
      ": requires that fillComplete() has already been called.");
    if (! isStorageOptimized ()) {
      resumeFill ();
      fillComplete (DoOptimizeStorage);
    }
  }

} // namespace Tpetra


//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_CRSGRAPH_INSTANT(LO,GO,NODE) \
  \
  template class CrsGraph< LO , GO , NODE >; \


#endif // TPETRA_CRSGRAPH_DEF_HPP
