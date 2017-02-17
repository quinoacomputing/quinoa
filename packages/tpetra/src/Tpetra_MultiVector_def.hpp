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

#ifndef TPETRA_MULTIVECTOR_DEF_HPP
#define TPETRA_MULTIVECTOR_DEF_HPP

#include <Kokkos_NodeTrace.hpp>

#include <Teuchos_Assert.hpp>
#include <Teuchos_as.hpp>

#include "Tpetra_Vector.hpp"

#ifdef DOXYGEN_USE_ONLY
  #include "Tpetra_MultiVector_decl.hpp"
#endif

namespace Tpetra {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
               size_t NumVectors,
               bool zeroOut) : /* default is true */
    DistObject<Scalar,LocalOrdinal,GlobalOrdinal,Node> (map),
    lclMV_ (map->getNode ())
  {
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;

    TEUCHOS_TEST_FOR_EXCEPTION(NumVectors < 1, std::invalid_argument,
      "Tpetra::MultiVector::MultiVector(): NumVectors must be strictly positive.");
    const size_t myLen = getLocalLength();
    if (myLen > 0) {
      RCP<Node> node = map->getNode();
      ArrayRCP<Scalar> data = node->template allocBuffer<Scalar>(myLen*NumVectors);
      MVT::initializeValues(lclMV_,myLen,NumVectors,data,myLen);
      if (zeroOut) {
        MVT::Init(lclMV_, Teuchos::ScalarTraits<Scalar>::zero());
      }
    }
    else {
      MVT::initializeValues(lclMV_,0,NumVectors,Teuchos::null,0);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  MultiVector (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& source) :
    DistObject<Scalar,LocalOrdinal,GlobalOrdinal,Node> (source),
    lclMV_ (MVT::getNode (source.lclMV_))
  {
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;

    // copy data from the source MultiVector into this multivector
    RCP<Node> node = MVT::getNode(source.lclMV_);
    const LocalOrdinal myLen = getLocalLength();
    const size_t numVecs = source.getNumVectors();
    if (myLen > 0) {
      // allocate data
      ArrayRCP<Scalar> data = node->template allocBuffer<Scalar>(myLen*numVecs);
      MVT::initializeValues(lclMV_,myLen,numVecs,data,myLen);
      // copy data
      {
        ArrayRCP<Scalar> dstdata = data;
        ArrayRCP<const Scalar> srcdata = MVT::getValues(source.lclMV_);
        for (size_t j = 0; j < numVecs; ++j) {
          ArrayRCP<const Scalar> srcj = source.getSubArrayRCP(srcdata,j);
          KOKKOS_NODE_TRACE("MultiVector::MultiVector(MV)")
          node->template copyBuffers<Scalar>(myLen,srcj,dstdata);
          dstdata += myLen;
        }
      }
    }
    else {
      MVT::initializeValues(lclMV_,0,numVecs,Teuchos::null,0);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
               const Teuchos::ArrayRCP<Scalar>& view,
               size_t LDA,
               size_t NumVectors,
               EPrivateHostViewConstructor /* dummy */) :
    DistObject<Scalar,LocalOrdinal,GlobalOrdinal,Node> (map),
    lclMV_ (map->getNode ())
  {
    using Teuchos::as;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;

    const std::string tfecfFuncName("MultiVector(view,LDA,NumVector)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(NumVectors < 1, std::invalid_argument,
      ": NumVectors must be strictly positive, but you specified NumVectors = "
      << NumVectors << ".");
    const size_t myLen = getLocalLength();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(LDA < myLen, std::invalid_argument,
      ": LDA must be large enough to accomodate the local entries.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      as<size_t> (view.size ()) < LDA*(NumVectors - 1) + myLen,
      std::invalid_argument,
      ": view does not contain enough data to specify the entries in this.");
    MVT::initializeValues (lclMV_, myLen, NumVectors, view, myLen);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map,
               Teuchos::ArrayRCP<Scalar> data,
               size_t LDA,
               size_t NumVectors,
               EPrivateComputeViewConstructor /* dummy */) :
    DistObject<Scalar,LocalOrdinal,GlobalOrdinal,Node> (map),
    lclMV_ (map->getNode ())
  {
    using Teuchos::as;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;

    const std::string tfecfFuncName("MultiVector(data,LDA,NumVector)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(NumVectors < 1, std::invalid_argument,
      ": NumVectors must be strictly positive.");
    const size_t myLen = getLocalLength();
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      as<size_t> (data.size ()) < LDA*(NumVectors - 1) + myLen,
      std::runtime_error,
      ": data does not contain enough data to specify the entries in this.");
#endif
    MVT::initializeValues(lclMV_,myLen,NumVectors,data,LDA);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map,
               Teuchos::ArrayRCP<Scalar> data,
               size_t LDA,
               Teuchos::ArrayView<const size_t> WhichVectors,
               EPrivateComputeViewConstructor /* dummy */) :
    DistObject<Scalar,LocalOrdinal,GlobalOrdinal,Node> (map),
    lclMV_ (map->getNode ()),
    whichVectors_ (WhichVectors)
  {
    using Teuchos::as;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;

    const std::string tfecfFuncName("MultiVector(data,LDA,WhichVectors)");
    const size_t myLen = getLocalLength();
    size_t maxVector = *std::max_element(WhichVectors.begin(), WhichVectors.end());
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(LDA < myLen, std::runtime_error,
      ": LDA must be large enough to accomodate the local entries.");
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      as<size_t> (data.size ()) < LDA * maxVector + myLen, std::runtime_error,
      ": data does not contain enough data to specify the entries in this.");
#endif
    if (WhichVectors.size () == 1) {
      // shift data so that desired vector is vector 0
      maxVector = 0;
      data += LDA * WhichVectors[0];
      // kill whichVectors_; we are constant stride
      whichVectors_.clear ();
    }
    if (myLen > 0) {
      MVT::initializeValues(lclMV_,myLen,maxVector+1,data,LDA);
    }
    else {
      MVT::initializeValues(lclMV_,0,WhichVectors.size(),Teuchos::null,0);
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
               const Teuchos::ArrayView<const Scalar>& A,
               size_t LDA,
               size_t NumVectors) :
    DistObject<Scalar,LocalOrdinal,GlobalOrdinal,Node> (map),
    lclMV_ (map->getNode ())
  {
    using Teuchos::as;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::RCP;

    const std::string tfecfFuncName("MultiVector(A,LDA)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(NumVectors < 1, std::invalid_argument,
      ": NumVectors must be strictly positive.");
    const size_t myLen = getLocalLength();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(LDA < myLen, std::runtime_error,
      ": LDA must be large enough to accomodate the local entries.");
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      as<size_t> (A.size ()) < LDA*(NumVectors - 1) + myLen,
      std::runtime_error,
      ": A does not contain enough data to specify the entries in this.");
#endif
    if (myLen > 0) {
      RCP<Node> node = MVT::getNode(lclMV_);
      ArrayRCP<Scalar> mydata = node->template allocBuffer<Scalar>(myLen*NumVectors);
      MVT::initializeValues(lclMV_,myLen,NumVectors,mydata,myLen);
      KOKKOS_NODE_TRACE("MultiVector::MultiVector(1D)")
      ArrayRCP<Scalar> myview = node->template viewBufferNonConst<Scalar>(Kokkos::WriteOnly,myLen*NumVectors,mydata);
      typename ArrayView<const Scalar>::iterator srcit = A.begin();
      for (size_t j = 0; j < NumVectors; ++j) {
        std::copy(srcit,srcit+myLen,myview);
        srcit += LDA;
        myview += myLen;
      }
      mydata = Teuchos::null;
      myview = Teuchos::null;
    }
    else {
      MVT::initializeValues(lclMV_,0,NumVectors,Teuchos::null,0);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
               const Teuchos::ArrayView<const ArrayView<const Scalar> >& ArrayOfPtrs,
               size_t NumVectors) :
    DistObject<Scalar,LocalOrdinal,GlobalOrdinal,Node> (map),
    lclMV_ (map->getNode ())
  {
    using Teuchos::as;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::RCP;

    const std::string tfecfFuncName("MultiVector(ArrayOfPtrs)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      NumVectors < 1 || NumVectors != as<size_t>(ArrayOfPtrs.size()),
      std::runtime_error,
      ": ArrayOfPtrs.size() must be strictly positive and as large as ArrayOfPtrs.");
    const size_t myLen = getLocalLength();
    if (myLen > 0) {
      RCP<Node> node = MVT::getNode(lclMV_);
      ArrayRCP<Scalar> mydata = node->template allocBuffer<Scalar>(myLen*NumVectors);
      MVT::initializeValues(lclMV_,myLen,NumVectors,mydata,myLen);
      KOKKOS_NODE_TRACE("MultiVector::MultiVector(2D)")
      ArrayRCP<Scalar> myview = node->template viewBufferNonConst<Scalar>(Kokkos::WriteOnly,myLen*NumVectors,mydata);
      for (size_t j = 0; j < NumVectors; ++j) {
#ifdef HAVE_TPETRA_DEBUG
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          as<size_t>(ArrayOfPtrs[j].size()) != getLocalLength(),
          std::runtime_error,
          ": ArrayOfPtrs[" << j << "].size() (== " << ArrayOfPtrs[j].size()
          << ") is not equal to getLocalLength() (== " << getLocalLength());
#endif
        typename ArrayView<const Scalar>::iterator src = ArrayOfPtrs[j].begin();
        std::copy(src,src+myLen,myview);
        myview += myLen;
      }
      myview = Teuchos::null;
      mydata = Teuchos::null;
    }
    else {
      MVT::initializeValues(lclMV_,0,NumVectors,Teuchos::null,0);
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::~MultiVector() {
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::isConstantStride() const {
    return whichVectors_.empty();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getLocalLength() const {
    return this->getMap()->getNodeNumElements();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  global_size_t MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getGlobalLength() const {
    return this->getMap()->getGlobalNumElements();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getStride() const {
    if (isConstantStride()) {
      return MVT::getStride(lclMV_);
    }
    return 0;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  checkSizes(const DistObject<Scalar,LocalOrdinal,GlobalOrdinal,Node> &sourceObj)
  {
    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;
    using Teuchos::rcpFromRef;

    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
    // rcp_dynamic_cast gives us superior cast failure output to dynamic_cast.
    RCP<const MV> A = rcp_dynamic_cast<const MV> (rcpFromRef (sourceObj), true);

    // This method is called in DistObject::doTransfer().  By that
    // point, we've already constructed an Import or Export object
    // using the two multivectors' Maps, which means that (hopefully)
    // we've already checked other attributes of the multivectors.
    // Thus, all we need to do here is check the number of columns.
    return (A->getNumVectors() == this->getNumVectors());
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::copyAndPermute(
                          const DistObject<Scalar,LocalOrdinal,GlobalOrdinal,Node> & sourceObj,
                          size_t numSameIDs,
                          const Teuchos::ArrayView<const LocalOrdinal> &permuteToLIDs,
                          const Teuchos::ArrayView<const LocalOrdinal> &permuteFromLIDs)
  {
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::RCP;

    const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &sourceMV = dynamic_cast<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &>(sourceObj);
    typename ArrayView<const LocalOrdinal>::iterator pTo, pFrom;
    // any other error will be caught by Teuchos
    TEUCHOS_TEST_FOR_EXCEPTION(permuteToLIDs.size() != permuteFromLIDs.size(), std::runtime_error,
      "Tpetra::MultiVector::copyAndPermute(): permuteToLIDs and permuteFromLIDs must have the same size.");
    RCP<Node> node = MVT::getNode(lclMV_);
    const size_t numCols = getNumVectors();

    // Copy the vectors one at a time.  This works whether or not the
    // multivectors have constant stride.
    for (size_t j = 0; j < numCols; ++j) {
      // Get host views of the current source and destination column.
      //
      // FIXME (mfh 10 Mar 2012) Copying should really be done on the
      // device.  There's no need to bring everything back to the
      // host, copy there, and then copy back.
      ArrayRCP<const Scalar> srcptr = sourceMV.getSubArrayRCP(sourceMV.cview_,j);
      ArrayRCP<      Scalar> dstptr =          getSubArrayRCP(ncview_,j);

      // The first numImportIDs IDs are the same between source and
      // target, so we can just copy the data.  (This favors faster
      // contiguous access whenever we can do it.)
      std::copy(srcptr,srcptr+numSameIDs,dstptr);

      // For the remaining GIDs, execute the permutations.  This may
      // involve noncontiguous access of both source and destination
      // vectors, depending on the LID lists.
      //
      // FIXME (mfh 09 Apr 2012) Permutation should be done on the
      // device, not on the host.
      //
      // FIXME (mfh 20 June 2012) For an Export with duplicate GIDs on
      // the same process, this merges their values by replacement of
      // the last encountered GID, not by the specified merge rule
      // (such as ADD).
      for (pTo = permuteToLIDs.begin(), pFrom = permuteFromLIDs.begin();
           pTo != permuteToLIDs.end(); ++pTo, ++pFrom) {
        dstptr[*pTo] = srcptr[*pFrom];
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::packAndPrepare(
          const DistObject<Scalar,LocalOrdinal,GlobalOrdinal,Node> & sourceObj,
          const ArrayView<const LocalOrdinal> &exportLIDs,
          Array<Scalar> &exports,
          const ArrayView<size_t> &numExportPacketsPerLID,
          size_t& constantNumPackets,
          Distributor & /* distor */ )
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::as;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

    const MV& sourceMV = dynamic_cast<const MV&> (sourceObj);

    /* The layout in the export for MultiVectors is as follows:
       exports = { all of the data from row exportLIDs.front() ;
                   ....
                   all of the data from row exportLIDs.back() }
      This doesn't have the best locality, but is necessary because
      the data for a Packet (all data associated with an LID) is
      required to be contiguous. */
    TEUCHOS_TEST_FOR_EXCEPTION(as<int> (numExportPacketsPerLID.size ()) != exportLIDs.size (),
      std::runtime_error, "Tpetra::MultiVector::packAndPrepare(): size of num"
      "ExportPacketsPerLID buffer should be the same as exportLIDs.  numExport"
      "PacketsPerLID.size() = " << numExportPacketsPerLID.size() << ", but "
      "exportLIDs.size() = " << exportLIDs.size() << ".");
    const KMV& srcData = sourceMV.lclMV_;
    const size_t numCols = sourceMV.getNumVectors ();
    const size_t stride = MVT::getStride (srcData);
    constantNumPackets = numCols;
    exports.resize (numCols * exportLIDs.size ());
    typename ArrayView<const LocalOrdinal>::iterator idptr;
    typename Array<Scalar>::iterator expptr;
    expptr = exports.begin();

    if (sourceMV.isConstantStride ()) {
      size_t i = 0;
      for (idptr = exportLIDs.begin(); idptr != exportLIDs.end(); ++idptr, ++i) {
        for (size_t j = 0; j < numCols; ++j) {
          *expptr++ = sourceMV.cview_[j*stride + (*idptr)];
        }
        //we shouldn't need to set numExportPacketsPerLID[i] since we have set
        //constantNumPackets to a nonzero value. But we'll set it anyway, since
        //I'm not sure if the API will remain the way it is.
        numExportPacketsPerLID[i] = numCols;
      }
    }
    else {
      size_t i = 0;
      for (idptr = exportLIDs.begin(); idptr != exportLIDs.end(); ++idptr, ++i) {
        for (size_t j = 0; j < numCols; ++j) {
          *expptr++ = sourceMV.cview_[sourceMV.whichVectors_[j]*stride + (*idptr)];
        }
        numExportPacketsPerLID[i] = numCols;
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  unpackAndCombine (const ArrayView<const LocalOrdinal> &importLIDs,
                    const ArrayView<const Scalar> &imports,
                    const ArrayView<size_t> &numPacketsPerLID,
                    size_t constantNumPackets,
                    Distributor & /* distor */,
                    CombineMode CM)
  {
    using Teuchos::ArrayView;
    using Teuchos::as;
    typedef Teuchos::ScalarTraits<Scalar> SCT;

    const std::string tfecfFuncName("unpackAndCombine()");
    /* The layout in the export for MultiVectors is as follows:
       imports = { all of the data from row exportLIDs.front() ;
                   ....
                   all of the data from row exportLIDs.back() }
      This doesn't have the best locality, but is necessary because
      the data for a Packet (all data associated with an LID) is
      required to be contiguous. */
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(as<size_t>(imports.size()) != getNumVectors()*importLIDs.size(), std::runtime_error,
        ": Imports buffer size must be consistent with the amount of data to be exported.");
#endif
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(as<size_t>(constantNumPackets) == as<size_t>(0), std::runtime_error,
        ": constantNumPackets input argument must be nonzero.");

    const size_t myStride = MVT::getStride(lclMV_);
    const size_t numVecs  = getNumVectors();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(as<size_t>(numPacketsPerLID.size()) != as<size_t>(importLIDs.size()), std::runtime_error,
        ": numPacketsPerLID must have the same length as importLIDs.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(as<size_t>(numVecs) != as<size_t>(constantNumPackets), std::runtime_error,
        ": constantNumPackets must equal numVecs.");

    if (numVecs > 0 && importLIDs.size()) {
      typename ArrayView<const       Scalar>::iterator impptr;
      typename ArrayView<const LocalOrdinal>::iterator  idptr;
      impptr = imports.begin();
      // NOTE (mfh 10 Mar 2012) If you want to implement custom
      // combine modes, start editing here.  Also, if you trust
      // inlining, it would be nice to condense this code by using a
      // binary function object f:
      //
      // ncview_[...] = f (ncview_[...], *impptr++);
      if (CM == INSERT || CM == REPLACE) {
        if (isConstantStride()) {
          for (idptr = importLIDs.begin(); idptr != importLIDs.end(); ++idptr) {
            for (size_t j = 0; j < numVecs; ++j) {
              ncview_[myStride*j + *idptr] = *impptr++;
            }
          }
        }
        else {
          for (idptr = importLIDs.begin(); idptr != importLIDs.end(); ++idptr) {
            for (size_t j = 0; j < numVecs; ++j) {
              ncview_[myStride*whichVectors_[j] + *idptr] = *impptr++;
            }
          }
        }
      }
      else if (CM == ADD) {
        if (isConstantStride()) {
          for (idptr = importLIDs.begin(); idptr != importLIDs.end(); ++idptr) {
            for (size_t j = 0; j < numVecs; ++j) {
              ncview_[myStride*j + *idptr] += *impptr++;
            }
          }
        }
        else {
          for (idptr = importLIDs.begin(); idptr != importLIDs.end(); ++idptr) {
            for (size_t j = 0; j < numVecs; ++j) {
              ncview_[myStride*whichVectors_[j] + *idptr] += *impptr++;
            }
          }
        }
      }
      else if (CM == ABSMAX) {
        if (isConstantStride()) {
          for (idptr = importLIDs.begin(); idptr != importLIDs.end(); ++idptr) {
            for (size_t j = 0; j < numVecs; ++j) {
              Scalar &curval       = ncview_[myStride*j + *idptr];
              const Scalar &newval = *impptr++;
              curval = std::max( SCT::magnitude(curval), SCT::magnitude(newval) );
            }
          }
        }
        else {
          for (idptr = importLIDs.begin(); idptr != importLIDs.end(); ++idptr) {
            for (size_t j = 0; j < numVecs; ++j) {
              Scalar &curval       = ncview_[myStride*whichVectors_[j] + *idptr];
              const Scalar &newval = *impptr++;
              curval = std::max( SCT::magnitude(curval), SCT::magnitude(newval) );
            }
          }
        }
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(CM != ADD && CM != REPLACE && CM != INSERT && CM != ABSMAX, std::invalid_argument,
            ": Invalid CombineMode: " << CM);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  inline size_t MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getNumVectors() const {
    if (isConstantStride()) {
      return MVT::getNumCols(lclMV_);
    }
    else {
      return whichVectors_.size();
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  dot (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A,
       const Teuchos::ArrayView<Scalar> &dots) const
  {
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::as;
    using Teuchos::arcp_const_cast;

    const std::string tfecfFuncName("dot()");
    const size_t myLen   = getLocalLength();
    const size_t numVecs = getNumVectors();
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( !this->getMap()->isCompatible(*A.getMap()), std::runtime_error,
        ": MultiVectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << *this->getMap()
        << "A.getMap(): " << std::endl << *A.getMap() << std::endl);
#else
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( getLocalLength() != A.getLocalLength(), std::runtime_error,
        ": MultiVectors do not have the same local length.");
#endif
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(A.getNumVectors() != numVecs, std::runtime_error,
        ": MultiVectors must have the same number of vectors.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(Teuchos::as<size_t>(dots.size()) != numVecs, std::runtime_error,
        ": dots.size() must be as large as the number of vectors in *this and A.");
    if (isConstantStride() && A.isConstantStride()) {
      MVT::Dot(lclMV_,A.lclMV_,dots);
    }
    else {
      KMV v(MVT::getNode(lclMV_)), a(MVT::getNode(lclMV_));
      ArrayRCP<Scalar> vptr  = arcp_const_cast<Scalar>(MVT::getValues(lclMV_));
      ArrayRCP<Scalar> avptr = arcp_const_cast<Scalar>(MVT::getValues(A.lclMV_));
      for (size_t j=0; j < numVecs; ++j) {
        ArrayRCP<Scalar> vj  =   getSubArrayRCP( vptr,j);
        ArrayRCP<Scalar> avj = A.getSubArrayRCP(avptr,j);
        MVT::initializeValues(a,myLen, 1, avj, myLen);
        MVT::initializeValues(v,myLen, 1,  vj, myLen);
        dots[j] = MVT::Dot((const KMV&)v,(const KMV &)a);
      }
    }
    if (this->isDistributed()) {
      Array<Scalar> ldots(dots);
      Teuchos::reduceAll(*this->getMap()->getComm(),Teuchos::REDUCE_SUM,as<int>(numVecs),ldots.getRawPtr(),dots.getRawPtr());
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  norm2 (const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& norms) const
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::as;
    using Teuchos::reduceAll;
    using Teuchos::REDUCE_SUM;
    using Teuchos::ScalarTraits;
    typedef typename ScalarTraits<Scalar>::magnitudeType MT;
    typedef ScalarTraits<MT> STM;

    const size_t numVecs = this->getNumVectors();
    TEUCHOS_TEST_FOR_EXCEPTION(as<size_t>(norms.size()) != numVecs,
      std::runtime_error,
      "Tpetra::MultiVector::norm2(norms): norms.size() must be as large as the number of vectors in *this.");
    if (isConstantStride ()) {
      MVT::Norm2Squared (lclMV_,norms);
    }
    else {
      KMV v (MVT::getNode (lclMV_));
      ArrayRCP<Scalar> vi;
      for (size_t i=0; i < numVecs; ++i) {
        vi = arcp_const_cast<Scalar> (MVT::getValues (lclMV_, whichVectors_[i]));
        MVT::initializeValues (v, MVT::getNumRows (lclMV_), 1, vi, MVT::getStride (lclMV_));
        norms[i] = MVT::Norm2Squared (v);
      }
    }
    if (this->isDistributed ()) {
      Array<MT> lnorms (norms);
      // FIXME (mfh 25 Apr 2012) Somebody explain to me why we're
      // calling Teuchos::reduceAll when MultiVector has a perfectly
      // good reduce() function.
      reduceAll (*this->getMap ()->getComm (), REDUCE_SUM, as<int> (numVecs),
                 lnorms.getRawPtr (), norms.getRawPtr ());
    }
    for (typename ArrayView<MT>::iterator n = norms.begin(); n != norms.begin()+numVecs; ++n) {
      (*n) = STM::squareroot (*n);
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  normWeighted (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& weights,
                const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::as;
    using Teuchos::reduceAll;
    using Teuchos::REDUCE_SUM;
    using Teuchos::ScalarTraits;
    typedef ScalarTraits<Scalar> SCT;
    typedef typename SCT::magnitudeType Mag;

    const std::string tfecfFuncName("normWeighted()");

    const Mag OneOverN = ScalarTraits<Mag>::one() / getGlobalLength();
    bool OneW = false;
    const size_t numVecs = this->getNumVectors();
    if (weights.getNumVectors() == 1) {
      OneW = true;
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(weights.getNumVectors() != numVecs, std::runtime_error,
        ": MultiVector of weights must contain either one vector or the same number of vectors as this.");
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( !this->getMap()->isCompatible(*weights.getMap()), std::runtime_error,
      ": MultiVectors do not have compatible Maps:" << std::endl
      << "this->getMap(): " << std::endl << *this->getMap()
      << "weights.getMap(): " << std::endl << *weights.getMap() << std::endl);
#else
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( getLocalLength() != weights.getLocalLength(), std::runtime_error,
      ": MultiVectors do not have the same local length.");
#endif
    const size_t myLen = getLocalLength();

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(as<size_t>(norms.size()) != numVecs, std::runtime_error,
      ": norms.size() must be as large as the number of vectors in *this.");
    if (isConstantStride() && weights.isConstantStride()) {
      MVT::WeightedNorm(lclMV_,weights.lclMV_,norms);
    }
    else {
      KMV v(MVT::getNode(lclMV_)), w(MVT::getNode(lclMV_));
      ArrayRCP<Scalar> vptr = arcp_const_cast<Scalar> (MVT::getValues (lclMV_));
      ArrayRCP<Scalar> wptr = arcp_const_cast<Scalar> (MVT::getValues (weights.lclMV_));
      ArrayRCP<Scalar> wj   = wptr.persistingView (0,myLen);
      MVT::initializeValues (w,myLen, 1, wj, myLen);
      for (size_t j=0; j < numVecs; ++j) {
        ArrayRCP<Scalar> vj = getSubArrayRCP (vptr, j);
        MVT::initializeValues(v,myLen, 1,  vj, myLen);
        if (! OneW) {
          wj = weights.getSubArrayRCP(wptr,j);
          MVT::initializeValues (w, myLen, 1, wj, myLen);
        }
        norms[j] = MVT::WeightedNorm ((const KMV&)v, (const KMV &)w);
      }
    }
    if (this->isDistributed ()) {
      Array<Mag> lnorms(norms);
      reduceAll (*this->getMap ()->getComm (), REDUCE_SUM, as<int>(numVecs),
                 lnorms.getRawPtr (), norms.getRawPtr ());
    }
    for (typename ArrayView<Mag>::iterator n = norms.begin(); n != norms.begin()+numVecs; ++n) {
      *n = ScalarTraits<Mag>::squareroot(*n * OneOverN);
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  norm1 (const ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const
  {
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::arcp_const_cast;
    using Teuchos::as;
    using Teuchos::reduceAll;
    using Teuchos::REDUCE_SUM;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Mag;

    const size_t numVecs = this->getNumVectors();
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(norms.size()) != numVecs, std::runtime_error,
        "Tpetra::MultiVector::norm1(norms): norms.size() must be as large as the number of vectors in *this.");
    if (isConstantStride()) {
      MVT::Norm1(lclMV_,norms);
    }
    else {
      KMV v(MVT::getNode(lclMV_));
      ArrayRCP<Scalar> vj;
      for (size_t j=0; j < numVecs; ++j) {
        vj = arcp_const_cast<Scalar> (MVT::getValues(lclMV_,whichVectors_[j]) );
        MVT::initializeValues (v, MVT::getNumRows(lclMV_), 1, vj, MVT::getStride (lclMV_));
        norms[j] = MVT::Norm1 ((const KMV&)v);
      }
    }
    if (this->isDistributed ()) {
      Array<Mag> lnorms (norms);
      reduceAll (*this->getMap ()->getComm (), REDUCE_SUM, as<int> (numVecs), lnorms.getRawPtr (), norms.getRawPtr ());
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  normInf (const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const
  {
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::arcp_const_cast;
    using Teuchos::as;
    using Teuchos::reduceAll;
    using Teuchos::REDUCE_MAX;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Mag;

    const size_t numVecs = this->getNumVectors();
    TEUCHOS_TEST_FOR_EXCEPTION(as<size_t>(norms.size()) != numVecs, std::runtime_error,
      "Tpetra::MultiVector::normInf(norms): norms.size() must be as large as the number of vectors in *this.");
    if (isConstantStride ()) {
      MVT::NormInf(lclMV_,norms);
    }
    else {
      KMV v(MVT::getNode(lclMV_));
      ArrayRCP<Scalar> vj;
      for (size_t j=0; j < numVecs; ++j) {
        vj = Teuchos::arcp_const_cast<Scalar>( MVT::getValues(lclMV_,whichVectors_[j]) );
        MVT::initializeValues(v,MVT::getNumRows(lclMV_), 1, vj, MVT::getStride(lclMV_));
        norms[j] = MVT::NormInf((const KMV&)v);
      }
    }
    if (this->isDistributed()) {
      Array<Mag> lnorms(norms);
      reduceAll (*this->getMap ()->getComm (), REDUCE_MAX, as<int> (numVecs), lnorms.getRawPtr (), norms.getRawPtr ());
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  meanValue (const Teuchos::ArrayView<Scalar> &means) const
  {
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::arcp_const_cast;
    using Teuchos::as;
    using Teuchos::arcp_const_cast;
    using Teuchos::reduceAll;
    using Teuchos::REDUCE_SUM;
    typedef Teuchos::ScalarTraits<Scalar> SCT;

    const size_t numVecs = getNumVectors();
    const size_t myLen   = getLocalLength();
    TEUCHOS_TEST_FOR_EXCEPTION(as<size_t>(means.size()) != numVecs, std::runtime_error,
      "Tpetra::MultiVector::meanValue(): means.size() must be as large as the number of vectors in *this.");
    // compute local components of the means
    // sum these across all nodes
    if (isConstantStride()) {
      MVT::Sum(lclMV_,means);
    }
    else {
      KMV v(MVT::getNode(lclMV_));
      ArrayRCP<Scalar> vptr = arcp_const_cast<Scalar>(MVT::getValues(lclMV_));
      for (size_t j=0; j < numVecs; ++j) {
        ArrayRCP<Scalar> vj =   getSubArrayRCP( vptr,j);
        MVT::initializeValues(v,myLen, 1,  vj, myLen);
        means[j] = MVT::Sum((const KMV &)v);
      }
    }
    if (this->isDistributed()) {
      Array<Scalar> lmeans(means);
      // only combine if we are a distributed MV
      reduceAll (*this->getMap ()->getComm (), REDUCE_SUM, as<int> (numVecs),
                 lmeans.getRawPtr (), means.getRawPtr ());
    }
    // mfh 12 Apr 2012: Don't take out the cast from the ordinal type
    // to the magnitude type, since operator/ (std::complex<T>, int)
    // isn't necessarily defined.
    const Scalar OneOverN =
      SCT::one() / as<typename SCT::magnitudeType> (getGlobalLength ());
    for (typename ArrayView<Scalar>::iterator i = means.begin();
         i != means.begin() + numVecs;
         ++i) {
      (*i) = (*i) * OneOverN;
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  randomize()
  {
    if (isConstantStride ()) {
      MVT::Random (lclMV_);
    }
    else {
      const size_t numVecs = this->getNumVectors ();
      KMV v (MVT::getNode (lclMV_));
      Teuchos::ArrayRCP<Scalar> vj;
      for (size_t j = 0; j < numVecs; ++j) {
        vj = MVT::getValuesNonConst (lclMV_, whichVectors_[j]);
        MVT::initializeValues (v, MVT::getNumRows (lclMV_), 1, vj, MVT::getStride (lclMV_));
        MVT::Random (v);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  putScalar (const Scalar &alpha)
  {
    const size_t numVecs = getNumVectors();
    if (isConstantStride()) {
      MVT::Init(lclMV_,alpha);
    }
    else {
      KMV v(MVT::getNode(lclMV_));
      Teuchos::ArrayRCP<Scalar> vj;
      for (size_t j=0; j < numVecs; ++j) {
        vj = MVT::getValuesNonConst(lclMV_,whichVectors_[j]);
        MVT::initializeValues(v,MVT::getNumRows(lclMV_), 1, vj, MVT::getStride(lclMV_));
        MVT::Init(v,alpha);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  replaceMap (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map)
  {
#ifdef HAVE_TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(! this->getMap ()->isCompatible (*map),
      std::invalid_argument, "Tpetra::MultiVector::replaceMap(): The input map "
      "is not compatible with this multivector's current map.  The replaceMap() "
      "method is not for data redistribution; use Import or Export for that.");
#endif // HAVE_TEUCHOS_DEBUG
    this->map_ = map;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  scale (const Scalar &alpha)
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;

    // NOTE: can't substitute putScalar(0.0) for scale(0.0), because
    //       the former will overwrite NaNs present in the MultiVector, while the
    //       semantics of this call require multiplying them by 0, which IEEE requires to be NaN
    const size_t numVecs = getNumVectors();
    if (alpha == Teuchos::ScalarTraits<Scalar>::one()) {
      // do nothing
    }
    else if (isConstantStride()) {
      MVT::Scale(lclMV_,alpha);
    }
    else {
      KMV v(MVT::getNode(lclMV_));
      ArrayRCP<Scalar> vj;
      for (size_t j=0; j < numVecs; ++j) {
        vj = arcp_const_cast<Scalar> (MVT::getValues(lclMV_, whichVectors_[j]));
        MVT::initializeValues (v, MVT::getNumRows (lclMV_), 1, vj, MVT::getStride (lclMV_));
        MVT::Scale (v, alpha);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  scale (Teuchos::ArrayView<const Scalar> alphas)
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;
    using Teuchos::as;

    const size_t numVecs = this->getNumVectors();
    TEUCHOS_TEST_FOR_EXCEPTION(as<size_t>(alphas.size()) != numVecs, std::runtime_error,
      "Tpetra::MultiVector::scale(alphas): alphas.size() must be as large as the number of vectors in *this.");
    KMV vec(MVT::getNode(lclMV_));
    const size_t myLen = MVT::getNumRows(lclMV_);
    if (myLen == 0) return;
    ArrayRCP<Scalar> mybuf = MVT::getValuesNonConst(lclMV_);
    for (size_t j = 0; j < numVecs; ++j) {
      if (alphas[j] == Teuchos::ScalarTraits<Scalar>::one()) {
        // do nothing: NaN * 1.0 == NaN, Number*1.0 == Number
      }
      else {
        ArrayRCP<Scalar> mybufj = getSubArrayRCP(mybuf,j);
        MVT::initializeValues(vec,myLen,1,mybufj,myLen);
        MVT::Scale(vec,alphas[j]);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  scale (const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A)
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;
    using Teuchos::as;

    const std::string tfecfFuncName("scale(alpha,A)");

    const size_t numVecs = getNumVectors(),
                 myLen   = getLocalLength();
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(! this->getMap ()->isCompatible (*A.getMap ()),
      std::runtime_error, ": MultiVectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << *(this->getMap ())
        << "A.getMap(): " << std::endl << *(A.getMap ()) << std::endl);
#else
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( getLocalLength() != A.getLocalLength(), std::runtime_error,
      ": MultiVectors do not have the same local length.");
#endif
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(A.getNumVectors() != numVecs, std::runtime_error,
      ": MultiVectors must have the same number of vectors.");
    if (isConstantStride() && A.isConstantStride()) {
      // set me == alpha*A
      MVT::Scale(lclMV_,alpha,(const KMV&)A.lclMV_);
    }
    else {
      KMV v(MVT::getNode(lclMV_)), a(MVT::getNode(lclMV_));
      ArrayRCP<Scalar> vptr  = MVT::getValuesNonConst (lclMV_);
      ArrayRCP<Scalar> avptr = arcp_const_cast<Scalar> (MVT::getValues (A.lclMV_));
      for (size_t j=0; j < numVecs; ++j) {
        ArrayRCP<Scalar>  vj =   getSubArrayRCP (vptr,  j);
        ArrayRCP<Scalar> avj = A.getSubArrayRCP (avptr, j);
        MVT::initializeValues (a,myLen, 1, avj, myLen);
        MVT::initializeValues (v,myLen, 1,  vj, myLen);
        MVT::Scale (v, alpha, (const KMV &)a);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  reciprocal (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A)
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;
    using Teuchos::as;

    const std::string tfecfFuncName("reciprocal()");

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( !this->getMap()->isCompatible(*A.getMap()), std::runtime_error,
        ": MultiVectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << *this->getMap()
        << "A.getMap(): " << std::endl << *A.getMap() << std::endl);
#else
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( getLocalLength() != A.getLocalLength(), std::runtime_error,
        ": MultiVectors do not have the same local length.");
#endif
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(A.getNumVectors() != this->getNumVectors(), std::runtime_error,
        ": MultiVectors must have the same number of vectors.");
    const size_t numVecs = getNumVectors();
    const size_t myLen = getLocalLength();
    try {
      if (isConstantStride() && A.isConstantStride()) {
        MVT::Recip(lclMV_,(const KMV&)A.lclMV_);
      }
      else {
        KMV v(MVT::getNode(lclMV_)), a(MVT::getNode(lclMV_));
        ArrayRCP<Scalar> vptr = MVT::getValuesNonConst(lclMV_),
                        avptr = arcp_const_cast<Scalar>(MVT::getValues(A.lclMV_));
        for (size_t j=0; j < numVecs; ++j) {
          ArrayRCP<Scalar> vj =   getSubArrayRCP( vptr,j),
                          avj = A.getSubArrayRCP(avptr,j);
          MVT::initializeValues(a,myLen, 1, avj, myLen);
          MVT::initializeValues(v,myLen, 1,  vj, myLen);
          MVT::Recip(v,(const KMV &)a);
        }
      }
    }
    catch (std::runtime_error &e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(true,std::runtime_error,
          ": caught exception from Kokkos:" << std::endl
          << e.what() << std::endl);
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::abs(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A) {
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;
    using Teuchos::as;

    const std::string tfecfFuncName("abs()");
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( !this->getMap()->isCompatible(*A.getMap()), std::runtime_error,
        ": MultiVectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << *this->getMap()
        << "A.getMap(): " << std::endl << *A.getMap() << std::endl);
#else
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( getLocalLength() != A.getLocalLength(), std::runtime_error,
        ": MultiVectors do not have the same local length.");
#endif
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(A.getNumVectors() != this->getNumVectors(), std::runtime_error,
        ": MultiVectors must have the same number of vectors.");
    const size_t myLen = getLocalLength();
    const size_t numVecs = getNumVectors();
    if (isConstantStride() && A.isConstantStride()) {
      MVT::Abs(lclMV_,(const KMV&)A.lclMV_);
    }
    else {
      KMV v(MVT::getNode(lclMV_)), a(MVT::getNode(lclMV_));
      ArrayRCP<Scalar> vptr = MVT::getValuesNonConst(lclMV_),
                      avptr = arcp_const_cast<Scalar>(MVT::getValues(A.lclMV_));
      for (size_t j=0; j < numVecs; ++j) {
        ArrayRCP<Scalar> vj =   getSubArrayRCP( vptr,j),
                        avj = A.getSubArrayRCP(avptr,j);
        MVT::initializeValues(a,myLen, 1, avj, myLen);
        MVT::initializeValues(v,myLen, 1,  vj, myLen);
        MVT::Abs(v,(const KMV &)a);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::update(
                      const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A,
                      const Scalar &beta)
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;
    using Teuchos::as;

    const std::string tfecfFuncName("update()");
    // this = beta*this + alpha*A
    // must support case where &this == &A
    // can't short circuit on alpha==0.0 or beta==0.0, because 0.0*NaN != 0.0
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( !this->getMap()->isCompatible(*A.getMap()), std::runtime_error,
        ": MultiVectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << this->getMap()
        << "A.getMap(): " << std::endl << *A.getMap() << std::endl);
#else
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( getLocalLength() != A.getLocalLength(), std::runtime_error,
        ": MultiVectors do not have the same local length.");
#endif
    const size_t myLen = getLocalLength();
    const size_t numVecs = getNumVectors();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(A.getNumVectors() != this->getNumVectors(), std::runtime_error,
        ": MultiVectors must have the same number of vectors.");
    if (isConstantStride() && A.isConstantStride()) {
      MVT::GESUM(lclMV_,alpha,(const KMV&)A.lclMV_,beta);
    }
    else {
      KMV v(MVT::getNode(lclMV_)), a(MVT::getNode(lclMV_));
      ArrayRCP<Scalar> vptr = MVT::getValuesNonConst(lclMV_),
                      avptr = arcp_const_cast<Scalar>(MVT::getValues(A.lclMV_));
      for (size_t j=0; j < numVecs; ++j) {
        ArrayRCP<Scalar> vj =   getSubArrayRCP( vptr,j),
                        avj = A.getSubArrayRCP(avptr,j);
        MVT::initializeValues(a,myLen, 1, avj, myLen);
        MVT::initializeValues(v,myLen, 1,  vj, myLen);
        MVT::GESUM(v,alpha,(const KMV &)a,beta);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::update(
                      const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A,
                      const Scalar &beta,  const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B,
                      const Scalar &gamma)
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;
    using Teuchos::as;

    const std::string tfecfFuncName("update()");
    // this = alpha*A + beta*B + gamma*this
    // must support case where &this == &A or &this == &B
    // can't short circuit on alpha==0.0 or beta==0.0 or gamma==0.0, because 0.0*NaN != 0.0
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( !this->getMap()->isCompatible(*A.getMap()) || !this->getMap()->isCompatible(*B.getMap()),
        std::runtime_error,
        ": MultiVectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << *this->getMap()
        << "A.getMap(): " << std::endl << *A.getMap() << std::endl
        << "B.getMap(): " << std::endl << *B.getMap() << std::endl);
#else
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( getLocalLength() != A.getLocalLength() || getLocalLength() != B.getLocalLength(), std::runtime_error,
        ": MultiVectors do not have the same local length.");
#endif
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(A.getNumVectors() != this->getNumVectors() || B.getNumVectors() != this->getNumVectors(), std::runtime_error,
        ": MultiVectors must have the same number of vectors.");
    const size_t myLen = getLocalLength();
    const size_t numVecs = getNumVectors();
    if (isConstantStride() && A.isConstantStride() && B.isConstantStride()) {
      MVT::GESUM(lclMV_,alpha,(const KMV&)A.lclMV_,beta,(const KMV&)B.lclMV_,gamma);
    }
    else {
      KMV v(MVT::getNode(lclMV_)), a(MVT::getNode(lclMV_)), b(MVT::getNode(lclMV_));
      ArrayRCP<Scalar> vptr = MVT::getValuesNonConst(lclMV_),
                      avptr = arcp_const_cast<Scalar>(MVT::getValues(A.lclMV_)),
                      bvptr = arcp_const_cast<Scalar>(MVT::getValues(B.lclMV_));
      for (size_t j=0; j < numVecs; ++j) {
        ArrayRCP<Scalar> vj =   getSubArrayRCP( vptr,j),
                        avj = A.getSubArrayRCP(avptr,j),
                        bvj = B.getSubArrayRCP(bvptr,j);
        MVT::initializeValues(b,myLen, 1, bvj, myLen);
        MVT::initializeValues(a,myLen, 1, avj, myLen);
        MVT::initializeValues(v,myLen, 1,  vj, myLen);
        MVT::GESUM(v,alpha,(const KMV&)a,beta,(const KMV&)b,gamma);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<const Scalar>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  getData (size_t j) const
  {
    Teuchos::RCP<Node> node = MVT::getNode(lclMV_);
    KOKKOS_NODE_TRACE("MultiVector::getData()")
    return node->template viewBuffer<Scalar> (getLocalLength (), getSubArrayRCP (MVT::getValues(lclMV_), j));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<Scalar>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  getDataNonConst(size_t j)
  {
    Teuchos::RCP<Node> node = MVT::getNode(lclMV_);
    KOKKOS_NODE_TRACE("MultiVector::getDataNonConst()")
    return node->template viewBufferNonConst<Scalar> (Kokkos::ReadWrite, getLocalLength (), getSubArrayRCP (MVT::getValuesNonConst (lclMV_), j));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>&
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  operator= (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &source)
  {
    const std::string tfecfFuncName("operator=()");
    // Check for special case of this=Source, in which case we do nothing.
    if (this != &source) {
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( !this->getMap()->isCompatible(*source.getMap()), std::invalid_argument,
          ": MultiVectors do not have compatible Maps:" << std::endl
          << "this->getMap(): " << std::endl << *this->getMap()
          << "source.getMap(): " << std::endl << *source.getMap() << std::endl);
#else
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( getLocalLength() != source.getLocalLength(), std::invalid_argument,
          ": MultiVectors do not have the same local length.");
#endif
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(source.getNumVectors() != getNumVectors(), std::invalid_argument,
          ": MultiVectors must have the same number of vectors.");
      Teuchos::RCP<Node> node = MVT::getNode (lclMV_);
      const size_t numVecs = getNumVectors();
      if (isConstantStride() && source.isConstantStride() && getLocalLength()==getStride() && source.getLocalLength()==source.getStride()) {
        // Both multivectors' data are stored contiguously, so we can
        // copy in one call.
        KOKKOS_NODE_TRACE("MultiVector::operator=()")
        node->template copyBuffers<Scalar>(getLocalLength()*numVecs, MVT::getValues(source.lclMV_), MVT::getValuesNonConst(lclMV_) );
      }
      else {
        // We have to copy the columns one at a time.
        for (size_t j=0; j < numVecs; ++j) {
          KOKKOS_NODE_TRACE("MultiVector::operator=()")
          node->template copyBuffers<Scalar>(getLocalLength(), source.getSubArrayRCP(MVT::getValues(source.lclMV_),j),
                                                                      getSubArrayRCP(MVT::getValuesNonConst(lclMV_),j) );
        }
      }
    }
    return(*this);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  subCopy (const Teuchos::ArrayView<const size_t> &cols) const
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    TEUCHOS_TEST_FOR_EXCEPTION(cols.size() < 1, std::runtime_error,
        "Tpetra::MultiVector::subCopy(cols): cols must contain at least one column.");
    size_t numCopyVecs = cols.size();
    const bool zeroData = false;
    RCP<Node> node = MVT::getNode(lclMV_);
    RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > mv;
    // mv is allocated with constant stride
    mv = rcp (new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> (this->getMap (), numCopyVecs, zeroData));
    // copy data from *this into mv
    for (size_t j=0; j<numCopyVecs; ++j) {
      KOKKOS_NODE_TRACE("MultiVector::subCopy()")
      node->template copyBuffers<Scalar> (getLocalLength (),
                                          getSubArrayRCP (MVT::getValues (lclMV_), cols[j]),
                                          MVT::getValuesNonConst (mv->lclMV_, j));
    }
    return mv;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  subCopy (const Teuchos::Range1D &colRng) const
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    TEUCHOS_TEST_FOR_EXCEPTION(colRng.size() == 0, std::runtime_error,
        "Tpetra::MultiVector::subCopy(Range1D): range must include at least one vector.");
    size_t numCopyVecs = colRng.size();
    const bool zeroData = false;
    RCP<Node> node = MVT::getNode(lclMV_);
    RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > mv;
    // mv is allocated with constant stride
    mv = rcp (new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> (this->getMap (), numCopyVecs, zeroData));
    // copy data from *this into mv
    for (size_t js=colRng.lbound(), jd=0; jd<numCopyVecs; ++jd, ++js) {
      KOKKOS_NODE_TRACE("MultiVector::subCopy()")
      node->template copyBuffers<Scalar> (getLocalLength (),
                                          getSubArrayRCP (MVT::getValues (lclMV_), js),
                                          MVT::getValuesNonConst (mv->lclMV_, jd));
    }
    return mv;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  offsetView (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& subMap,
              size_t offset) const
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::rcp;

    typedef const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> CMV;
    TEUCHOS_TEST_FOR_EXCEPTION( subMap->getNodeNumElements() + offset > this->getLocalLength(), std::runtime_error,
        "Tpetra::MultiVector::offsetView(subMap,offset): sizes are not sane.\noffset == " << offset << "\nsubMap: " << subMap->description() << "\nthis->rowMap: " << this->getMap()->description());
    const size_t numVecs = this->getNumVectors(),
                myStride = MVT::getStride(lclMV_),
                  newLen = subMap->getNodeNumElements();
    // this is const, so the lclMV_ is const, so that we can only get const buffers
    // we will cast away the const; this is okay, because
    //   a) the constructor doesn't modify the data, and
    //   b) we are encapsulating in a const MV before returning
    ArrayRCP<const Scalar> cbuf = MVT::getValues(lclMV_);
    ArrayRCP<Scalar>      ncbuf = arcp_const_cast<Scalar> (cbuf);
    RCP<CMV> constViewMV;
    if (isConstantStride()) {
      // view goes from first entry of first vector to last entry of last vector
      ArrayRCP<Scalar> subdata = ncbuf.persistingView( offset, myStride * (numVecs-1) + newLen );
      constViewMV = rcp (new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> (subMap, subdata, myStride, numVecs, COMPUTE_VIEW_CONSTRUCTOR));
    }
    else {
      // use same which index, but with an offset start pointer
      size_t maxSubVecIndex = *std::max_element(whichVectors_.begin(), whichVectors_.end());
      ArrayRCP<Scalar> subdata = ncbuf.persistingView( offset, myStride * maxSubVecIndex + newLen );
      constViewMV = rcp (new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> (subMap, subdata, myStride, whichVectors_, COMPUTE_VIEW_CONSTRUCTOR));
    }
    return constViewMV;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  offsetViewNonConst (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& subMap,
                      size_t offset)
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::rcp;

    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
    TEUCHOS_TEST_FOR_EXCEPTION( subMap->getNodeNumElements() + offset > this->getLocalLength(), std::runtime_error,
      "Tpetra::MultiVector::offsetView(subMap,offset): sizes are not sane.\noffset == "
      << offset << "\nsubMap: " << subMap->description() << "\nthis->rowMap: "
      << this->getMap()->description());
    const size_t numVecs = this->getNumVectors(),
                myStride = MVT::getStride(lclMV_),
                  newLen = subMap->getNodeNumElements();
    // this is const, so the lclMV_ is const, so that we can only get const buffers
    // we will cast away the const; this is okay, because
    //   a) the constructor doesn't modify the data, and
    //   b) we are encapsulating in a const MV before returning
    ArrayRCP<Scalar> buf = MVT::getValuesNonConst(lclMV_);
    RCP<MV> subViewMV;
    if (isConstantStride()) {
      // view goes from first entry of first vector to last entry of last vector
      ArrayRCP<Scalar> subdata = buf.persistingView( offset, myStride * (numVecs-1) + newLen );
      subViewMV = rcp (new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> (subMap, subdata, myStride, numVecs, COMPUTE_VIEW_CONSTRUCTOR));
    }
    else {
      // use same which index, but with an offset start pointer
      size_t maxSubVecIndex = *std::max_element(whichVectors_.begin(), whichVectors_.end());
      ArrayRCP<Scalar> subdata = buf.persistingView( offset, myStride * maxSubVecIndex + newLen );
      subViewMV = rcp (new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> (subMap, subdata, myStride, whichVectors_, COMPUTE_VIEW_CONSTRUCTOR));
    }
    return subViewMV;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  subView (const ArrayView<const size_t> &cols) const
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_const_cast;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

    TEUCHOS_TEST_FOR_EXCEPTION(cols.size() == 0, std::runtime_error,
      "Tpetra::MultiVector::subView(ArrayView): range must include at least one vector.");
    // this is const, so the lclMV_ is const, so that we can only get const buffers
    // we will cast away the const; this is okay, because
    //   a) the constructor doesn't modify the data, and
    //   b) we are encapsulating in a const MV before returning
    const size_t myStride = MVT::getStride(lclMV_),
                 myLen    = MVT::getNumRows(lclMV_),
              numViewCols = cols.size();
    // use the smallest view possible of the buffer: from the first element of the minInd vector to the last element of the maxInd vector
    // this minimizes overlap between views, and keeps view of the minimum amount necessary, in order to allow node to achieve maximum efficiency.
    // adjust the indices appropriately; shift so that smallest index is 0
    ArrayRCP<const Scalar> cbuf = MVT::getValues(lclMV_);
    ArrayRCP<Scalar>      ncbuf = arcp_const_cast<Scalar> (cbuf);
    Array<size_t> newCols (numViewCols);
    size_t minInd = Teuchos::OrdinalTraits<size_t>::max();
    size_t maxInd = Teuchos::OrdinalTraits<size_t>::zero();
    if (isConstantStride()) {
      for (size_t j=0; j < numViewCols; ++j) {
        newCols[j] = cols[j];
        if (newCols[j] < minInd) minInd = newCols[j];
        if (maxInd < newCols[j]) maxInd = newCols[j];
      }
    }
    else {
      for (size_t j=0; j < numViewCols; ++j) {
        newCols[j] = whichVectors_[cols[j]];
        if (newCols[j] < minInd) minInd = newCols[j];
        if (maxInd < newCols[j]) maxInd = newCols[j];
      }
    }
    ArrayRCP<Scalar> minbuf = ncbuf.persistingView(minInd * myStride, myStride * (maxInd - minInd) + myLen);
    for (size_t j=0; j < numViewCols; ++j) {
      newCols[j] -= minInd;
    }
    RCP<const MV> constViewMV;
    return rcp_const_cast<const MV> (rcp (new MV (this->getMap (), minbuf, myStride, newCols (), COMPUTE_VIEW_CONSTRUCTOR)));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  subView (const Teuchos::Range1D &colRng) const
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

    TEUCHOS_TEST_FOR_EXCEPTION(colRng.size() == 0, std::runtime_error,
      "Tpetra::MultiVector::subView(Range1D): range must include at least one vector.");
    size_t numViewVecs = colRng.size();
    // this is const, so the lclMV_ is const, so that we can only get const buffers
    // we will cast away the const; this is okay, because
    //   a) the constructor doesn't modify the data, and
    //   b) we are encapsulating in a const MV before returning
    ArrayRCP<const Scalar> cbuf = MVT::getValues(lclMV_);
    ArrayRCP<Scalar>      ncbuf = arcp_const_cast<Scalar>(cbuf);
    // resulting MultiVector is constant stride only if *this is
    if (isConstantStride()) {
      // view goes from first entry of first vector to last entry of last vector
      ArrayRCP<Scalar> subdata = ncbuf.persistingView( MVT::getStride(lclMV_) * colRng.lbound(),
                                                       MVT::getStride(lclMV_) * (numViewVecs-1) + getLocalLength() );
      return rcp (new MV (this->getMap (), subdata, MVT::getStride (lclMV_), numViewVecs, COMPUTE_VIEW_CONSTRUCTOR));
    }
    // otherwise, use a subset of this whichVectors_ to construct new multivector
    Array<size_t> whchvecs( whichVectors_.begin()+colRng.lbound(), whichVectors_.begin()+colRng.ubound()+1 );
    return rcp (new MV (this->getMap (), ncbuf, MVT::getStride (lclMV_), whchvecs, COMPUTE_VIEW_CONSTRUCTOR));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  subViewNonConst (const ArrayView<const size_t> &cols)
  {
    using Teuchos::as;
    using Teuchos::arcp_const_cast;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

    TEUCHOS_TEST_FOR_EXCEPTION(cols.size() == 0, std::runtime_error,
      "Tpetra::MultiVector::subViewNonConst(ArrayView): range must include at least one vector.");
    if (isConstantStride()) {
      return rcp (new MV (this->getMap (), MVT::getValuesNonConst (lclMV_), MVT::getStride (lclMV_), cols, COMPUTE_VIEW_CONSTRUCTOR));
    }
    // else, lookup current whichVectors_ using cols
    Array<size_t> newcols(cols.size());
    for (size_t j = 0; j < as<size_t> (cols.size ()); ++j) {
      newcols[j] = whichVectors_[cols[j]];
    }
    return rcp (new MV (this->getMap (), MVT::getValuesNonConst (lclMV_), MVT::getStride (lclMV_), newcols (), COMPUTE_VIEW_CONSTRUCTOR));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  subViewNonConst (const Teuchos::Range1D &colRng)
  {
    using Teuchos::as;
    using Teuchos::arcp_const_cast;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

    TEUCHOS_TEST_FOR_EXCEPTION(colRng.size() == 0, std::runtime_error,
      "Tpetra::MultiVector::subViewNonConst(Range1D): range must include at least one vector.");
    size_t numViewVecs = colRng.size();
    // resulting MultiVector is constant stride only if *this is
    if (isConstantStride ()) {
      // view goes from first entry of first vector to last entry of last vector
      const size_t stride = MVT::getStride(lclMV_);
      ArrayRCP<Scalar> data = MVT::getValuesNonConst(lclMV_);
      ArrayRCP<Scalar> subdata = data.persistingView( stride * colRng.lbound(),
                                                      stride * (numViewVecs-1) + getLocalLength() );
      return rcp (new MV (this->getMap (), subdata, stride, numViewVecs, COMPUTE_VIEW_CONSTRUCTOR));
    }
    // otherwise, use a subset of this whichVectors_ to construct new multivector
    Array<size_t> whchvecs( whichVectors_.begin()+colRng.lbound(), whichVectors_.begin()+colRng.ubound()+1 );
    const size_t stride = MVT::getStride(lclMV_);
    ArrayRCP<Scalar> data = MVT::getValuesNonConst(lclMV_);
    return rcp (new MV (this->getMap (), data, stride, whchvecs, COMPUTE_VIEW_CONSTRUCTOR));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  getVector (size_t j) const
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;
    using Teuchos::rcp;
    using Teuchos::rcp_const_cast;
    typedef Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> V;

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( vectorIndexOutOfRange(j), std::runtime_error,
        "Tpetra::MultiVector::getVector(j): index j (== " << j << ") exceeds valid column range for this multivector.");
#endif
    // this is const, so lclMV_ is const, so we get const buff
    // it is safe to cast away the const because we will wrap it in a const Vector below
    ArrayRCP<Scalar> ncbuff;
    if (getLocalLength() > 0) {
      ArrayRCP<const Scalar> cbuff = getSubArrayRCP(MVT::getValues(lclMV_),j);
      ncbuff = arcp_const_cast<Scalar>(cbuff);
    }
    return rcp_const_cast<const V> (rcp (new V (this->getMap (), ncbuff, COMPUTE_VIEW_CONSTRUCTOR)));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getVectorNonConst(size_t j)
  {
    using Teuchos::ArrayRCP;
    using Teuchos::rcp;
    typedef Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> V;

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( vectorIndexOutOfRange(j), std::runtime_error,
        "Tpetra::MultiVector::getVectorNonConst(j): index j (== " << j << ") exceeds valid column range for this multivector.");
#endif
    ArrayRCP<Scalar> ncbuff;
    if (getLocalLength() > 0) {
      ncbuff = getSubArrayRCP (MVT::getValuesNonConst (lclMV_), j);
    }
    return rcp (new V (this->getMap (), ncbuff, COMPUTE_VIEW_CONSTRUCTOR));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  get1dCopy (Teuchos::ArrayView<Scalar> A, size_t LDA) const
  {
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::RCP;
    using Teuchos::rcp;

    const std::string tfecfFuncName("get1dCopy(A,LDA)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(LDA < getLocalLength(), std::runtime_error,
      ": specified stride is not large enough for local vector length.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(Teuchos::as<size_t>(A.size()) < LDA*(getNumVectors()-1)+getLocalLength(), std::runtime_error,
      ": specified stride/storage is not large enough for the number of vectors.");
    RCP<Node> node = MVT::getNode(lclMV_);
    const size_t myStride = MVT::getStride(lclMV_),
                  numCols = getNumVectors(),
                  myLen   = getLocalLength();
    if (myLen > 0) {
      ArrayRCP<const Scalar> mydata = MVT::getValues(lclMV_);
      KOKKOS_NODE_TRACE("MultiVector::getVectorNonConst()")
      ArrayRCP<const Scalar> myview = node->template viewBuffer<Scalar>(myStride*(numCols-1)+myLen,mydata);
      typename ArrayView<Scalar>::iterator Aptr = A.begin();
      for (size_t j=0; j<numCols; j++) {
        ArrayRCP<const Scalar> myviewj = getSubArrayRCP(myview,j);
        std::copy(myviewj,myviewj+myLen,Aptr);
        Aptr += LDA;
      }
      myview = Teuchos::null;
      mydata = Teuchos::null;
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  get2dCopy (Teuchos::ArrayView<const Teuchos::ArrayView<Scalar> > ArrayOfPtrs) const
  {
    using Teuchos::as;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::RCP;
    using Teuchos::rcp;

    const std::string tfecfFuncName("get2dCopy(ArrayOfPtrs)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(Teuchos::as<size_t>(ArrayOfPtrs.size()) != getNumVectors(), std::runtime_error,
        ": Array of pointers must contain as many pointers as the MultiVector has rows.");
    RCP<Node> node = MVT::getNode(lclMV_);
    const size_t numCols = getNumVectors(),
                   myLen = getLocalLength();
    if (myLen > 0) {
      ArrayRCP<const Scalar> mybuff = MVT::getValues(lclMV_);
      KOKKOS_NODE_TRACE("MultiVector::get2dCopy()")
      ArrayRCP<const Scalar> myview = node->template viewBuffer<Scalar>(mybuff.size(), mybuff);
      for (size_t j=0; j<numCols; ++j) {
#ifdef HAVE_TPETRA_DEBUG
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(as<size_t>(ArrayOfPtrs[j].size()) != getLocalLength(), std::runtime_error,
          ": The ArrayView provided in ArrayOfPtrs[" << j << "] was not large enough to contain the local entries.");
#endif
        ArrayRCP<const Scalar> myviewj = getSubArrayRCP(myview,j);
        std::copy(myviewj,myviewj+myLen,ArrayOfPtrs[j].begin());
      }
      myview = Teuchos::null;
      mybuff = Teuchos::null;
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<const Scalar>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::get1dView () const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(!isConstantStride(), std::runtime_error,
      "Tpetra::MultiVector::get1dView() requires that this MultiVector have constant stride.");
    Teuchos::RCP<Node> node = MVT::getNode(lclMV_);
    KOKKOS_NODE_TRACE("MultiVector::get1dView()")
    return node->template viewBuffer<Scalar>( getStride()*(getNumVectors()-1)+getLocalLength(), MVT::getValues(lclMV_) );
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<Scalar>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::get1dViewNonConst ()
  {
    TEUCHOS_TEST_FOR_EXCEPTION(!isConstantStride(), std::runtime_error,
      "Tpetra::MultiVector::get1dViewNonConst(): requires that this MultiVector have constant stride.");
    Teuchos::RCP<Node> node = MVT::getNode(lclMV_);
    KOKKOS_NODE_TRACE("MultiVector::get1dViewNonConst()")
    return node->template viewBufferNonConst<Scalar>( Kokkos::ReadWrite, getStride()*(getNumVectors()-1)+getLocalLength(), MVT::getValuesNonConst(lclMV_) );
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::get2dViewNonConst()
  {
    using Teuchos::arcp;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;

    RCP<Node> node = MVT::getNode(lclMV_);
    ArrayRCP<ArrayRCP<Scalar> > views = arcp<ArrayRCP<Scalar> > (getNumVectors ());
    if (isConstantStride ()) {
      const size_t myStride = MVT::getStride(lclMV_),
                    numCols = getNumVectors(),
                    myLen   = getLocalLength();
      if (myLen > 0) {
        KOKKOS_NODE_TRACE("MultiVector::get2dViewNonConst()")
        ArrayRCP<Scalar> myview = node->template viewBufferNonConst<Scalar>(Kokkos::ReadWrite,myStride*(numCols-1)+myLen,MVT::getValuesNonConst(lclMV_));
        for (size_t j=0; j<numCols; ++j) {
          views[j] = myview.persistingView(0,myLen);
          myview += myStride;
        }
      }
    }
    else {
      const size_t myStride = MVT::getStride(lclMV_),
                    numCols = MVT::getNumCols(lclMV_),
                     myCols = getNumVectors(),
                     myLen  = MVT::getNumRows(lclMV_);
      if (myLen > 0) {
        KOKKOS_NODE_TRACE("MultiVector::get2dViewNonConst()")
        ArrayRCP<Scalar> myview = node->template viewBufferNonConst<Scalar>(Kokkos::ReadWrite,myStride*(numCols-1)+myLen,MVT::getValuesNonConst(lclMV_));
        for (size_t j=0; j<myCols; ++j) {
          views[j] = myview.persistingView(whichVectors_[j]*myStride,myLen);
        }
      }
    }
    return views;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::get2dView() const
  {
    using Teuchos::arcp;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;

    RCP<Node> node = MVT::getNode(lclMV_);
    ArrayRCP< ArrayRCP<const Scalar> > views = arcp<ArrayRCP<const Scalar> >(getNumVectors());
    if (isConstantStride()) {
      const size_t myStride = MVT::getStride(lclMV_),
                    numCols = getNumVectors(),
                    myLen   = getLocalLength();
      if (myLen > 0) {
        KOKKOS_NODE_TRACE("MultiVector::get2dView()")
        ArrayRCP<const Scalar> myview = node->template viewBuffer<Scalar>(myStride*(numCols-1)+myLen,MVT::getValues(lclMV_));
        for (size_t j=0; j<numCols; ++j) {
          views[j] = myview.persistingView(0,myLen);
          myview += myStride;
        }
      }
    }
    else {
      const size_t myStride = MVT::getStride(lclMV_),
                    numCols = MVT::getNumCols(lclMV_),
                     myCols = getNumVectors(),
                     myLen  = MVT::getNumRows(lclMV_);
      if (myLen > 0) {
        KOKKOS_NODE_TRACE("MultiVector::get2dView()")
        ArrayRCP<const Scalar> myview = node->template viewBuffer<Scalar>(myStride*(numCols-1)+myLen,MVT::getValues(lclMV_));
        for (size_t j=0; j<myCols; ++j) {
          views[j] = myview.persistingView(whichVectors_[j]*myStride,myLen);
        }
      }
    }
    return views;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  multiply (Teuchos::ETransp transA,
            Teuchos::ETransp transB,
            const Scalar &alpha,
            const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A,
            const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& B,
            const Scalar &beta)
  {
    using Teuchos::NO_TRANS;      // enums
    using Teuchos::TRANS;
    using Teuchos::CONJ_TRANS;
    using Teuchos::null;
    using Teuchos::ScalarTraits;  // traits
    using Teuchos::as;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

    // This routine performs a variety of matrix-matrix multiply operations, interpreting
    // the MultiVector (this-aka C , A and B) as 2D matrices.  Variations are due to
    // the fact that A, B and C can be local replicated or global distributed
    // MultiVectors and that we may or may not operate with the transpose of
    // A and B.  Possible cases are:
    //                                       Num
    //      OPERATIONS                        cases  Notes
    //  1) C(local) = A^X(local) * B^X(local)  4    (X=Trans or Not, No comm needed)
    //  2) C(local) = A^T(distr) * B  (distr)  1    (2D dot product, replicate C)
    //  3) C(distr) = A  (distr) * B^X(local)  2    (2D vector update, no comm needed)
    //
    // The following operations are not meaningful for 1D distributions:
    //
    // u1) C(local) = A^T(distr) * B^T(distr)  1
    // u2) C(local) = A  (distr) * B^X(distr)  2
    // u3) C(distr) = A^X(local) * B^X(local)  4
    // u4) C(distr) = A^X(local) * B^X(distr)  4
    // u5) C(distr) = A^T(distr) * B^X(local)  2
    // u6) C(local) = A^X(distr) * B^X(local)  4
    // u7) C(distr) = A^X(distr) * B^X(local)  4
    // u8) C(local) = A^X(local) * B^X(distr)  4
    //
    // Total of 32 case (2^5).

    const std::string errPrefix("Tpetra::MultiVector::multiply(transOpA,transOpB,alpha,A,B,beta): ");

    TEUCHOS_TEST_FOR_EXCEPTION( ScalarTraits<Scalar>::isComplex && (transA == TRANS || transB == TRANS), std::invalid_argument,
        errPrefix << "non-conjugate transpose not supported for complex types.");
    transA = (transA == NO_TRANS ? NO_TRANS : CONJ_TRANS);
    transB = (transB == NO_TRANS ? NO_TRANS : CONJ_TRANS);

    // Compute effective dimensions, w.r.t. transpose operations on
    size_t A_nrows = (transA==CONJ_TRANS) ? A.getNumVectors() : A.getLocalLength();
    size_t A_ncols = (transA==CONJ_TRANS) ? A.getLocalLength() : A.getNumVectors();
    size_t B_nrows = (transB==CONJ_TRANS) ? B.getNumVectors() : B.getLocalLength();
    size_t B_ncols = (transB==CONJ_TRANS) ? B.getLocalLength() : B.getNumVectors();

    Scalar beta_local = beta; // local copy of beta; might be reassigned below

    TEUCHOS_TEST_FOR_EXCEPTION( getLocalLength() != A_nrows || getNumVectors() != B_ncols || A_ncols != B_nrows, std::runtime_error,
        errPrefix << "dimension of *this, op(A) and op(B) must be consistent.");

    bool A_is_local = !A.isDistributed();
    bool B_is_local = !B.isDistributed();
    bool C_is_local = !this->isDistributed();
    bool Case1 = ( C_is_local &&  A_is_local &&  B_is_local);                                           // Case 1: C(local) = A^X(local) * B^X(local)
    bool Case2 = ( C_is_local && !A_is_local && !B_is_local && transA==CONJ_TRANS && transB==NO_TRANS); // Case 2: C(local) = A^T(distr) * B  (distr)
    bool Case3 = (!C_is_local && !A_is_local &&  B_is_local && transA==NO_TRANS  );                     // Case 3: C(distr) = A  (distr) * B^X(local)

    // Test that we are considering a meaningful cases
    TEUCHOS_TEST_FOR_EXCEPTION( !Case1 && !Case2 && !Case3, std::runtime_error,
        errPrefix << "multiplication of op(A) and op(B) into *this is not a supported use case.");

    if (beta != ScalarTraits<Scalar>::zero() && Case2)
    {
      // if Case2, then C is local and contributions must be summed across all nodes
      // however, if beta != 0, then accumulate beta*C into the sum
      // when summing across all nodes, we only want to accumulate this once, so
      // set beta == 0 on all nodes except node 0
      int MyPID = this->getMap()->getComm()->getRank();
      if (MyPID!=0) beta_local = ScalarTraits<Scalar>::zero();
    }

    // Check if A, B, C have constant stride, if not then make temp copy (strided)
    RCP<const MV> Atmp, Btmp;
    RCP<MV>       Ctmp;
    if (isConstantStride() == false) Ctmp = rcp (new MV (*this));
    else Ctmp = rcp(this,false);

    if (A.isConstantStride() == false) Atmp = rcp (new MV (A));
    else Atmp = rcp(&A,false);

    if (B.isConstantStride() == false) Btmp = rcp (new MV (B));
    else Btmp = rcp(&B,false);

#ifdef HAVE_TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(!Ctmp->isConstantStride() || !Btmp->isConstantStride() || !Atmp->isConstantStride(), std::logic_error,
        errPrefix << "failed making temporary strided copies of input multivectors.");
#endif

    KMV &C_mv = Ctmp->lclMV_;
    {
      // get local multivectors
      const KMV &A_mv = Atmp->lclMV_;
      const KMV &B_mv = Btmp->lclMV_;
      // do the multiply (GEMM)
      MVT::GEMM(C_mv,transA,transB,alpha,A_mv,B_mv,beta_local);
    }

    // Dispose of (possibly) extra copies of A, B
    Atmp = null;
    Btmp = null;

    RCP<Node> node = MVT::getNode(lclMV_);
    // If *this was not strided, copy the data from the strided version and then delete it
    if (! isConstantStride ()) {
      // *this is not strided, we must put data from Ctmp into *this
      TEUCHOS_TEST_FOR_EXCEPT(&C_mv != &lclMV_);
      const size_t numVecs = MVT::getNumCols(lclMV_);
      for (size_t j=0; j < numVecs; ++j) {
        KOKKOS_NODE_TRACE("MultiVector::multiply()")
        node->template copyBuffers<Scalar>(getLocalLength(),MVT::getValues(C_mv,j),MVT::getValuesNonConst(lclMV_,whichVectors_[j]));
      }
    }

    // If Case 2 then sum up *this and distribute it to all processors.
    if (Case2) {
      this->reduce();
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  elementWiseMultiply (Scalar scalarAB,
                       const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A,
                       const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& B,
                       Scalar scalarThis)
  {
    using Teuchos::arcp_const_cast;

    const std::string tfecfFuncName("elementWiseMultiply()");

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( !this->getMap()->isCompatible(*A.getMap()) || !this->getMap()->isCompatible(*B.getMap()), std::runtime_error,
      ": MultiVectors do not have compatible Maps:" << std::endl
      << "this->getMap(): " << std::endl << *this->getMap()
      << "A.getMap(): " << std::endl << *A.getMap() << std::endl
      << "B.getMap(): " << std::endl << *B.getMap() << std::endl);
#else
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( getLocalLength() != A.getLocalLength() || getLocalLength() != B.getLocalLength(), std::runtime_error,
        ": MultiVectors do not have the same local length.");
#endif
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(B.getNumVectors() != this->getNumVectors(), std::runtime_error,
        ": MultiVectors 'this' and B must have the same number of vectors.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(A.getNumVectors() != 1, std::runtime_error,
        ": MultiVectors A must have just 1 vector.");
    try {
      MVT::ElemMult(lclMV_,scalarThis,scalarAB,(const KMV &)A.lclMV_,(const KMV &)B.lclMV_);
    }
    catch (std::runtime_error &e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(true,std::runtime_error,
          ": caught exception from Kokkos:" << std::endl
          << e.what() << std::endl);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::reduce()
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    using Teuchos::as;
    using Teuchos::RCP;
    using Teuchos::reduceAll;
    using Teuchos::REDUCE_SUM;

    // This method should be called only for locally replicated
    // MultiVectors (! isDistributed ()).
    TEUCHOS_TEST_FOR_EXCEPTION(this->isDistributed (), std::runtime_error,
      "Tpetra::MultiVector::reduce() should only be called with locally "
      "replicated or otherwise not distributed MultiVector objects.");
    RCP<const Teuchos::Comm<int> > comm = this->getMap()->getComm();
    if (comm->getSize() == 1) return;
    RCP<Node> node = MVT::getNode(lclMV_);
    // sum the data across all multivectors
    // need to have separate packed buffers for send and receive
    // if we're packed, we'll just set the receive buffer as our data, the send as a copy
    // if we're not packed, we'll use allocated buffers for both.
    ArrayView<Scalar> target;
    const size_t myStride = MVT::getStride(lclMV_),
                 numCols  = MVT::getNumCols(lclMV_),
                 myLen    = MVT::getNumRows(lclMV_);
    Array<Scalar> sourceBuffer(numCols*myLen), tmparr(0);
    bool packed = isConstantStride() && (myStride == myLen);
    KOKKOS_NODE_TRACE("MultiVector::reduce()")
    ArrayRCP<Scalar> bufView = node->template viewBufferNonConst<Scalar>(
                                          Kokkos::ReadWrite,myStride*(numCols-1)+myLen,
                                          MVT::getValuesNonConst(lclMV_) );
    if (packed) {
      // copy data from view to sourceBuffer, reduce into view below
      target = bufView(0,myLen*numCols);
      std::copy(target.begin(),target.end(),sourceBuffer.begin());
    }
    else {
      // copy data into sourceBuffer, reduce from sourceBuffer into tmparr below, copy back to view after that
      tmparr.resize(myLen*numCols);
      Scalar *sptr = sourceBuffer.getRawPtr();
      ArrayRCP<const Scalar> vptr = bufView;
      for (size_t j=0; j<numCols; ++j) {
        std::copy(vptr,vptr+myLen,sptr);
        sptr += myLen;
        vptr += myStride;
      }
      target = tmparr();
    }
    // reduce
    reduceAll<int,Scalar> (*comm, REDUCE_SUM, as<int> (numCols*myLen),
                           sourceBuffer.getRawPtr (), target.getRawPtr ());
    if (! packed) {
      // copy tmparr back into view
      const Scalar *sptr = tmparr.getRawPtr();
      ArrayRCP<Scalar> vptr = bufView;
      for (size_t j=0; j<numCols; ++j) {
        std::copy (sptr, sptr+myLen, vptr);
        sptr += myLen;
        vptr += myStride;
      }
    }
    bufView = Teuchos::null;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  replaceLocalValue (LocalOrdinal MyRow,
                     size_t VectorIndex,
                     const Scalar &ScalarValue)
  {
#ifdef HAVE_TPETRA_DEBUG
    const std::string tfecfFuncName("replaceLocalValue()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(MyRow < this->getMap()->getMinLocalIndex() || MyRow > this->getMap()->getMaxLocalIndex(),
      std::runtime_error, ": row index is invalid.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( vectorIndexOutOfRange(VectorIndex),
      std::runtime_error, ": vector index is invalid.");
#endif
    getDataNonConst(VectorIndex)[MyRow] = ScalarValue;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  sumIntoLocalValue (LocalOrdinal MyRow,
                     size_t VectorIndex,
                     const Scalar &ScalarValue)
  {
#ifdef HAVE_TPETRA_DEBUG
    const std::string tfecfFuncName("sumIntoLocalValue");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(MyRow < this->getMap()->getMinLocalIndex() || MyRow > this->getMap()->getMaxLocalIndex(), std::runtime_error,
      ": row index is invalid.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( vectorIndexOutOfRange(VectorIndex), std::runtime_error,
      ": vector index is invalid.");
#endif
    getDataNonConst(VectorIndex)[MyRow] += ScalarValue;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  replaceGlobalValue (GlobalOrdinal GlobalRow, size_t VectorIndex, const Scalar &ScalarValue)
  {
    LocalOrdinal MyRow = this->getMap()->getLocalElement(GlobalRow);
#ifdef HAVE_TPETRA_DEBUG
    const std::string tfecfFuncName("replaceGlobalValue()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(MyRow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error,
        ": row index is not present on this processor.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( vectorIndexOutOfRange(VectorIndex), std::runtime_error,
        ": vector index is invalid.");
#endif
    getDataNonConst(VectorIndex)[MyRow] = ScalarValue;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::sumIntoGlobalValue(GlobalOrdinal GlobalRow, size_t VectorIndex, const Scalar &ScalarValue)
  {
    LocalOrdinal MyRow = this->getMap()->getLocalElement(GlobalRow);
#ifdef HAVE_TEUCHOS_DEBUG
    const std::string tfecfFuncName("sumIntoGlobalValue()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(MyRow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error,
      ": row index is not present on this processor.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( vectorIndexOutOfRange(VectorIndex), std::runtime_error,
      ": vector index is invalid.");
#endif
    getDataNonConst(VectorIndex)[MyRow] += ScalarValue;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  template <class T>
  Teuchos::ArrayRCP<T>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  getSubArrayRCP (Teuchos::ArrayRCP<T> arr,
                  size_t j) const
  {
    const size_t stride = MVT::getStride (lclMV_);
    const size_t myLen = getLocalLength ();
    Teuchos::ArrayRCP<T> ret;
    if (isConstantStride()) {
      ret = arr.persistingView (j*stride, myLen);
    }
    else {
      ret = arr.persistingView (whichVectors_[j]*stride, myLen);
    }
    return ret;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  const Kokkos::MultiVector<Scalar,Node>&
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getLocalMV() const {
    return lclMV_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Kokkos::MultiVector<Scalar,Node>&
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getLocalMVNonConst() {
    return lclMV_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::description() const
  {
    std::ostringstream oss;
    oss << Teuchos::Describable::description();
    oss << "{length = "<<getGlobalLength()
        << ", getNumVectors = "<<getNumVectors()
        << ", isConstantStride() = "<<isConstantStride()
        << "}";
    return oss.str();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  describe (Teuchos::FancyOStream &out,
            const Teuchos::EVerbosityLevel verbLevel) const
  {
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;
    using std::endl;
    using std::setw;

    // Set default verbosity if applicable.
    const Teuchos::EVerbosityLevel vl =
      (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;

    RCP<const Teuchos::Comm<int> > comm = this->getMap()->getComm();
    const int myImageID = comm->getRank();
    const int numImages = comm->getSize();

    if (vl != VERB_NONE) {
      // Don't set the tab level unless we're printing something.
      Teuchos::OSTab tab (out);

      if (myImageID == 0) { // >= VERB_LOW prints description()
        out << this->description() << endl;
      }
      for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
        if (myImageID == imageCtr) {
          if (vl != VERB_LOW) {
            // At verbosity > VERB_LOW, each process prints something.
            out << "Process " << myImageID << ":" << endl;

            Teuchos::OSTab procTab (out);
            // >= VERB_MEDIUM: print the local vector length.
            out << "local length=" << getLocalLength();
            if (vl != VERB_MEDIUM) {
              // >= VERB_HIGH: print isConstantStride() and getStride()
              if (isConstantStride()) {
                out << ", constant stride=" << getStride() << endl;
              }
              else {
                out << ", not constant stride" << endl;
              }
              if (vl == VERB_EXTREME) {
                // VERB_EXTREME: print all the values in the multivector.
                out << "Values:" << endl;
                ArrayRCP<ArrayRCP<const Scalar> > X = this->get2dView();
                for (size_t i = 0; i < getLocalLength(); ++i) {
                  for (size_t j = 0; j < getNumVectors(); ++j) {
                    out << X[j][i];
                    if (j < getNumVectors() - 1) {
                      out << " ";
                    }
                  } // for each column
                  out << endl;
                } // for each row
              } // if vl == VERB_EXTREME
            } // if (vl != VERB_MEDIUM)
            else { // vl == VERB_LOW
              out << endl;
            }
          } // if vl != VERB_LOW
        } // if it is my process' turn to print
      } // for each process in the communicator
    } // if vl != VERB_NONE
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  createViews() const
  {
    TPETRA_EFFICIENCY_WARNING(! cview_.is_null (), std::runtime_error,
      "::createViews(): The const view "
      "has already been created and is therefore not null.  (For"
      "Tpetra developers: cview_.total_count() = " << cview_.total_count ()
      << ".  This "
      "means that MultiVector is either creating a view unnecessarily, or "
      "hanging on to a view beyond its needed scope.  This "
      "probably does not affect correctness, it but does affect total memory "
      "use.  Please report this performance bug to the Tpetra developers.");

    Teuchos::RCP<Node> node = this->getMap ()->getNode ();
    if (cview_.is_null () && getLocalLength () > 0) {
      Teuchos::ArrayRCP<const Scalar> buff = MVT::getValues (lclMV_);
      KOKKOS_NODE_TRACE("MultiVector::createViews()")
      cview_ = node->template viewBuffer<Scalar> (buff.size (), buff);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  createViewsNonConst (Kokkos::ReadWriteOption rwo)
  {

    TPETRA_EFFICIENCY_WARNING(! ncview_.is_null (), std::runtime_error,
      "::createViewsNonConst(): The nonconst view "
      "has already been created and is therefore not null.  (For"
      "Tpetra developers: ncview_.total_count() = " << ncview_.total_count ()
      << ".  This "
      "means that MultiVector is either creating a view unnecessarily, or "
      "hanging on to a view beyond its needed scope.  This "
      "probably does not affect correctness, it but does affect total memory "
      "use.  Please report this performance bug to the Tpetra developers.");

    Teuchos::RCP<Node> node = this->getMap ()->getNode ();
    if (ncview_.is_null () && getLocalLength () > 0) {
      Teuchos::ArrayRCP<Scalar> buff = MVT::getValuesNonConst (lclMV_);
      KOKKOS_NODE_TRACE("MultiVector::createViewsNonConst()")
      ncview_ = node->template viewBufferNonConst<Scalar> (rwo, buff.size (), buff);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::releaseViews () const
  {
    const int constViewCount = cview_.total_count ();
    const int nonconstViewCount = ncview_.total_count ();

    // This prevents unused variable compiler warnings, in case Tpetra
    // efficiency warnings aren't enabled.
    (void) constViewCount;
    (void) nonconstViewCount;

    TPETRA_EFFICIENCY_WARNING(constViewCount > 1 || nonconstViewCount > 1,
      std::runtime_error, "::releaseViews(): Either the const view or the "
      "nonconst view has a reference count greater than 1.  (For Tpetra "
      "developers: cview_.total_count() = " << constViewCount << " and ncview_."
      "total_count() = " << nonconstViewCount << ".  This means that release"
      "Views() won't actually free memory.  This probably does not affect "
      "correctness, it but does affect total memory use.  Please report this "
      "performance bug to the Tpetra developers.");

    cview_ = Teuchos::null;
    ncview_ = Teuchos::null;
  }
} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_MULTIVECTOR_INSTANT(SCALAR,LO,GO,NODE) \
  \
  template class MultiVector< SCALAR , LO , GO , NODE >; \


#endif // TPETRA_MULTIVECTOR_DEF_HPP
