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

#ifndef TPETRA_DISTOBJECT_HPP
#define TPETRA_DISTOBJECT_HPP

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Distributor.hpp"

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>

namespace Tpetra {

  /// \class DistObject 
  /// \brief Base class for distributed Tpetra objects that support data redistribution.
  ///
  /// \c DistObject is a base class for all Tpetra distributed global
  /// objects, including \c CrsMatrix and \c MultiVector.  It provides
  /// the basic mechanisms and interface specifications for importing
  /// and exporting operations using \c Import and \c Export objects.
  ///
  /// \tparam LocalOrdinal The type of local IDs.  Same as \c Map's \c
  ///   LocalOrdinal template parameter.  This should be an integer
  ///   type, preferably signed.
  ///
  /// \tparam GlobalOrdinal The type of global IDs.  Same as Map's \c
  ///   GlobalOrdinal template parameter.  Defaults to the same type
  ///   as LocalOrdinal.  This should also be an integer type,
  ///   preferably signed.
  ///
  /// \tparam Node Same as Map's \c Node template parameter.  Defaults
  ///   to the default Kokkos Node type.
  ///
  /// Most Tpetra users will create subclasses of \c DistObject (like
  /// \c CrsMatrix or \c MultiVector), describe data distribution via
  /// \c Map instances, create data redistribution plans via \c Import
  /// or \c Export, and invoke the \c doImport() or \c doExport()
  /// methods of \c DistObject with an \c Import or \c Export object
  /// to redistribute data.  Thus, the only methods of \c DistObject
  /// of interest to most Tpetra users are \c doImport() and \c
  /// doExport().
  ///
  /// If you want to implement your own \c DistObject subclass, you
  /// should start by implementing the four pure virtual methods: \c
  /// checkSizes(), \c copyAndPermute(), \c packAndPrepare(), and \c
  /// unpackAndCombine().  The implementation of \c doTransfer()
  /// includes documentation that explains how \c DistObject uses
  /// those methods to do data redistribution.
  ///
  /// If you are writing a \c DistObject class that uses Kokkos
  /// compute buffers and aims to work for any Kokkos Node type, you
  /// should also implement the three hooks that create and release
  /// views: \c createViews(), \c createViewsNonConst(), and \c
  /// releaseViews().  The default implementation of these hooks does
  /// nothing.  The documentation of these methods explains different
  /// ways you might choose to implement them.
  ///
  /// Distributed Tpetra objects may be either "distributed global" or
  /// "replicated local."  Distributed global objects are partitioned
  /// across multiple processes in a communicator.  Each process owns
  /// at least one element in the object's Map that is not owned by
  /// another process.  For replicated local objects, each element in
  /// the object's Map is owned redundantly by all processes in the
  /// object's communicator.  Some algorithms use objects that are too
  /// small to be distributed across all processes.  The upper
  /// Hessenberg matrix in a GMRES iterative solve is a good example.
  /// In other cases, such as with block iterative methods, block dot
  /// product functions produce small dense matrices that are required
  /// by all images.  Replicated local objects handle these
  /// situations.
  template <class Packet, 
            class LocalOrdinal = int, 
            class GlobalOrdinal = LocalOrdinal, 
            class Node = Kokkos::DefaultNode::DefaultNodeType>
  class DistObject : virtual public Teuchos::Describable {
  public:
    //! @name Constructors and destructor
    //@{ 

    //! Constructor.
    explicit DistObject (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map);

    //! Copy constructor.
    DistObject (const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>& source);

    //! Destructor (virtual for memory safety of derived classes).
    virtual ~DistObject ();

    //@}
    //! @name Public methods for redistributing data
    //@{ 

    //! Import using an Import object ("forward mode").
    void 
    doImport (const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>& source, 
              const Import<LocalOrdinal,GlobalOrdinal,Node>& importer, 
              CombineMode CM);

    //! Export using an Export object ("forward mode").
    void 
    doExport (const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node> &dest, 
              const Export<LocalOrdinal,GlobalOrdinal,Node>& exporter, 
              CombineMode CM);

    //! Import using an Export object ("reverse mode").
    void 
    doImport (const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>& source,
              const Export<LocalOrdinal,GlobalOrdinal,Node>& exporter, 
              CombineMode CM);

    //! Export using an Import object ("reverse mode").
    void 
    doExport (const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>& dest,
              const Import<LocalOrdinal,GlobalOrdinal,Node>& importer, 
              CombineMode CM);

    //@}
    //! @name Attribute accessor methods
    //@{ 

    /// \brief Whether this is a globally distributed object.
    ///
    /// For a definition of "globally distributed" (and its opposite,
    /// "locally replicated"), see the documentation of Map's \c
    /// isDistributed() method.
    inline bool isDistributed () const;

    //! The Map with which this DistObject was constructed.
    inline const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& 
    getMap() const { return map_; }

    //@}
    //! @name I/O methods
    //@{ 

    /// \brief Print this object to the given output stream.
    ///
    /// We generally assume that all MPI processes can print to the
    /// given stream.
    void print (std::ostream &os) const;

    //@} 


    //@}
    //! @name Implementation of \c Teuchos::Describable
    //@{ 

    /// \brief One-line descriptiion of this object.
    ///
    /// We declare this method virtual so that subclasses of
    /// \c DistObject may override it.
    virtual std::string description () const;


    /// \brief Print a descriptiion of this object to the given output stream.
    ///
    /// We declare this method virtual so that subclasses of
    /// \c DistObject may override it.
    virtual void 
    describe (Teuchos::FancyOStream &out, 
              const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;
    //@} 

  protected:

    /// \enum ReverseOption
    /// \brief Whether the data transfer should be performed in forward or reverse mode.
    ///
    /// "Reverse mode" means calling \c doExport() with an \c Import
    /// object, or calling \c doImport() with an \c Export object.
    /// "Forward mode" means calling \c doExport() with an \c Export
    /// object, or calling \c doImport() with an \c Import object.
    enum ReverseOption {
      DoForward, //*!< Perform the transfer in forward mode.
      DoReverse  //*!< Perform the transfer in reverse mode.
    };

    /// \brief Redistribute data across memory images.
    ///
    /// \param source [in] The source object, to redistribute into
    ///   the destination object, which is \c *this object.
    ///
    /// \param CM [in] The combine mode that describes how to combine
    ///   values that map to the same global ID on the same process.
    ///
    /// \param permuteToLIDs [in] See \c copyAndPermute().
    ///
    /// \param permuteFromLIDs [in] See \c copyAndPermute().
    ///
    /// \param remoteLIDs [in] List of entries (as local IDs) in the
    ///   destination object to receive from other processes.
    ///
    /// \param exportLIDs [in] See \c packAndPrepare().
    ///
    /// \param distor [in/out] The Distributor object that knows how
    ///   to redistribute data.
    ///
    /// \param revOp [in] Whether to do a forward or reverse mode
    ///   redistribution.
    virtual void 
    doTransfer (const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node> &source,
                CombineMode CM,
                size_t numSameIDs,
                const Teuchos::ArrayView<const LocalOrdinal> &permuteToLIDs,
                const Teuchos::ArrayView<const LocalOrdinal> &permuteFromLIDs,
                const Teuchos::ArrayView<const LocalOrdinal> &remoteLIDs,
                const Teuchos::ArrayView<const LocalOrdinal> &exportLIDs,
                Distributor &distor,
                ReverseOption revOp);

    /// \name Methods implemented by subclasses and used by \c doTransfer().
    /// 
    /// The \c doTransfer() method uses the subclass' implementations
    /// of these methods to implement data transfer.  Subclasses of
    /// DistObject must implement these methods.  This is an instance
    /// of the <a
    /// href="http://en.wikipedia.org/wiki/Template_method_pattern">Template
    /// Method Pattern</a>.  ("Template" here doesn't mean "C++
    /// template"; it means "pattern with holes that are filled in by
    /// the subclass' method implementations.")
    //@{ 

    /// \brief Compare the source and target (\e this) objects for compatibility.
    ///
    /// \return True if they are compatible, else false.
    virtual bool 
    checkSizes (const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>& source) = 0;

    /// \brief Perform copies and permutations that are local to this process.
    ///
    /// \param source [in] On entry, the source object, from which we
    ///   are distributing.  We distribute to the destination object,
    ///   which is \c *this object.
    /// \param numSameIDs [in] The umber of elements that
    ///   are the same on the source and destination (this) objects.
    ///   These elements are owned by the same process in both the
    ///   source and destination objects.  No permutation occurs.
    /// \param numPermuteIDs [in] The number of elements that are
    ///   locally permuted between the source and destination objects.
    /// \param permuteToLIDs [in] List of the elements that are
    ///   permuted.  They are listed by their LID in the destination
    ///   object.
    /// \param permuteFromLIDs [in] List of the elements that are
    ///   permuted.  They are listed by their LID in the source
    ///   object.
    virtual void 
    copyAndPermute (const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>& source,
                    size_t numSameIDs,
                    const Teuchos::ArrayView<const LocalOrdinal>& permuteToLIDs,
                    const Teuchos::ArrayView<const LocalOrdinal>& permuteFromLIDs) = 0;

    /// \brief Perform any packing or preparation required for communication.
    ///
    /// \param source [in] Source object for the redistribution.
    ///
    /// \param exportLIDs [in] List of the entries (as local IDs in
    ///   the source object) we will be sending to other images.
    ///
    /// \param exports [out] On exit, the buffer for data to send.
    ///
    /// \param numPacketsPerLID [out] On exit, numPacketsPerLID[i]
    ///   contains the number of packets to be exported for
    ///   exportLIDs[i].  If constantNumPackets is nonzero, you should
    ///   use that instead, and not rely on numPacketsPerLID[i] being
    ///   filled.
    ///
    /// \param constantNumPackets [out] On exit, 0 if numPacketsPerLID
    ///   has variable contents (different size for each LID).  If
    ///   nonzero, then it is expected that num-packets-per-LID is
    ///   constant, and constantNumPackets holds that value.
    ///
    /// \param distor [in] The Distributor object we are using.
    virtual void 
    packAndPrepare (const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>& source,
                    const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
                    Teuchos::Array<Packet>& exports,
                    const Teuchos::ArrayView<size_t>& numPacketsPerLID,
                    size_t& constantNumPackets,
                    Distributor &distor) = 0;

    /// \brief Perform any unpacking and combining after communication.
    ///
    /// \param importLIDs [in] List of the entries (as LIDs in the
    ///   destination object) we received from other images.
    ///
    /// \param imports [in] Buffer containing data we received.
    ///
    /// \param numPacketsPerLID [in] numPacketsPerLID[i] contains the
    ///   number of packets imported for importLIDs[i].
    ///
    /// \param constantNumPackets [in] If nonzero, then
    ///   numPacketsPerLID is constant (same value in all entries) and
    ///   constantNumPackets is that value.  If zero, use
    ///   numPacketsPerLID[i] instead.
    ///
    /// \param distor [in] The Distributor object we are using.
    ///
    /// \param CM [in] The combine mode to use when combining the
    ///   imported entries with existing entries.
    virtual void 
    unpackAndCombine (const Teuchos::ArrayView<const LocalOrdinal> &importLIDs,
                      const Teuchos::ArrayView<const Packet> &imports,
                      const Teuchos::ArrayView<size_t> &numPacketsPerLID,
                      size_t constantNumPackets,
                      Distributor &distor,
                      CombineMode CM) = 0;
    //@} 
    
    /// \brief Hook for creating a const view.
    ///
    /// \c doTransfer() calls this on the source object.  By default,
    /// it does nothing, but the source object can use this as a hint
    /// to fetch data from a compute buffer on an off-CPU device (such
    /// as a GPU) into host memory.
    virtual void createViews () const {}

    /// \brief Hook for creating a nonconst view.
    ///
    /// \c doTransfer() calls this on the destination (\c *this)
    /// object.  By default, it does nothing, but the destination
    /// object can use this as a hint to fetch data from a compute
    /// buffer on an off-CPU device (such as a GPU) into host memory.
    ///
    /// \param rwo [in] Whether to create a write-only or a
    ///   read-and-write view.  For Kokkos Node types where compute
    ///   buffers live in a separate memory space (e.g., in the device
    ///   memory of a discrete accelerator like a GPU), a write-only
    ///   view only requires copying from host memory to the compute
    ///   buffer, whereas a read-and-write view requires copying both
    ///   ways (once to read, from the compute buffer to host memory,
    ///   and once to write, back to the compute buffer).
    virtual void createViewsNonConst (Kokkos::ReadWriteOption rwo) {}

    /// \brief Hook for releasing views.
    ///
    /// \c doTransfer() calls this on both the source and destination
    /// objects, once it no longer needs to access that object's data.
    /// By default, this method does nothing.  Implementations may use
    /// this as a hint to free host memory which is a view of a
    /// compute buffer, once the host memory view is no longer needed.
    /// Some implementations may prefer to mirror compute buffers in
    /// host memory; for these implementations, \c releaseViews() may
    /// do nothing.
    virtual void releaseViews () const {}

    //! The Map over which this object is distributed.
    Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > map_;

  private:
    //! Buffer into which packed data are imported (received from other processes).
    Teuchos::Array<Packet> imports_;

    /// \brief Number of packets to receive for each receive operation.
    /// 
    /// This array is used in \c Distributor::doPosts() (and \c
    /// doReversePosts()) when starting the ireceive operation.  
    ///
    /// This may be ignored in \c doTransfer() if constantNumPackets
    /// is nonzero, indicating a constant number of packets per LID.
    /// (For example, MultiVector sets the constantNumPackets output
    /// argument of \c packAndPrepare() to the number of columns in
    /// the multivector.)
    Teuchos::Array<size_t> numImportPacketsPerLID_;

    //! Buffer from which packed data are exported (sent to other processes).
    Teuchos::Array<Packet> exports_;

    /// \brief Number of packets to send for each send operation.
    ///
    /// This array is used in \c Distributor::doPosts() (and \c
    /// doReversePosts()) for preparing for the send operation.
    ///
    /// This may be ignored in \c doTransfer() if constantNumPackets
    /// is nonzero, indicating a constant number of packets per LID.
    /// (For example, MultiVector sets the constantNumPackets output
    /// argument of \c packAndPrepare() to the number of columns in
    /// the multivector.)
    Teuchos::Array<size_t> numExportPacketsPerLID_;
  }; // class DistObject

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::
  DistObject (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map)
    : map_ (map)
  {}

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::
  DistObject (const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>& source)
    : map_ (source.map_)
  {}

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::~DistObject() 
  {}
  
  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::description () const
  {
    using Teuchos::TypeNameTraits;

    std::ostringstream os;
    os << "Tpetra::DistObject<"
       << TypeNameTraits<Packet>::name ()
       << ", " << TypeNameTraits<LocalOrdinal>::name ()
       << ", " << TypeNameTraits<GlobalOrdinal>::name ()
       << ", " << TypeNameTraits<Node>::name ()
       << ">";
    return os.str ();
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void 
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::
  describe (Teuchos::FancyOStream &out, 
            const Teuchos::EVerbosityLevel verbLevel) const 
  {
    using Teuchos::rcpFromRef;
    using std::endl;

    const Teuchos::EVerbosityLevel vl = (verbLevel == Teuchos::VERB_DEFAULT) ? 
      Teuchos::VERB_LOW : verbLevel;

    if (vl != Teuchos::VERB_NONE) {
      out << this->description () << endl;
      Teuchos::OSTab tab (rcpFromRef (out));
      out << "Export buffer size (in packets): " << exports_.size() << endl
          << "Import buffer size (in packets): " << imports_.size() << endl
          << "Map over which this object is distributed:" << endl;
      map_->describe (out, vl);
    }
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void 
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::
  doImport (const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node> & A, 
            const Import<LocalOrdinal,GlobalOrdinal,Node> & importer, 
            CombineMode CM) 
  {
    TEUCHOS_TEST_FOR_EXCEPTION(*getMap() != *importer.getTargetMap(), 
      std::invalid_argument, "doImport: The target DistObject's Map is not "
      "identical to the Import's target Map.");
    TEUCHOS_TEST_FOR_EXCEPTION(*A.getMap() != *importer.getSourceMap(), 
      std::invalid_argument, "doImport: The source DistObject's Map is not "
      "identical to the Import's source Map.");
    size_t numSameIDs = importer.getNumSameIDs();

    typedef ArrayView<const LocalOrdinal> view_type;
    const view_type exportLIDs      = importer.getExportLIDs();
    const view_type remoteLIDs      = importer.getRemoteLIDs();
    const view_type permuteToLIDs   = importer.getPermuteToLIDs();
    const view_type permuteFromLIDs = importer.getPermuteFromLIDs();
    this->doTransfer (A, CM, numSameIDs, permuteToLIDs, permuteFromLIDs, 
                      remoteLIDs, exportLIDs, importer.getDistributor (), 
                      DoForward);
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void 
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::
  doExport (const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node> & A, 
            const Export<LocalOrdinal,GlobalOrdinal,Node> & exporter, 
            CombineMode CM) 
  {
    TEUCHOS_TEST_FOR_EXCEPTION(   *getMap() != *exporter.getTargetMap(), std::invalid_argument, 
      "doExport: The target DistObject's Map is not identical to the Export's target Map.");
    TEUCHOS_TEST_FOR_EXCEPTION( *A.getMap() != *exporter.getSourceMap(), std::invalid_argument, 
      "doExport: The source DistObject's Map is not identical to the Export's source Map.");
    size_t numSameIDs = exporter.getNumSameIDs();

    typedef ArrayView<const LocalOrdinal> view_type;
    view_type exportLIDs      = exporter.getExportLIDs();
    view_type remoteLIDs      = exporter.getRemoteLIDs();
    view_type permuteToLIDs   = exporter.getPermuteToLIDs();
    view_type permuteFromLIDs = exporter.getPermuteFromLIDs();
    doTransfer (A, CM, numSameIDs, permuteToLIDs, permuteFromLIDs, remoteLIDs, 
                exportLIDs, exporter.getDistributor (), DoForward);
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void 
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::
  doImport (const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node> & A,
            const Export<LocalOrdinal,GlobalOrdinal,Node> & exporter, 
            CombineMode CM) 
  {
    TEUCHOS_TEST_FOR_EXCEPTION(  * getMap() != *exporter.getSourceMap(), std::invalid_argument,
      "doImport (with Export): The target DistObject's Map is not identical to the Export's source Map.");
    TEUCHOS_TEST_FOR_EXCEPTION( *A.getMap() != *exporter.getTargetMap(), std::invalid_argument,
      "doImport (with Export): The source DistObject's Map is not identical to the Export's target Map.");
    size_t numSameIDs = exporter.getNumSameIDs();

    typedef ArrayView<const LocalOrdinal> view_type;
    view_type exportLIDs      = exporter.getRemoteLIDs();
    view_type remoteLIDs      = exporter.getExportLIDs();
    view_type permuteToLIDs   = exporter.getPermuteFromLIDs();
    view_type permuteFromLIDs = exporter.getPermuteToLIDs();
    doTransfer (A, CM, numSameIDs, permuteToLIDs, permuteFromLIDs, remoteLIDs, 
                exportLIDs, exporter.getDistributor (), DoReverse);
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void 
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::
  doExport (const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node> & A,
            const Import<LocalOrdinal,GlobalOrdinal,Node> & importer, 
            CombineMode CM) 
  {
    TEUCHOS_TEST_FOR_EXCEPTION( *getMap() != *importer.getSourceMap(), 
      std::invalid_argument, "doExport (with Import): The target object's Map "
      "is not identical to the Import's source Map.");
    TEUCHOS_TEST_FOR_EXCEPTION( *A.getMap() != *importer.getTargetMap(), 
      std::invalid_argument, "doExport (with Import): The source object's Map "
      "is not identical to the Import's target Map.");
    size_t numSameIDs = importer.getNumSameIDs();

    typedef ArrayView<const LocalOrdinal> view_type;
    view_type exportLIDs      = importer.getRemoteLIDs();
    view_type remoteLIDs      = importer.getExportLIDs();
    view_type permuteToLIDs   = importer.getPermuteFromLIDs();
    view_type permuteFromLIDs = importer.getPermuteToLIDs();
    doTransfer (A, CM, numSameIDs, permuteToLIDs, permuteFromLIDs, remoteLIDs, 
                exportLIDs, importer.getDistributor (), DoReverse);
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool 
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::isDistributed() const {
    return map_->isDistributed ();
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void 
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::
  doTransfer (const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>& source,
              CombineMode CM,
              size_t numSameIDs, 
              const Teuchos::ArrayView<const LocalOrdinal>& permuteToLIDs, 
              const Teuchos::ArrayView<const LocalOrdinal>& permuteFromLIDs,
              const Teuchos::ArrayView<const LocalOrdinal>& remoteLIDs,    
              const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
              Distributor &distor, 
              ReverseOption revOp) 
  {
    using Teuchos::as;

    TEUCHOS_TEST_FOR_EXCEPTION( ! checkSizes(source), std::invalid_argument, 
      "Tpetra::DistObject::doTransfer(): checkSizes() indicates that the "
      "destination object is not a legal target for redistribution from the "
      "source object.  This probably means that they do not have the same "
      "dimensions.  For example, MultiVectors must have the same number of "
      "rows and columns.");
    Kokkos::ReadWriteOption rwo = Kokkos::ReadWrite;
    if (CM == INSERT || CM == REPLACE) {
      const size_t numIDsToWrite = 
        numSameIDs + permuteToLIDs.size() + remoteLIDs.size();
      if (numIDsToWrite == this->getMap()->getNodeNumElements()) {
        // We're overwriting all of our local data in the destination
        // object, so a write-only view suffices.
        //
        // FIXME (mfh 10 Apr 2012) This doesn't make sense for a
        // CrsMatrix with a dynamic graph.  INSERT mode could mean
        // that we're adding new entries to the object, but we don't
        // want to get rid of the old ones.
        rwo = Kokkos::WriteOnly;
      }
    }
    // Tell the source to create a read-only view of its data.  On a
    // discrete accelerator such as a GPU, this brings EVERYTHING from
    // device memory to host memory.
    //
    // FIXME (mfh 23 Mar 2012) By passing in the list of GIDs (or
    // rather, local LIDs to send) and packet counts, createViews()
    // could create a "sparse view" that only brings in the necessary
    // data from device to host memory.
    source.createViews();

    // Tell the target to create a view of its data.  Depending on
    // rwo, this could be a write-only view or a read-and-write view.
    // On a discrete accelerator such as a GPU, a write-only view only
    // requires a transfer from host to device memory.  A
    // read-and-write view requires a two-way transfer.  This has the
    // same problem as createViews(): it transfers EVERYTHING, not
    // just the necessary data.
    //
    // FIXME (mfh 23 Mar 2012) By passing in the list of GIDs (or
    // rather, local LIDs into which to receive) and packet counts,
    // createViewsNonConst() could create a "sparse view" that only
    // transfers the necessary data.
    this->createViewsNonConst(rwo); 

    if (numSameIDs + permuteToLIDs.size()) {
      // There is at least one GID to copy or permute.
      copyAndPermute (source, numSameIDs, permuteToLIDs, permuteFromLIDs);
    }
    size_t constantNumPackets = 0;
    numExportPacketsPerLID_.resize(exportLIDs.size());
    numImportPacketsPerLID_.resize(remoteLIDs.size());

    // Ask the source to pack data.  Also ask it whether there are a
    // constant number of packets per element (constantNumPackets is
    // an output argument).  If there are, constantNumPackets will
    // come back nonzero.  Otherwise, the source will fill the
    // numExportPacketsPerLID_ array.
    packAndPrepare (source, exportLIDs, exports_, numExportPacketsPerLID_(), 
                    constantNumPackets, distor);

    // We don't need the source's data anymore, so it can let go of
    // its views.  On a discrete accelerator, this frees host memory,
    // since device memory has the "master" version of the data.
    source.releaseViews();

    if (constantNumPackets != 0) {
      // There are a constant number of packets per element.  We
      // already know (from the number of "remote" (incoming)
      // elements) how many incoming elements we expect, so we can
      // resize the buffer accordingly.
      const size_t rbufLen = remoteLIDs.size() * constantNumPackets;
      if (as<size_t> (imports_.size()) != rbufLen) {
        imports_.resize (rbufLen);
      }
    }
    if ((isDistributed() && revOp == DoReverse) || 
        (source.isDistributed() && revOp == DoForward)) {
      // call one of the doPostsAndWaits functions
      if (revOp == DoReverse) {
        if (constantNumPackets == 0) { //variable num-packets-per-LID:
          distor.doReversePostsAndWaits (numExportPacketsPerLID_().getConst(), 1,
                                         numImportPacketsPerLID_());
          size_t totalImportPackets = 0;
          for (Array_size_type i = 0; i < numImportPacketsPerLID_.size(); ++i) {
            totalImportPackets += numImportPacketsPerLID_[i];
          }
          imports_.resize(totalImportPackets);
          distor.doReversePostsAndWaits (exports_().getConst(),
                                         numExportPacketsPerLID_(),
                                         imports_(), 
                                         numImportPacketsPerLID_());
        }
        else {
          distor.doReversePostsAndWaits (exports_().getConst(),
                                         constantNumPackets,
                                         imports_());
        }
      }
      else { // revOp == DoForward
        if (constantNumPackets == 0) { //variable num-packets-per-LID:
          distor.doPostsAndWaits (numExportPacketsPerLID_().getConst(), 1,
                                  numImportPacketsPerLID_());
          size_t totalImportPackets = 0;
          for (Array_size_type i = 0; i < numImportPacketsPerLID_.size(); ++i) {
            totalImportPackets += numImportPacketsPerLID_[i];
          }
          imports_.resize(totalImportPackets);
          distor.doPostsAndWaits (exports_().getConst(), 
                                  numExportPacketsPerLID_(),
                                  imports_(), 
                                  numImportPacketsPerLID_());
        }
        else {
          distor.doPostsAndWaits (exports_().getConst(), 
                                  constantNumPackets, 
                                  imports_());
        }
      }
      unpackAndCombine (remoteLIDs, imports_(), numImportPacketsPerLID_(), 
                        constantNumPackets, distor, CM);
    }
    this->releaseViews();
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void 
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::print (std::ostream &os) const
  {
    using Teuchos::FancyOStream;
    using Teuchos::getFancyOStream;
    using Teuchos::RCP;
    using Teuchos::rcpFromRef;
    using std::endl;

    RCP<FancyOStream> out = getFancyOStream (rcpFromRef (os));
    this->describe (*out, Teuchos::VERB_DEFAULT);
  }

} // namespace Tpetra

#endif /* TPETRA_DISTOBJECT_HPP */
