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

#ifndef TPETRA_CRSGRAPH_DECL_HPP
#define TPETRA_CRSGRAPH_DECL_HPP

#include <Teuchos_Describable.hpp>
#include <Teuchos_CompileTimeAssert.hpp>

#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_DefaultKernels.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_RowGraph.hpp"
#include "Tpetra_DistObject.hpp"
#include "Tpetra_Util.hpp"

namespace Tpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration
  template <class S, class LO, class GO, class N, class SpMatOps>
  class CrsMatrix;
#endif

  struct RowInfo {
    size_t localRow;
    size_t allocSize;
    size_t numEntries;
    size_t offset1D;
  };

  enum ELocalGlobal {
    LocalIndices,
    GlobalIndices
  };

  //! \brief A graph accessed by rows and stored sparsely.
  /*!

   \tparam LocalOrdinal A ordinal type for lists of local
     indices. This specifies the \c LocalOrdinal type for Map objects
     used by this graph.

   \tparam GlobalOrdinal A ordinal type for lists of global
     indices. This specifies the \c GlobalOrdinal type for Map objects
     used by this graph.

   \tparam Node A shared-memory node class, fulfilling the \ref
     kokkos_node_api "Kokkos Node API"

   \tparam LocalMatOps A local sparse matrix operations class,
     fulfiling the \ref kokkos_crs_ops "Kokkos CRS Ops API".

   This class allows the construction of sparse graphs with row-access. 
   
   <b>Local vs. Global</b>
   
   Graph entries can be added using either local or global coordinates
   for the indices. The accessors isGloballyIndexed() and
   isLocallyIndexed() indicate whether the indices are currently
   stored as global or local indices. Many of the class methods are
   divided into global and local versions, which differ only in
   whether they accept/return indices in the global or local
   coordinate space. Some of these methods may only be used if the
   graph coordinates are in the appropriate coordinates.  For example,
   getGlobalRowView() returns a View to the indices in global
   coordinates; if the indices are not in global coordinates, then no
   such View can be created.
    
   The global/local distinction does distinguish between operation on
   the global/local graph. Almost all methods operate on the local
   graph, i.e., the rows of the graph associated with the local node,
   per the distribution specified by the row map. Access to non-local
   rows requires performing an explicit communication via the
   import/export capabilities of the CrsGraph object; see
   DistObject. However, the method insertGlobalValues() is an
   exception to this rule, as non-local rows are allowed to be added
   via the local graph. These rows are stored in the local graph and
   communicated to the appropriate node on the next call to
   globalAssemble() or fillComplete() (the latter calls the former).
   */
  template <class LocalOrdinal, 
            class GlobalOrdinal = LocalOrdinal, 
            class Node = Kokkos::DefaultNode::DefaultNodeType,
            class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps >
  class CrsGraph : 
    public RowGraph<LocalOrdinal,GlobalOrdinal,Node>,
    public DistObject<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Node>,
    public Teuchos::ParameterListAcceptorDefaultBase
  {
    template <class S, class LO, class GO, class N, class SpMatOps>
    friend class CrsMatrix;

  public: 
    typedef LocalOrdinal                         local_ordinal_type;
    typedef GlobalOrdinal                        global_ordinal_type;
    typedef Node                                 node_type;
    typedef Map<LocalOrdinal,GlobalOrdinal,Node> map_type;

    //! @name Constructor/Destructor Methods
    //@{ 

    /// \brief Constructor specifying fixed number of entries for each row.
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param maxNumEntriesPerRow [in] Maximum number of graph
    ///   entries per row.  If pftype==DynamicProfile, this is only a
    ///   hint, and you can set this to zero without affecting
    ///   correctness.  If pftype==StaticProfile, this sets the amount
    ///   of storage allocated, and you cannot exceed this number of
    ///   entries in any row.
    ///
    /// \param pftype [in] Whether to allocate storage dynamically
    ///   (DynamicProfile) or statically (StaticProfile).
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsGraph (const RCP<const map_type>& rowMap, 
              size_t maxNumEntriesPerRow, 
              ProfileType pftype = DynamicProfile,
              const RCP<ParameterList>& params = null);

    /// \brief Constructor specifying (possibly different) number of entries in each row.
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param NumEntriesPerRowToAlloc [in] Maximum number of graph
    ///   entries to allocate for each row.  If
    ///   pftype==DynamicProfile, this is only a hint.  If
    ///   pftype==StaticProfile, this sets the amount of storage
    ///   allocated, and you cannot exceed the allocated number of
    ///   entries for any row.
    ///
    /// \param pftype [in] Whether to allocate storage dynamically
    ///   (DynamicProfile) or statically (StaticProfile).
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsGraph (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rowMap, 
              const ArrayRCP<const size_t>& NumEntriesPerRowToAlloc, 
              ProfileType pftype = DynamicProfile,
              const RCP<ParameterList>& params = null);

    /// \brief Constructor specifying column Map and fixed number of entries for each row.
    ///
    /// The column Map will be used to filter any graph indices
    /// inserted using insertLocalIndices() or insertGlobalIndices().
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param colMap [in] Distribution of columns of the graph.
    ///
    /// \param maxNumEntriesPerRow [in] Maximum number of graph
    ///   entries per row.  If pftype==DynamicProfile, this is only a
    ///   hint, and you can set this to zero without affecting
    ///   correctness.  If pftype==StaticProfile, this sets the amount
    ///   of storage allocated, and you cannot exceed this number of
    ///   entries in any row.
    ///
    /// \param pftype [in] Whether to allocate storage dynamically
    ///   (DynamicProfile) or statically (StaticProfile).
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsGraph (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rowMap, 
              const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& colMap, 
              size_t maxNumEntriesPerRow, 
              ProfileType pftype = DynamicProfile,
              const RCP<ParameterList>& params = null);

    /// \brief Constructor specifying column Map and number of entries in each row.
    ///
    /// The column Map will be used to filter any graph indices
    /// inserted using insertLocalIndices() or insertGlobalIndices().
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param colMap [in] Distribution of columns of the graph.
    ///
    /// \param NumEntriesPerRowToAlloc [in] Maximum number of graph
    ///   entries to allocate for each row.  If
    ///   pftype==DynamicProfile, this is only a hint.  If
    ///   pftype==StaticProfile, this sets the amount of storage
    ///   allocated, and you cannot exceed the allocated number of
    ///   entries for any row.
    ///
    /// \param pftype [in] Whether to allocate storage dynamically
    ///   (DynamicProfile) or statically (StaticProfile).
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsGraph (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rowMap, 
              const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& colMap, 
              const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, 
              ProfileType pftype = DynamicProfile,
              const RCP<ParameterList>& params = null);

    //! Destructor.
    virtual ~CrsGraph();

    //@}

    //! @name Implementation of Teuchos::ParameterListAcceptor
    //@{ 

    //! Set the given list of parameters (must be nonnull).
    void setParameterList (const RCP<ParameterList>& params);

    //! Default parameter list suitable for validation.
    RCP<const ParameterList> getValidParameters () const;

    //@}

      //! @name Insertion/Removal Methods
      //@{ 

      //! Insert graph indices, using global IDs.
      /** All index values must be in the global space. 
          \pre \c globalRow exists as an ID in the global row map
          \pre <tt>isLocallyIndexed() == false</tt>
          \pre <tt>isStorageOptimized() == false</tt>

          \post <tt>indicesAreAllocated() == true</tt>
          \post <tt>isGloballyIndexed() == true</tt>

          \note If \c globalRow does not belong to the graph on this node, then it will be communicated to the appropriate node when globalAssemble() is called (which will, at the latest, occur during the next call to fillComplete().) Otherwise, the entries will be inserted in the local graph. 
          \note If the graph row already contains entries at the indices corresponding to values in \c indices, then the redundant indices will be eliminated; this may happen at insertion or during the next call to fillComplete().
        */
      void insertGlobalIndices(GlobalOrdinal globalRow, const ArrayView<const GlobalOrdinal> &indices);

      //! Insert graph indices, using local IDs.
      /**
          \pre \c localRow is a local row belonging to the graph on this node
          \pre <tt>isGloballyIndexed() == false</tt>
          \pre <tt>isStorageOptimized() == false</tt>
          \pre <tt>hasColMap() == true</tt>

          \post <tt>indicesAreAllocated() == true</tt>
          \post <tt>isLocallyIndexed() == true</tt>

          \note If the graph row already contains entries at the indices corresponding to values in \c indices, then the redundant indices will be eliminated; this may happen at insertion or during the next call to fillComplete().
        */
      void insertLocalIndices(LocalOrdinal localRow, const ArrayView<const LocalOrdinal> &indices);

      //! Remove all graph indices from the specified local row.
      /**
          \pre \c localRow is a local row of this graph.
          \pre <tt>isGloballyIndexed() == false</tt>
          \pre <tt>isStorageOptimized() == false</tt>

          \post <tt>getNumEntriesInLocalRow(localRow) == 0</tt>
          \post <tt>indicesAreAllocated() == true</tt>
          \post <tt>isLocallyIndexed() == true</tt>
        */
      void removeLocalIndices(LocalOrdinal localRow);

      //@}

      //! @name Transformational Methods
      /** 
          Each of the methods in this group is a global collective. It is
          necessary to call these mehtods on all nodes participating in the
          communicator associated with this graph.
        */
      //@{ 

      //! \brief Communicate non-local contributions to other nodes.
      void globalAssemble();

      /*! Resume fill operations.
          After calling fillComplete(), resumeFill() must be called before initiating any changes to the graph.

          resumeFill() may be called repeatedly. 

          \post  <tt>isFillActive() == true<tt>
          \post  <tt>isFillComplete() == false<tt>
       */
      void resumeFill(const RCP<ParameterList> &params = null);

      /*! \brief Signal that data entry is complete, specifying domain and range maps. 

          Off-node indices are distributed (via globalAssemble()), indices are sorted, redundant indices are eliminated, and global indices are transformed to local indices.

          \pre  <tt>isFillActive() == true<tt>
          \pre <tt>isFillComplete()() == false<tt>

          \post <tt>isFillActive() == false<tt>
          \post <tt>isFillComplete() == true<tt>
          \post if <tt>os == DoOptimizeStorage<tt>, then <tt>isStorageOptimized() == true</tt>. See isStorageOptimized() for consequences.
       */ 
      void fillComplete(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap, const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap, const RCP<ParameterList> &params = null);

      /*! \brief Signal that data entry is complete. 

          Off-node entries are distributed (via globalAssemble()), repeated entries are summed, and global indices are transformed to local indices.

          \note This method calls fillComplete( getRowMap(), getRowMap(), os ).

          \pre  <tt>isFillActive() == true<tt>
          \pre <tt>isFillComplete()() == false<tt>

          \post <tt>isFillActive() == false<tt>
          \post <tt>isFillComplete() == true<tt>
          \post if <tt>os == DoOptimizeStorage<tt>, then <tt>isStorageOptimized() == true</tt>. See isStorageOptimized() for consequences.
       */
      void fillComplete(const RCP<ParameterList> &params = null);

      //@}

      //! @name Methods implementing RowGraph.
      //@{ 

      //! Returns the communicator.
      const RCP<const Comm<int> > & getComm() const;

      //! Returns the underlying node.
      RCP<Node> getNode() const;

      //! Returns the Map that describes the row distribution in this graph.
      const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getRowMap() const;

      //! \brief Returns the Map that describes the column distribution in this graph.
      const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getColMap() const;

      //! Returns the Map associated with the domain of this graph.
      const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getDomainMap() const;

      //! Returns the Map associated with the domain of this graph.
      const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getRangeMap() const;

      //! Returns the importer associated with this graph.
      RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > getImporter() const;

      //! Returns the exporter associated with this graph.
      RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> > getExporter() const;

      //! Returns the number of global rows in the graph.
      /** Undefined if isFillActive().
        */
      global_size_t getGlobalNumRows() const;

      //! \brief Returns the number of global columns in the graph.
      /** Returns the number of entries in the domain map of the matrix.
          Undefined if isFillActive().
        */
      global_size_t getGlobalNumCols() const;

      //! Returns the number of graph rows owned on the calling node.
      size_t getNodeNumRows() const;

      //! Returns the number of columns connected to the locally owned rows of this graph.
      /** Throws std::runtime_error if <tt>hasColMap() == false</tt>
        */
      size_t getNodeNumCols() const;

      //! Returns the index base for global indices for this graph. 
      GlobalOrdinal getIndexBase() const;

      //! Returns the global number of entries in the graph.
      /** Undefined if isFillActive().
        */
      global_size_t getGlobalNumEntries() const;

      //! Returns the local number of entries in the graph.
      size_t getNodeNumEntries() const;

      //! \brief Returns the current number of entries on this node in the specified global row.
      /*! Returns OrdinalTraits<size_t>::invalid() if the specified global row does not belong to this graph. */
      size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const;

      //! Returns the current number of entries on this node in the specified local row.
      /*! Returns OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this graph. */
      size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const;

      //! \brief Returns the total number of indices allocated for the graph, across all rows on this node.
      /*! This is the allocation available to the user. Actual allocation may be larger, for example, after 
          calling fillComplete(), and thus this does not necessarily reflect the memory consumption of the 
          this graph.  

          This quantity is computed during the actual allocation. Therefore, if <tt>indicesAreAllocated() == false</tt>,
          this method returns <tt>OrdinalTraits<size_t>::invalid()</tt>.
      */
      size_t getNodeAllocationSize() const;

      //! \brief Returns the current number of allocated entries for this node in the specified global row .
      /** Throws exception std::runtime_error if the specified global row does not belong to this node. */
      size_t getNumAllocatedEntriesInGlobalRow(GlobalOrdinal globalRow) const;

      //! Returns the current number of allocated entries on this node in the specified local row.
      /** Throws exception std::runtime_error if the specified local row is not valid for this node. */
      size_t getNumAllocatedEntriesInLocalRow(LocalOrdinal localRow) const;

      //! \brief Returns the number of global diagonal entries, based on global row/column index comparisons. 
      /** Undefined if isFillActive().
        */
      global_size_t getGlobalNumDiags() const;

      //! \brief Returns the number of local diagonal entries, based on global row/column index comparisons. 
      /** Undefined if isFillActive().
        */
      size_t getNodeNumDiags() const;

      //! \brief Returns the maximum number of entries across all rows/columns on all nodes. 
      /** Undefined if isFillActive().
        */
      size_t getGlobalMaxNumRowEntries() const;

      //! \brief Returns the maximum number of entries across all rows/columns on this node. 
      /** Undefined if isFillActive().
        */
      size_t getNodeMaxNumRowEntries() const;

      //! \brief Indicates whether the graph has a well-defined column map. 
      bool hasColMap() const; 

      //! \brief Indicates whether the graph is lower triangular.
      /** Undefined if isFillActive().
        */
      bool isLowerTriangular() const;

      //! \brief Indicates whether the graph is upper triangular.
      /** Undefined if isFillActive().
        */
      bool isUpperTriangular() const;

      //! \brief If graph indices are in the local range, this function returns true. Otherwise, this function returns false. */
      bool isLocallyIndexed() const;

      //! \brief If graph indices are in the global range, this function returns true. Otherwise, this function returns false. */
      bool isGloballyIndexed() const;

      //! Returns \c true if fillComplete() has been called and the graph is in compute mode.
      bool isFillComplete() const;

      //! Returns \c true if resumeFill() has been called and the graph is in edit mode.
      bool isFillActive() const;

      //! Indicates whether the graph indices in all rows are known to be sorted.
      /** A fill-complete graph is always sorted, as is a newly constructed graph. A graph is sorted immediately after 
         calling resumeFill(), but any changes to the graph may result in the sorting status becoming unknown (and therefore, presumed unsorted.)
         */
      bool isSorted() const;

      //! \brief Returns \c true if storage has been optimized.
      /**
        Optimized storage means that the allocation of each row is equal to the
        number of entries. The effect is that a pass through the matrix, i.e.,
        during a mat-vec, requires minimal memory traffic. One limitation of
        optimized storage is that no new indices can be added to the graph.
        */
      bool isStorageOptimized() const;

      //! Returns \c true if the graph was allocated with static data structures.
      ProfileType getProfileType() const;

      //! Extract a list of elements in a specified global row of the graph. Put into pre-allocated storage.
      /*!
        \param LocalRow - (In) Global row number for which indices are desired.
        \param Indices - (Out) Global column indices corresponding to values.
        \param NumIndices - (Out) Number of indices.

         Note: A std::runtime_error exception is thrown indices is not large enough to hold the column indices associated
         with row \c GlobalRow. If \c GlobalRow does not belong to this node, then \c indices is unchanged and \c NumIndices is 
         returned as OrdinalTraits<size_t>::invalid().
       */
      void getGlobalRowCopy(GlobalOrdinal GlobalRow, 
                            const ArrayView<GlobalOrdinal> &Indices, 
                            size_t &NumIndices
                            ) const;

      //! Extract a list of elements in a specified local row of the graph. Put into storage allocated by calling routine.
      /*!
        \param LocalRow - (In) Local row number for which indices are desired.
        \param Indices - (Out) Local column indices corresponding to values.
        \param NumIndices - (Out) Number of indices.

         Note: A std::runtime_error exception is thrown indices is not large enough to hold the column indices associated
         with row \c LocalRow. If \c LocalRow is not valid for this node, then \c indices is unchanged and \c NumIndices is 
         returned as OrdinalTraits<size_t>::invalid().

        \pre <tt>isLocallyIndexed()==true</tt> or <tt>hasColMap() == true</tt>
       */
      void getLocalRowCopy(LocalOrdinal LocalRow, 
                           const ArrayView<LocalOrdinal> &indices, 
                           size_t &NumIndices
                           ) const;

      //! Extract a const, non-persisting view of global indices in a specified row of the graph.
      /*!
        \param GlobalRow - (In) Global row number for which indices are desired.
        \param Indices   - (Out) Global column indices corresponding to values.
        \pre <tt>isLocallyIndexed() == false</tt>
        \post <tt>indices.size() == getNumEntriesInGlobalRow(GlobalRow)</tt>

         Note: If \c GlobalRow does not belong to this node, then \c indices is set to null.
       */
      void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView<const GlobalOrdinal> &Indices) const;

      //! Extract a const, non-persisting view of local indices in a specified row of the graph.
      /*!
        \param LocalRow - (In) Local row number for which indices are desired.
        \param Indices  - (Out) Global column indices corresponding to values.
        \pre <tt>isGloballyIndexed() == false</tt>
        \post <tt>indices.size() == getNumEntriesInLocalRow(LocalRow)</tt>

         Note: If \c LocalRow does not belong to this node, then \c indices is set to null.
       */
      void getLocalRowView(LocalOrdinal LocalRow, ArrayView<const LocalOrdinal> &indices) const;

      //@}

      //! @name Overridden from Teuchos::Describable 
      //@{

      /** \brief Return a simple one-line description of this object. */
      std::string description() const;

      /** \brief Print the object with some verbosity level to an FancyOStream object. */
      void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

      //@}

      //! @name Methods implementing Tpetra::DistObject
      //@{

      bool checkSizes(const DistObject<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Node>& source);

      void copyAndPermute(const DistObject<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Node> & source,
                          size_t numSameIDs,
                          const ArrayView<const LocalOrdinal> &permuteToLIDs,
                          const ArrayView<const LocalOrdinal> &permuteFromLIDs);

      void packAndPrepare(const DistObject<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Node> & source,
                          const ArrayView<const LocalOrdinal> &exportLIDs,
                          Array<GlobalOrdinal> &exports,
                          const ArrayView<size_t> & numPacketsPerLID,
                          size_t& constantNumPackets,
                          Distributor &distor);

      void unpackAndCombine(const ArrayView<const LocalOrdinal> &importLIDs,
                            const ArrayView<const GlobalOrdinal> &imports,
                            const ArrayView<size_t> &numPacketsPerLID,
                            size_t constantNumPackets,
                            Distributor &distor,
                            CombineMode CM);
      //@}

      //! \name Advanced methods, at increased risk of deprecation.
      //@{

      //! Get an ArrayRCP of the row-offsets.
      /*!  The returned buffer exists in host-memory.
       */
      ArrayRCP<const size_t> getNodeRowPtrs() const;

      //! Get an ArrayRCP of the packed column-indices.
      /*!  The returned buffer exists in host-memory.
       */
      ArrayRCP<const LocalOrdinal> getNodePackedIndices() const;

      //@}

      //! \name Deprecated methods; will be removed at some point in the near future.
      //@{

      /** \brief Re-allocate the data into contiguous storage.

          This method is deprecated and will be removed in a future version of Tpetra, as 
          the implementation of storage optimization has been below Tpetra to Kokkos.

          Currently, the implementation simply calls resumeFill() and then fillComplete(OptimizeStorage). As such, it is 
          required to be called by all nodes that participate in the associated communicator.
       */
      TPETRA_DEPRECATED void optimizeStorage();

      //! Deprecated. Get a persisting const view of the elements in a specified global row of the graph.
      TPETRA_DEPRECATED ArrayRCP<const GlobalOrdinal> getGlobalRowView(GlobalOrdinal GlobalRow) const;

      //! Deprecated. Get a persisting const view of the elements in a specified local row of the graph.
      TPETRA_DEPRECATED ArrayRCP<const LocalOrdinal> getLocalRowView(LocalOrdinal LocalRow) const;

      //! \brief Signal that data entry is complete, specifying domain and range maps. 
      TPETRA_DEPRECATED void fillComplete(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap, const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap, OptimizeOption os);

      //! \brief Signal that data entry is complete. 
      TPETRA_DEPRECATED void fillComplete(OptimizeOption os);

      //@}

    private:
      // We forbid copy construction by declaring this method private
      // and not implementing it.
      CrsGraph (const CrsGraph<LocalOrdinal,GlobalOrdinal,Node> &Source);

      // We forbid assignment (operator=) by declaring this method
      // private and not implementing it.
      CrsGraph<LocalOrdinal,GlobalOrdinal,Node>& 
      operator= (const CrsGraph<LocalOrdinal,GlobalOrdinal,Node> &rhs);

    protected:
      typedef typename LocalMatOps::template graph<LocalOrdinal,Node>::graph_type local_graph_type;
      
      // these structs are conveniences, to cut down on the number of
      // arguments to some of the methods below.
      struct SLocalGlobalViews {
        ArrayView<const GlobalOrdinal> ginds;
        ArrayView<const LocalOrdinal>  linds;
      };
      struct SLocalGlobalNCViews {
        ArrayView<GlobalOrdinal>       ginds;
        ArrayView<LocalOrdinal>        linds;
      };
      //
      // Allocation
      //
      bool indicesAreAllocated() const;

      void allocateIndices (ELocalGlobal lg);

      template <class T>
      ArrayRCP<T> allocateValues1D () const;
      template <class T>
      ArrayRCP<ArrayRCP<T> > allocateValues2D () const;

      template <ELocalGlobal lg>
      RowInfo updateAlloc (RowInfo rowinfo, size_t allocSize);
      template <ELocalGlobal lg, class T> 
      RowInfo updateAllocAndValues (RowInfo rowinfo, size_t allocSize, ArrayRCP<T> &rowVals);
      //
      // Local versus global indices
      //
      void computeIndexState();
      void makeColMap();
      void makeIndicesLocal();
      void makeImportExport();
      //
      // insert/suminto/replace
      //
      template<ELocalGlobal lg>
      size_t filterIndices (const SLocalGlobalNCViews &inds) const;

      template<ELocalGlobal lg, class T>
      size_t filterIndicesAndValues (const SLocalGlobalNCViews &inds, const ArrayView<T> &vals) const;

      template<ELocalGlobal lg, ELocalGlobal I> 
      size_t insertIndices (RowInfo rowInfo, const SLocalGlobalViews &newInds);

      template<ELocalGlobal lg, ELocalGlobal I, class IterO, class IterN> 
      void insertIndicesAndValues (RowInfo rowInfo, const SLocalGlobalViews &newInds, IterO rowVals, IterN newVals);

      template<ELocalGlobal lg, class IterO, class IterN, class BinaryFunction> 
      void transformValues (RowInfo rowInfo, const SLocalGlobalViews &inds, IterO rowVals, IterN newVals, BinaryFunction f) const;
      //
      // Sorting and merging
      //
      bool                       isMerged() const;
      void                       setSorted(bool sorted);
      void                       setMerged(bool merged);
      void                       sortAllIndices();
      void                       sortRowIndices(RowInfo rowinfo);
      template <class Scalar> void sortRowIndicesAndValues(RowInfo rowinfo, ArrayView<Scalar> values);
      void                                             mergeAllIndices();
      void                                             mergeRowIndices(RowInfo rowinfo);
      template <class Iter, class BinaryFunction> void mergeRowIndicesAndValues(RowInfo rowinfo, Iter rowValueIter, BinaryFunction f);
      // 
      void setDomainRangeMaps(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap, const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap);
      void staticAssertions() const;
      // global consts
      void clearGlobalConstants();
      void computeGlobalConstants();
      // graph data accessors
      RowInfo                         getRowInfo(size_t myRow) const;
      ArrayView<const LocalOrdinal>   getLocalView(RowInfo rowinfo) const;
      ArrayView<LocalOrdinal>         getLocalViewNonConst(RowInfo rowinfo);
      ArrayView<const GlobalOrdinal>  getGlobalView(RowInfo rowinfo) const;
      ArrayView<GlobalOrdinal>        getGlobalViewNonConst(RowInfo rowinfo);
      size_t                          findLocalIndex(RowInfo rowinfo, LocalOrdinal ind) const;
      size_t                          findGlobalIndex(RowInfo rowinfo, GlobalOrdinal ind) const;
      // local Kokkos objects
      void fillLocalGraph(const RCP<ParameterList> &params);
      const RCP<const local_graph_type> getLocalGraph() const;
      const RCP<local_graph_type> getLocalGraphNonConst();
      // debugging
      void checkInternalState() const;

      // Tpetra support objects
      RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > rowMap_, colMap_, rangeMap_, domainMap_;
      RCP<Import<LocalOrdinal,GlobalOrdinal,Node> > importer_;
      RCP<Export<LocalOrdinal,GlobalOrdinal,Node> > exporter_;

      // local data, stored in a Kokkos::CrsGraph. only initialized after fillComplete()
      RCP<local_graph_type> lclGraph_;

      // Local and Global Counts
      // nodeNumEntries_ and nodeNumAllocated_ are required to be always consistent
      // nodeMaxNumEntries_, nodeNumDiags_ and the global quantities are computed during fillComplete() and only valid when isFillComplete()
      global_size_t globalNumEntries_, globalNumDiags_, globalMaxNumRowEntries_;
      size_t          nodeNumEntries_,   nodeNumDiags_,   nodeMaxNumRowEntries_, nodeNumAllocated_;

      // allocate static or dynamic?
      ProfileType pftype_;
      // requested allocation sizes; we have to preserve these, because we perform late-allocation
      // number of non-zeros to allocate per row; set to null after they are allocated.
      ArrayRCP<const size_t> numAllocPerRow_;
      // number of non-zeros to allocate for all row; either this or numAllocPerRow_ is used, but not both.
      size_t numAllocForAllRows_;

      // graph indices. before allocation, all are null. 
      // after allocation, except during makeIndicesLocal(), one of local or global is null.
      // we will never have 1D and 2D structures being non-null
      // this is host memory
      // 1D == StaticAllocation, 2D == DynamicAllocation
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //
      // 1D/Static structures
      //
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //! lclInds1D_ are the indices for all rows
      ArrayRCP< LocalOrdinal>                     lclInds1D_;
      //! gblInds1D_ are the indices for all rows
      ArrayRCP<GlobalOrdinal>                     gblInds1D_;
      // offset to the beg entries of each row. only used for 1D (Static) allocation.
      // i.e., indices for row R are lclInds1D_[i] for i in [b,e) where b = rowPtrs_[R] and e = rowPtrs_[R+1]
      // only the first numRowEntries_[R] of these are valid
      // both of these are null for 2D (Dynamic) allocations
      // rowPtrs_ has length N+1, while numRowEntries_ has length N
      // TODO: it is possible that making these size_t is overkill; revisit
      ArrayRCP<size_t> rowPtrs_;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //
      // 2D/Dynamic structures. 
      //
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //! <tt>lclInds2D_[r]</tt> are the indices for row \c r. 
      ArrayRCP<ArrayRCP< LocalOrdinal> > lclInds2D_;
      //! <tt>gblInds2D_[r]</tt> are the indices for row \c r. 
      ArrayRCP<ArrayRCP<GlobalOrdinal> > gblInds2D_;

      //! The number valid entries in the row.
      ArrayRCP<size_t>       numRowEntries_;

      // TODO: these might be useful in the future
      // ArrayRCP< typedef ArrayRCP<const GlobalOrdinal>::iterator > gRowPtrs_;
      // ArrayRCP< typedef ArrayRCP<GlobalOrdinal>::iterator > gRowPtrsNC_;
      // ArrayRCP< typedef ArrayRCP<const LocalOrdinal>::iterator > lRowPtrs_;
      // ArrayRCP< typedef ArrayRCP<LocalOrdinal>::iterator > lRowPtrsNC_;

      bool indicesAreAllocated_,
           indicesAreLocal_,
           indicesAreGlobal_,
           fillComplete_, 
           lowerTriangular_,
           upperTriangular_,
           indicesAreSorted_,
           noRedundancies_,
           haveGlobalConstants_;

      // non-local data
      std::map<GlobalOrdinal, std::deque<GlobalOrdinal> > nonlocals_;

  }; // class CrsGraph

} // namespace Tpetra

#endif
