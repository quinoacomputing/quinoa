//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef KOKKOS_PACKEDCRSOPERATORADAPTOR_HPP
#define KOKKOS_PACKEDCRSOPERATORADAPTOR_HPP

namespace Kokkos {

  template <class Scalar, class Ordinal, class Node, class S, class M>
  class PackedCRSOperatorAdaptor {
  public:
    //@{ 
    //! @name Typedefs and structs

    //!
    typedef Scalar  ScalarType;
    //!
    typedef Ordinal OrdinalType;
    //!
    typedef Node    NodeType;

    /** \brief Rebind struct, for specifying type information for a different scalar.
          
        Each subclass should specify one of these.
      */
    template <class OtherScalarn>
    struct rebind {
      typedef PackedCRSOperatorAdaptor<OtherScalarn,Ordinal,Node,S,M> other;
    };

    //@}
    //! @name Constructors/Destructor
    //@{

    //! \brief Constructor accepting and retaining a node object.
    PackedCRSOperatorAdaptor(const RCP<Node> &node, S s, M m);

    //! PackedCRSOperatorAdaptor Destructor
    ~PackedCRSOperatorAdaptor();

    //@}
    //! @name Accessor routines.
    //@{ 
    
    //! Node accessor.
    RCP<Node> getNode() const;

    //@}
    //! @name Initialization of structure
    //@{

    //! Initialize structure of matrix, using CrsGraphHostCompute
    void initializeStructure(const CrsGraph<Ordinal,Node,PackedCRSOperatorAdaptor<void,Ordinal,Node,S,M> > &graph);

    //! Initialize values of matrix, using CrsMatrixHostCompute
    void initializeValues(const CrsGraph<Scalar,Ordinal,Node,PackedCRSOperatorAdaptor<void,Ordinal,Node,S,M> > &matrix);

    //! Clear all matrix structure and values.
    void clear();

    //@}
    //! @name Computational methods
    //@{

    //! Applies the matrix to a MultiVector, overwriting Y.
    template <class DomainScalar, class RangeScalar>
    void multiply(Teuchos::ETransp trans, RangeScalar alpha, const MultiVector<DomainScalar,Node> &X, MultiVector<RangeScalar,Node> &Y) const;

    //! Applies the matrix to a MultiVector, accumulating into Y.
    template <class DomainScalar, class RangeScalar>
    void multiply(Teuchos::ETransp trans, 
                  RangeScalar alpha, const MultiVector<DomainScalar,Node> &X, RangeScalar beta, MultiVector<RangeScalar,Node> &Y) const;

    //! Solves the matrix for a given set of right-hand-sides.
    template <class DomainScalar, class RangeScalar>
    void solve(Teuchos::ETransp trans, Teuchos::EUplo uplo, Teuchos::EDiag diag, 
               const MultiVector<DomainScalar,Node> &Y, MultiVector<RangeScalar,Node> &X) const;

    //@}

  protected:
    //! Copy constructor (protected and unimplemented)
    PackedCRSOperatorAdaptor(const PackedCRSOperatorAdaptor& source);

    RCP<Node> node_;

    // we do this one of two ways: 
    // 1D/packed: arrays of offsets, array of ordinals, array of values.
    // may or may not use these, will keep them
    ArrayRCP<const Ordinal> inds1D_;
    ArrayRCP<const size_t>  begs1D_, ends1D_;
    // 2D: array of arrays
    // don't use these, but must hold on to them in order to support sequential calls to initializeValues() 
    ArrayRCP<const ArrayRCP<Ordinal> > inds2D_;
    ArrayRCP<const size_t>         numEntries_;

    size_t numRows_;
    bool indsInit_, valsInit_, isEmpty_, isOpt_;

    M mult_;
    S solv_;

    static void packFrom2D(const ArrayRCP<const size_t>  &begs, 
                           const ArrayRCP<const size_t>  &ends,
                           const ArrayRCP<const Ordinal> &inds,
                           const ArrayRCP<const Scalar > &vals,
                                 ArrayRCP<size_t>  &rowOffsets,
                                 ArrayRCP<Ordinal> &packedInds,
                                 ArrayRCP<Scalar>  &packedVals);

    static void packFrom2D(const ArrayRCP<const size_t>  &begs, 
                           const ArrayRCP<const size_t>  &ends,
                           const ArrayRCP<ArrayRCP<Ordinal> > &inds,
                           const ArrayRCP<ArrayRCP<Scalar > > &vals,
                                 ArrayRCP<size_t>  &rowOffsets,
                                 ArrayRCP<Ordinal> &packedInds,
                                 ArrayRCP<Scalar>  &packedVals);
  };


  template<class Scalar, class Ordinal, class Node>
  PackedCRSOperatorAdaptor<Scalar,Ordinal,Node>::PackedCRSOperatorAdaptor(const RCP<Node> &node, M m, S s)
  : node_(node) 
  , mult_(m)
  , solv_(s)
  {
    clear();
  }

  template<class Scalar, class Ordinal, class Node>
  PackedCRSOperatorAdaptor<Scalar,Ordinal,Node>::~PackedCRSOperatorAdaptor() {
  }

  template <class Scalar, class Ordinal, class Node>
  RCP<Node> PackedCRSOperatorAdaptor<Scalar,Ordinal,Node>::getNode() const {
    return node_; 
  }

  template <class Scalar, class Ordinal, class Node>
  void PackedCRSOperatorAdaptor<Scalar,Ordinal,Node>::clear() {
    begs1D_     = null;
    ends1D_     = null;
    inds1D_     = null;
    inds2D_     = null;
    numEntries_ = null;
    numRows_  = 0;
    indsInit_ = false;
    valsInit_ = false;
    isEmpty_  = false;
    isOpt_    = false;
  }

  template <class Scalar, class Ordinal, class Node>
  void PackedCRSOperatorAdaptor<Scalar,Ordinal,Node>::initializeStructure(const CrsGraphHostCompute<Ordinal,Node,PackedCRSOperatorAdaptor<void,Ordinal,Node> > &graph) 
  {
    using Teuchos::arcp;
    TEST_FOR_EXCEPTION(indsInit_ == true || valsInit_ == true, std::runtime_error, 
        Teuchos::typeName(*this) << "::initializeStructure(): structure already initialized.");
    numRows_ = graph.getNumRows();
    if (graph.isEmpty() || numRows_ == 0) {
      isEmpty_ = true;
      isOpt_   = false;
    }
    else if (graph.is1DStructure()) {
      isEmpty_ = false;
      isOpt_   = graph.isOptimized();
      ArrayRCP<Ordinal> inds;
      ArrayRCP<size_t> begs, ends;
      const_cast<CrsGraphHostCompute<Ordinal,Node,PackedCRSOperatorAdaptor<void,Ordinal,Node> > &>(graph).get1DStructure( inds, begs, ends );
      inds1D_ = inds;
      begs1D_ = begs;
      ends1D_ = ends;
    }
    else {
      isEmpty_ = false;
      isOpt_   = false;
      ArrayRCP<ArrayRCP<Ordinal> > inds;
      ArrayRCP<size_t>  sizes;
      const_cast<CrsGraphHostCompute<Ordinal,Node,PackedCRSOperatorAdaptor<void,Ordinal,Node> > &>(graph).get2DStructure(inds,sizes);
      inds2D_     = inds;
      numEntries_ = sizes;
    }
    indsInit_ = true;
  }


  template <class Scalar, class Ordinal, class Node>
  void PackedCRSOperatorAdaptor<Scalar,Ordinal,Node>::initializeValues(const CrsMatrixHostCompute<Scalar,Ordinal,Node,PackedCRSOperatorAdaptor<void,Ordinal,Node> > &matrix) {
    using Teuchos::arcp;
    TEST_FOR_EXCEPTION(indsInit_ == false, std::runtime_error, 
        Teuchos::typeName(*this) << "::initializeValues(): must initialize values after graph.");
    TEST_FOR_EXCEPTION(numRows_ != matrix.getNumRows() || isEmpty_ != matrix.isEmpty() || 
                       (inds2D_ != null && matrix.is1DStructure()) || (inds1D_ != null && matrix.is2DStructure()),
                       std::runtime_error, Teuchos::typeName(*this) << "::initializeValues(): matrix not compatible with previously supplied graph.");
    if (!isEmpty_) {
      if (matrix.is1DStructure()) {
        ArrayRCP<Scalar> vals;
        const_cast<CrsMatrixHostCompute<Scalar,Ordinal,Node,PackedCRSOperatorAdaptor<void,Ordinal,Node> > &>(matrix).get1DValues( vals );
        vals1D_ = vals;
      }
      else {
        {
          ArrayRCP<ArrayRCP<Scalar> > vals;
          const_cast<CrsMatrixHostCompute<Scalar,Ordinal,Node,PackedCRSOperatorAdaptor<void,Ordinal,Node> > &>(matrix).get2DValues(vals);
          vals2D_ = vals;
        }
      }
    }
    valsInit_ = true;
  }


  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void PackedCRSOperatorAdaptor<Scalar,Ordinal,Node>::solve(
                      Teuchos::ETransp trans, Teuchos::EUplo uplo, Teuchos::EDiag diag, 
                      const MultiVector<DomainScalar,Node> &Y,
                            MultiVector<RangeScalar,Node> &X) const 
  {
    TEST_FOR_EXCEPTION(indsInit_ == false || valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::solve(): this solve was not fully initialized.");
    TEST_FOR_EXCEPTION(X.getNumCols() != Y.getNumCols(), std::runtime_error,
        Teuchos::typeName(*this) << "::solve(): Left hand side and right hand side multivectors have differing numbers of vectors.");
    TEST_FOR_EXCEPTION(X.getNumRows() < numRows_, std::runtime_error,
        Teuchos::typeName(*this) << "::solve(): Left-hand-side multivector does not have enough rows. Likely cause is that the column map was not provided to the Tpetra::CrsMatrix in the case of an implicit unit diagonal.");
    ReadyBufferHelper<Node> rbh(node_);
    if (numRows_ == 0) {
      // null op
    }
    else if (isEmpty_) {
      TEST_FOR_EXCEPTION(diag != Teuchos::UNIT_DIAG, std::runtime_error,
          Teuchos::typeName(*this) << "::solve(): solve of empty matrix only valid for an implicit unit diagonal.");
      // solve I * X = Y for X = Y
      DefaultArithmetic<MultiVector<RangeScalar,Node> >::Assign(X,Y);
    }
    else {
        TEST_FOR_EXCEPT(true);
        rbh.begin();
        // RangeScalar  *x = rbh.template addNonConstBuffer<DomainScalar>(X.getValuesNonConst());
        // DomainScalar *y = rbh.template addConstBuffer<RangeScalar>(Y.getValues());
        rbh.end();
        // wdp.numRows = numRows_;
        // wdp.xstride = X.getStride();
        // wdp.ystride = Y.getStride();
        // trans == Teuchos::NO_TRANS 
        // diag == Teuchos::UNIT_DIAG 
        // uplo == Teuchos::UPPER_TRI 
      }
    }
    return;
  }


  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void PackedCRSOperatorAdaptor<Scalar,Ordinal,Node>::multiply(
                                Teuchos::ETransp trans, 
                                RangeScalar alpha,
                                const MultiVector<DomainScalar,Node> &X, 
                                      MultiVector<RangeScalar ,Node> &Y) const 
  {
    typedef DefaultSparseMultiplyOp1<Scalar,Ordinal,DomainScalar,RangeScalar, 1>  Op1D;
    typedef DefaultSparseMultiplyOp2<Scalar,Ordinal,DomainScalar,RangeScalar, 1>  Op2D;
    typedef DefaultSparseTransposeMultiplyOp1<Scalar,Ordinal,DomainScalar,RangeScalar, 1> TOp1D;
    typedef DefaultSparseTransposeMultiplyOp2<Scalar,Ordinal,DomainScalar,RangeScalar, 1> TOp2D;
    TEST_FOR_EXCEPTION(indsInit_ == false || valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::multiply(): operation not fully initialized.");
    TEST_FOR_EXCEPT(X.getNumCols() != Y.getNumCols());
    ReadyBufferHelper<Node> rbh(node_);
    if (isEmpty_ == true) {
      // Y <= 0 * X
      //   <= 0
      DefaultArithmetic<MultiVector<RangeScalar,Node> >::Init(Y,Teuchos::ScalarTraits<RangeScalar>::zero());
    }
    else if (begs1D_ != null) {
      if (trans == Teuchos::NO_TRANS) {
        Op1D wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = Teuchos::ScalarTraits<RangeScalar>::zero(); // not used
        wdp.numRows = numRows_;
        wdp.begs    = rbh.template addConstBuffer<size_t>(begs1D_);
        wdp.ends    = rbh.template addConstBuffer<size_t>(ends1D_);
        wdp.inds    = rbh.template addConstBuffer<Ordinal>(inds1D_);
        wdp.vals    = rbh.template addConstBuffer<Scalar>(vals1D_);
        wdp.x       = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.y       = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<Op1D>(0,numRows_*numRHS,wdp);
      }
      else {
        TOp1D wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = Teuchos::ScalarTraits<RangeScalar>::zero(); // not used
        wdp.numRows = numRows_;
        wdp.numCols = Y.getNumRows();
        wdp.begs    = rbh.template addConstBuffer<size_t>(begs1D_);
        wdp.ends    = rbh.template addConstBuffer<size_t>(ends1D_);
        wdp.inds    = rbh.template addConstBuffer<Ordinal>(inds1D_);
        wdp.vals    = rbh.template addConstBuffer<Scalar>(vals1D_);
        wdp.x       = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.y       = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<TOp1D>(0,numRHS,wdp);
      }
    }
    else {
      if (trans == Teuchos::NO_TRANS) {
        Op2D wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = Teuchos::ScalarTraits<RangeScalar>::zero(); // not used
        wdp.numRows = numRows_;
        wdp.numEntries = rbh.template addConstBuffer<size_t>(numEntries_);
        wdp.inds_beg   = rbh.template addConstBuffer<const Ordinal *>(indPtrs_);
        wdp.vals_beg   = rbh.template addConstBuffer<const  Scalar *>(valPtrs_);
        wdp.x          = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.y          = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        rbh.end();
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<Op2D>(0,numRows_*numRHS,wdp);
      }
      else {
        TOp2D wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = Teuchos::ScalarTraits<RangeScalar>::zero(); // not used
        wdp.numRows = numRows_;
        wdp.numCols = Y.getNumRows();
        wdp.numEntries = rbh.template addConstBuffer<size_t>(numEntries_);
        wdp.inds_beg   = rbh.template addConstBuffer<const Ordinal *>(indPtrs_);
        wdp.vals_beg   = rbh.template addConstBuffer<const  Scalar *>(valPtrs_);
        wdp.x          = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.y          = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<TOp2D>(0,numRHS,wdp);
      }
    }
    return;
  }


  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void PackedCRSOperatorAdaptor<Scalar,Ordinal,Node>::multiply(
                                Teuchos::ETransp trans, 
                                RangeScalar alpha, const MultiVector<DomainScalar,Node> &X, 
                                RangeScalar beta, MultiVector<RangeScalar,Node> &Y) const 
  {
    // the 0 parameter means that beta is considered, and the output multivector enjoys accumulate semantics
    typedef DefaultSparseMultiplyOp1<Scalar,Ordinal,DomainScalar,RangeScalar, 0>  Op1D;
    typedef DefaultSparseMultiplyOp2<Scalar,Ordinal,DomainScalar,RangeScalar, 0>  Op2D;
    typedef DefaultSparseTransposeMultiplyOp1<Scalar,Ordinal,DomainScalar,RangeScalar, 0> TOp1D;
    typedef DefaultSparseTransposeMultiplyOp2<Scalar,Ordinal,DomainScalar,RangeScalar, 0> TOp2D;
    TEST_FOR_EXCEPTION(indsInit_ == false || valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::multiply(): operation not fully initialized.");
    TEST_FOR_EXCEPT(X.getNumCols() != Y.getNumCols());
    ReadyBufferHelper<Node> rbh(node_);
    if (isEmpty_ == true) {
      // Y <= alpha * 0 * X + beta * Y 
      //   <= beta * Y
      // NOTE: this neglects NaNs in X, which don't satisfy 0*NaN == 0
      //       however, X and Y may be of different size, and we need entries to determine how to mix those potential NaNs in X into Y
      //       therefore, the best we can do is scale Y to zero. Setting Y to zero would destroy NaNs in Y, which violates the semantics of the call.
      DefaultArithmetic<MultiVector<RangeScalar,Node> >::Scale(Y,beta);
    }
    else if (begs1D_ != null) {
      if (trans == Teuchos::NO_TRANS) {
        Op1D wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = beta;
        wdp.numRows = numRows_;
        wdp.begs    = rbh.template addConstBuffer<size_t>(begs1D_);
        wdp.ends    = rbh.template addConstBuffer<size_t>(ends1D_);
        wdp.inds    = rbh.template addConstBuffer<Ordinal>(inds1D_);
        wdp.vals    = rbh.template addConstBuffer< Scalar>(vals1D_);
        wdp.x       = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.y       = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<Op1D>(0,numRows_*numRHS,wdp);
      }
      else {
        TOp1D wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = beta;
        wdp.numRows = numRows_;
        wdp.numCols = Y.getNumRows();
        wdp.begs    = rbh.template addConstBuffer<size_t>(begs1D_);
        wdp.ends    = rbh.template addConstBuffer<size_t>(ends1D_);
        wdp.inds    = rbh.template addConstBuffer<Ordinal>(inds1D_);
        wdp.vals    = rbh.template addConstBuffer< Scalar>(vals1D_);
        wdp.x       = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.y       = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<TOp1D>(0,numRHS,wdp);
      }
    }
    else {
      if (trans == Teuchos::NO_TRANS) {
        Op2D wdp;
        rbh.begin();
        wdp.numRows = numRows_;
        wdp.numEntries = rbh.template addConstBuffer<size_t>(numEntries_);
        wdp.inds_beg   = rbh.template addConstBuffer<const Ordinal *>(indPtrs_);
        wdp.vals_beg   = rbh.template addConstBuffer<const  Scalar *>(valPtrs_);
        wdp.x          = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.y          = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        rbh.end();
        wdp.alpha   = alpha;
        wdp.beta    = beta;
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<Op2D>(0,numRows_*numRHS,wdp);
      }
      else {
        TOp2D wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = beta;
        wdp.numRows = numRows_;
        wdp.numCols = Y.getNumRows();
        wdp.numEntries = rbh.template addConstBuffer<size_t>(numEntries_);
        wdp.inds_beg   = rbh.template addConstBuffer<const Ordinal *>(indPtrs_);
        wdp.vals_beg   = rbh.template addConstBuffer<const  Scalar *>(valPtrs_);
        wdp.x          = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.y          = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<TOp2D>(0,numRHS,wdp);
      }
    }
    return;
  }


  /** \brief Default implementation of sparse matrix-vector multiplication and solve routines, for device-based nodes.
      \ingroup kokkos_crs_ops
    */
  template <class Scalar, class Ordinal, class Node>
  class DefaultDeviceSparseOps {
  public:
    typedef Scalar  ScalarType;
    typedef Ordinal OrdinalType;
    typedef Node    NodeType;

    template <class S2>
    struct rebind {
      typedef DefaultDeviceSparseOps<S2,Ordinal,Node> other;
    };

    //! @name Constructors/Destructor

    //@{

    //! DefaultDeviceSparseOps constuctor with variable number of indices per row.
    DefaultDeviceSparseOps(const RCP<Node> &node);

    //! DefaultDeviceSparseOps Destructor
    ~DefaultDeviceSparseOps();

    //@}

    //! @name Accessor routines.
    //@{ 
    
    //! Node accessor.
    RCP<Node> getNode() const;

    //@}

    //! @name Initialization of structure

    //@{

    //! Initialize structure of matrix, using CrsGraphDeviceCompute
    void initializeStructure(const CrsGraphDeviceCompute<Ordinal,Node,DefaultDeviceSparseOps<void,Ordinal,Node> > &graph);

    //! Initialize values of matrix, using CrsMatrixDeviceCompute
    void initializeValues(const CrsMatrixDeviceCompute<Scalar,Ordinal,Node,DefaultDeviceSparseOps<void,Ordinal,Node> > &matrix);

    //! Clear all matrix structure and values.
    void clear();

    //@}

    //! @name Computational methods

    //@{

    //! Applies the matrix to a MultiVector, overwriting Y.
    template <class DomainScalar, class RangeScalar>
    void multiply(Teuchos::ETransp trans, RangeScalar alpha, const MultiVector<DomainScalar,Node> &X, MultiVector<RangeScalar,Node> &Y) const;

    //! Applies the matrix to a MultiVector, accumulating into Y.
    template <class DomainScalar, class RangeScalar>
    void multiply(Teuchos::ETransp trans, 
                  RangeScalar alpha, const MultiVector<DomainScalar,Node> &X, RangeScalar beta, MultiVector<RangeScalar,Node> &Y) const;

    //! Solves the matrix for a given set of right-hand-sides.
    template <class DomainScalar, class RangeScalar>
    void solve(Teuchos::ETransp trans, Teuchos::EUplo uplo, Teuchos::EDiag diag, 
               const MultiVector<DomainScalar,Node> &Y, MultiVector<RangeScalar,Node> &X) const;

    //@}

  protected:
    //! Copy constructor (protected and unimplemented)
    DefaultDeviceSparseOps(const DefaultDeviceSparseOps& source);

    RCP<Node> node_;

    // we do this one of two ways: 
    // 1D: array of offsets, pointer for ordinals, pointer for values.
    ArrayRCP<const size_t>  pbuf_offsets1D_;
    ArrayRCP<const Ordinal> pbuf_inds1D_;
    ArrayRCP<const Scalar>  pbuf_vals1D_;

    size_t numRows_;
    bool indsInit_, valsInit_, isEmpty_;
  };

  template<class Scalar, class Ordinal, class Node>
  DefaultDeviceSparseOps<Scalar,Ordinal,Node>::DefaultDeviceSparseOps(const RCP<Node> &node)
  : node_(node) 
  {
    clear();
  }

  template<class Scalar, class Ordinal, class Node>
  DefaultDeviceSparseOps<Scalar,Ordinal,Node>::~DefaultDeviceSparseOps() {
  }

  template <class Scalar, class Ordinal, class Node>
  RCP<Node> DefaultDeviceSparseOps<Scalar,Ordinal,Node>::getNode() const { 
    return node_; 
  }

  template <class Scalar, class Ordinal, class Node>
  void DefaultDeviceSparseOps<Scalar,Ordinal,Node>::clear() {
    pbuf_offsets1D_  = null;
    pbuf_inds1D_     = null;
    pbuf_vals1D_     = null;
    numRows_  = 0;
    indsInit_ = false;
    valsInit_ = false;
    isEmpty_  = false;
  }

  template <class Scalar, class Ordinal, class Node>
  void DefaultDeviceSparseOps<Scalar,Ordinal,Node>::initializeStructure(const CrsGraphDeviceCompute<Ordinal,Node,DefaultDeviceSparseOps<void,Ordinal,Node> > &graph) {
    TEST_FOR_EXCEPTION(indsInit_ == true || valsInit_ == true, std::runtime_error, 
        Teuchos::typeName(*this) << "::initializeStructure(): structure already initialized.");
    numRows_ = graph.getNumRows();
    if (graph.isEmpty() || numRows_ == 0) {
      isEmpty_ = true;
    }
    else {
      ArrayRCP<Ordinal> inds;
      ArrayRCP<size_t > offs;
      const_cast<CrsGraphDeviceCompute<Ordinal,Node,DefaultDeviceSparseOps<void,Ordinal,Node> > &>(graph).getDeviceBuffers(inds, offs);
      pbuf_inds1D_    = inds;
      pbuf_offsets1D_ = offs;
    }
    indsInit_ = true;
  }

  template <class Scalar, class Ordinal, class Node>
  void DefaultDeviceSparseOps<Scalar,Ordinal,Node>::initializeValues(const CrsMatrixDeviceCompute<Scalar,Ordinal,Node,DefaultDeviceSparseOps<void,Ordinal,Node> > &matrix) {
    TEST_FOR_EXCEPTION(indsInit_ == false, std::runtime_error, 
        Teuchos::typeName(*this) << "::initializeValues(): must initialize values after graph.");
    TEST_FOR_EXCEPTION(numRows_ != matrix.getNumRows() || isEmpty_ != matrix.isEmpty(), std::runtime_error,
        Teuchos::typeName(*this) << "::initializeValues(): matrix not compatible with previously supplied graph.");
    if (!isEmpty_) {
      ArrayRCP<Scalar> vals;
      const_cast<CrsMatrixDeviceCompute<Scalar,Ordinal,Node,DefaultDeviceSparseOps<void,Ordinal,Node> > &>(matrix).getDeviceBuffer(vals);
      pbuf_vals1D_ = vals;
    }
    valsInit_ = true;
  }


  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void DefaultDeviceSparseOps<Scalar,Ordinal,Node>::solve(
                      Teuchos::ETransp trans, Teuchos::EUplo uplo, Teuchos::EDiag diag, 
                      const MultiVector<DomainScalar,Node> &Y,
                            MultiVector<RangeScalar,Node> &X) const {
    typedef DefaultSparseSolveOp1<Scalar,Ordinal,DomainScalar,RangeScalar>  Op1D;
    typedef DefaultSparseTransposeSolveOp1<Scalar,Ordinal,DomainScalar,RangeScalar>  TOp1D;
    TEST_FOR_EXCEPTION(indsInit_ == false || valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::solve(): this solve was not fully initialized.");
    TEST_FOR_EXCEPTION(X.getNumCols() != Y.getNumCols(), std::runtime_error,
        Teuchos::typeName(*this) << "::solve(): Left hand side and right hand side multivectors have differing numbers of vectors.");
    TEST_FOR_EXCEPTION(X.getNumRows() < numRows_, std::runtime_error,
        Teuchos::typeName(*this) << "::solve(): Left-hand-side multivector does not have enough rows. Likely cause is that the column map was not provided to the Tpetra::CrsMatrix in the case of an implicit unit diagonal.");

    ReadyBufferHelper<Node> rbh(node_);
    if (numRows_ == 0) {
      // null op
    }
    else if (isEmpty_) {
      TEST_FOR_EXCEPTION(diag != Teuchos::UNIT_DIAG, std::runtime_error,
          Teuchos::typeName(*this) << "::solve(): solve of empty matrix only valid for an implicit unit diagonal.");
      // solve I * X = Y for X = Y
      DefaultArithmetic<MultiVector<RangeScalar,Node> >::Assign(X,Y);
    }
    else {
      if (trans == Teuchos::NO_TRANS) {
        Op1D wdp;
        rbh.begin();
        wdp.begs    = rbh.template addConstBuffer<size_t>(pbuf_offsets1D_);
        wdp.ends    = wdp.begs+1;
        wdp.inds    = rbh.template addConstBuffer<Ordinal>(pbuf_inds1D_);
        wdp.vals    = rbh.template addConstBuffer< Scalar>(pbuf_vals1D_);
        wdp.x       = rbh.template addNonConstBuffer<DomainScalar>(X.getValuesNonConst());
        wdp.y       = rbh.template addConstBuffer<RangeScalar>(Y.getValues());
        rbh.end();
        wdp.numRows = numRows_;
        wdp.unitDiag = (diag == Teuchos::UNIT_DIAG ? true : false);
        wdp.upper    = (uplo == Teuchos::UPPER_TRI ? true : false);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<Op1D>(0,numRHS,wdp);
      }
      else {
        TOp1D wdp;
        rbh.begin();
        wdp.begs    = rbh.template addConstBuffer<size_t>(pbuf_offsets1D_);
        wdp.ends    = wdp.begs+1;
        wdp.inds    = rbh.template addConstBuffer<Ordinal>(pbuf_inds1D_);
        wdp.vals    = rbh.template addConstBuffer< Scalar>(pbuf_vals1D_);
        wdp.x       = rbh.template addNonConstBuffer<DomainScalar>(X.getValuesNonConst());
        wdp.y       = rbh.template addConstBuffer<RangeScalar>(Y.getValues());
        rbh.end();
        wdp.numRows = numRows_;
        wdp.unitDiag = (diag == Teuchos::UNIT_DIAG ? true : false);
        wdp.upper    = (uplo == Teuchos::UPPER_TRI ? true : false);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<TOp1D>(0,numRHS,wdp);
      }
    }
    return;
  }


  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void DefaultDeviceSparseOps<Scalar,Ordinal,Node>::multiply(
                                Teuchos::ETransp trans, 
                                RangeScalar alpha,
                                const MultiVector<DomainScalar,Node> &X, 
                                MultiVector<RangeScalar,Node> &Y) const {
    typedef DefaultSparseMultiplyOp1<Scalar,Ordinal,DomainScalar,RangeScalar, 1>  Op1D;
    typedef DefaultSparseTransposeMultiplyOp1<Scalar,Ordinal,DomainScalar,RangeScalar, 1> TOp1D;
    TEST_FOR_EXCEPTION(indsInit_ == false || valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::multiply(): operation not fully initialized.");
    TEST_FOR_EXCEPT(X.getNumCols() != Y.getNumCols());
    ReadyBufferHelper<Node> rbh(node_);
    if (isEmpty_ == true) {
      // Y <= 0 * X
      //   <= 0
      DefaultArithmetic<MultiVector<RangeScalar,Node> >::Init(Y,Teuchos::ScalarTraits<RangeScalar>::zero());
    }
    else {
      if (trans == Teuchos::NO_TRANS) {
        Op1D wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = Teuchos::ScalarTraits<RangeScalar>::zero(); // not used
        wdp.numRows = numRows_;
        wdp.begs    = rbh.template addConstBuffer<size_t>(pbuf_offsets1D_);
        wdp.ends    = wdp.begs+1;
        wdp.inds    = rbh.template addConstBuffer<Ordinal>(pbuf_inds1D_);
        wdp.vals    = rbh.template addConstBuffer<Scalar>(pbuf_vals1D_);
        wdp.x       = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.y       = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<Op1D>(0,numRows_*numRHS,wdp);
      }
      else {
        TOp1D wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = Teuchos::ScalarTraits<RangeScalar>::zero(); // not used
        wdp.numRows = numRows_;
        wdp.numCols = Y.getNumRows();
        wdp.begs    = rbh.template addConstBuffer<size_t>(pbuf_offsets1D_);
        wdp.ends    = wdp.begs+1;
        wdp.inds    = rbh.template addConstBuffer<Ordinal>(pbuf_inds1D_);
        wdp.vals    = rbh.template addConstBuffer<Scalar>(pbuf_vals1D_);
        wdp.x       = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.y       = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<TOp1D>(0,numRHS,wdp);
      }
    }
    return;
  }


  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void DefaultDeviceSparseOps<Scalar,Ordinal,Node>::multiply(
                                Teuchos::ETransp trans, 
                                RangeScalar alpha, const MultiVector<DomainScalar,Node> &X, 
                                RangeScalar beta, MultiVector<RangeScalar,Node> &Y) const {
    // the 0 parameter means that beta is considered, and the output multivector enjoys accumulate semantics
    typedef DefaultSparseMultiplyOp1<Scalar,Ordinal,DomainScalar,RangeScalar, 0>  Op1D;
    typedef DefaultSparseTransposeMultiplyOp1<Scalar,Ordinal,DomainScalar,RangeScalar, 0> TOp1D;
    TEST_FOR_EXCEPTION(indsInit_ == false || valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::multiply(): operation not fully initialized.");
    TEST_FOR_EXCEPT(X.getNumCols() != Y.getNumCols());
    ReadyBufferHelper<Node> rbh(node_);
    if (isEmpty_ == true) {
      // Y <= alpha * 0 * X + beta * Y 
      //   <= beta * Y
      // NOTE: this neglects NaNs in X, which don't satisfy 0*NaN == 0
      //       however, X and Y may be of different size, and we need entries to determine how to mix those potential NaNs in X into Y
      //       therefore, the best we can do is scale Y to zero. Setting Y to zero would destroy NaNs in Y, which violates the semantics of the call.
      DefaultArithmetic<MultiVector<RangeScalar,Node> >::Scale(Y,beta);
    }
    else {
      if (trans == Teuchos::NO_TRANS) {
        Op1D wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = beta;
        wdp.numRows = numRows_;
        wdp.begs    = rbh.template addConstBuffer<size_t>(pbuf_offsets1D_);
        wdp.ends    = wdp.begs+1;
        wdp.inds    = rbh.template addConstBuffer<Ordinal>(pbuf_inds1D_);
        wdp.vals    = rbh.template addConstBuffer< Scalar>(pbuf_vals1D_);
        wdp.x       = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.y       = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<Op1D>(0,numRows_*numRHS,wdp);
      }
      else {
        TOp1D wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = beta;
        wdp.numRows = numRows_;
        wdp.numCols = Y.getNumRows();
        wdp.begs    = rbh.template addConstBuffer<size_t>(pbuf_offsets1D_);
        wdp.ends    = wdp.begs+1;
        wdp.inds    = rbh.template addConstBuffer<Ordinal>(pbuf_inds1D_);
        wdp.vals    = rbh.template addConstBuffer< Scalar>(pbuf_vals1D_);
        wdp.x       = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.y       = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<TOp1D>(0,numRHS,wdp);
      }
    }
    return;
  }

  template <class Scalar, class Ordinal, class Node, class S, class M>
  static void 
  PackedCRSOperatorAdaptor<Scalar,Ordinal,Node,S,N>::packFrom2D(
                      const ArrayRCP<const size_t>  &begs, 
                      const ArrayRCP<const size_t>  &ends,
                      const ArrayRCP<const Ordinal> &inds,
                      const ArrayRCP<const Scalar > &vals,
                            <size_t>  &rowOffsets,
                            ArrayRCP<Ordinal> &packedInds,
                            ArrayRCP<Scalar>  &packedVals) 
  { 
    TEST_FOR_EXECPT(true);
  }

  template <class Scalar, class Ordinal, class Node, class S, class M>
  static void 
  PackedCRSOperatorAdaptor<Scalar,Ordinal,Node,S,N>::packFrom2D(const ArrayRCP<const size_t>  &begs, 
                      const ArrayRCP<const size_t>  &ends,
                      const ArrayRCP<ArrayRCP<Ordinal> > &inds,
                      const ArrayRCP<ArrayRCP<Scalar > > &vals,
                            ArrayRCP<size_t>  &rowOffsets,
                            ArrayRCP<Ordinal> &packedInds,
                            ArrayRCP<Scalar>  &packedVals)
  {
    TEST_FOR_EXCEPT(true);
  }

}



#endif // KOKKOS_PACKEDCRSOPERATORADAPTOR_HPP
