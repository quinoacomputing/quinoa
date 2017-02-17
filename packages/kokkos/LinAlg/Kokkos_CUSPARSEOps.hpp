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

#ifndef KOKKOS_CUSPARSEOPS_HPP
#define KOKKOS_CUSPARSEOPS_HPP

#include <Teuchos_DataAccess.hpp>
#include <Teuchos_CompileTimeAssert.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <stdexcept>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_CrsMatrixBase.hpp"
#include "Kokkos_CrsGraphBase.hpp"
#include "Kokkos_CUDANodeUtils.hpp"

#include "Kokkos_MultiVector.hpp"
#include "Kokkos_NodeHelpers.hpp"
#include "Kokkos_DefaultArithmetic.hpp"

#include <cusparse_v2.h>

namespace Teuchos {
  template<>
  class TypeNameTraits<cusparseHandle_t> {
  public:
    static std::string name() { return std::string("cusparseHandle_t"); }
    static std::string concreteName( const cusparseHandle_t& ) { return std::string("cusparseHandle_t"); }
  };
  template<>
  class TypeNameTraits<cusparseMatDescr_t> {
  public:
    static std::string name() { return std::string("cusparseMatDescr_t"); }
    static std::string concreteName( const cusparseMatDescr_t& ) { return std::string("cusparseMatDescr_t"); }
  };
  template<>
  class TypeNameTraits<cusparseSolveAnalysisInfo_t> {
  public:
    static std::string name() { return std::string("cusparseSolveAnalysisInfo_t"); }
    static std::string concreteName( const cusparseSolveAnalysisInfo_t& ) { return std::string("cusparseSolveAnalysisInfo_t"); }
  };
}

namespace Kokkos {

  namespace CUSPARSEdetails {

    class Session {
      public:
        //! Must be called to initialize the CUSPARSE session.
        static void init();

        //! Accessor for CUSPARSE session.
        static RCP<const cusparseHandle_t> getHandle();

      private:
        static RCP<cusparseHandle_t> session_handle_;
    };

    RCP<cusparseMatDescr_t>           createMatDescr();
    RCP<cusparseSolveAnalysisInfo_t>  createSolveAnalysisInfo();

    class CUSPARSESessionDestroyer {
    public:
      CUSPARSESessionDestroyer();
      void free(cusparseHandle_t *ptr);
    };

    class CUSPARSEMatDescDestroyer {
    public:
      CUSPARSEMatDescDestroyer();
      void free(cusparseMatDescr_t *ptr);
    };

    class CUSPARSESolveAnalysisDestroyer {
    public:
      CUSPARSESolveAnalysisDestroyer();
      void free(cusparseSolveAnalysisInfo_t *ptr);
    };


    template <typename T>
    struct CUSPARSEUnsupportedScalar
    {
      public:
      static inline T notSupported() { return T::this_type_is_not_supported_by_CUSPARSE(); }
    };


    template <class T>
    class CUSPARSETemplateAdaptors {
      public:
      //
      static inline cusparseStatus_t
      CSRMM(cusparseHandle_t handle, cusparseOperation_t transA,
            int m, int n, int k, int nnz, const T *alpha,
            const cusparseMatDescr_t descrA, const T *csrValA,
            const int *csrRowPtrA, const int *csrColIndA,
            const T *B, int ldb,
            const T *beta, T *C, int ldc)
      { return CUSPARSEUnsupportedScalar<T>::notSupported(); }
      //
      static inline cusparseStatus_t
      CSRSM_analysis(cusparseHandle_t handle, cusparseOperation_t transA,
                     int m, int nnz, const cusparseMatDescr_t descrA,
                     const T *csrValA, const int *csrRowPtrA,
                     const int *csrColIndA, cusparseSolveAnalysisInfo_t info)
      { return CUSPARSEUnsupportedScalar<T>::notSupported(); }
      //
      static inline cusparseStatus_t
      CSRSM_solve(cusparseHandle_t handle, cusparseOperation_t transA,
                  int m, int n, const T *alpha,
                  const cusparseMatDescr_t descrA,
                  const T *csrValA, const int *csrRowPtrA,
                  const int *csrColIndA, cusparseSolveAnalysisInfo_t info,
                  const T *X, int ldx,
                  T *Y, int ldy)
      { return CUSPARSEUnsupportedScalar<T>::notSupported(); }
    };

    // float support
#ifdef HAVE_KOKKOS_CUDA_FLOAT
    template <>
    class CUSPARSETemplateAdaptors<float> {
      public:
      //
      static inline cusparseStatus_t
      CSRMM(cusparseHandle_t handle, cusparseOperation_t transA,
            int m, int n, int k, int nnz, const float *alpha,
            const cusparseMatDescr_t descrA, const float *csrValA,
            const int *csrRowPtrA, const int *csrColIndA,
            const float *B, int ldb,
            const float *beta, float *C, int ldc)
      { return cusparseScsrmm(handle,transA,m,n,k,nnz,
                              alpha,descrA,csrValA,csrRowPtrA,csrColIndA,
                              B,ldb,beta,C,ldc); }
      //
      static inline cusparseStatus_t
      CSRSM_analysis(cusparseHandle_t handle, cusparseOperation_t transA,
                     int m, int nnz, const cusparseMatDescr_t descrA,
                     const float *csrValA, const int *csrRowPtrA,
                     const int *csrColIndA, cusparseSolveAnalysisInfo_t info)
      { return cusparseScsrsm_analysis(handle,transA,m,nnz,descrA,
                                       csrValA,csrRowPtrA,csrColIndA,info); }
      //
      static inline cusparseStatus_t
      CSRSM_solve(cusparseHandle_t handle, cusparseOperation_t transA,
                  int m, int n, const float *alpha,
                  const cusparseMatDescr_t descrA,
                  const float *csrValA, const int *csrRowPtrA,
                  const int *csrColIndA, cusparseSolveAnalysisInfo_t info,
                  const float *X, int ldx,
                  float *Y, int ldy)
      { return cusparseScsrsm_solve(handle,transA,m,n,
                                    alpha,descrA,csrValA,csrRowPtrA,csrColIndA,
                                    info,X,ldx,Y,ldy); }
    };
#endif

#ifdef HAVE_KOKKOS_CUDA_DOUBLE
    // double support
    template <>
    class CUSPARSETemplateAdaptors<double> {
      public:
      //
      static inline cusparseStatus_t
      CSRMM(cusparseHandle_t handle, cusparseOperation_t transA,
            int m, int n, int k, int nnz, const double *alpha,
            const cusparseMatDescr_t descrA, const double *csrValA,
            const int *csrRowPtrA, const int *csrColIndA,
            const double *B, int ldb,
            const double *beta, double *C, int ldc)
      { return cusparseDcsrmm(handle,transA,m,n,k,nnz,
                              alpha,descrA,csrValA,csrRowPtrA,csrColIndA,
                              B,ldb,beta,C,ldc); }
      //
      static inline cusparseStatus_t
      CSRSM_analysis(cusparseHandle_t handle, cusparseOperation_t transA,
                     int m, int nnz, const cusparseMatDescr_t descrA,
                     const double *csrValA, const int *csrRowPtrA,
                     const int *csrColIndA, cusparseSolveAnalysisInfo_t info)
      { return cusparseDcsrsm_analysis(handle,transA,m,nnz,descrA,
                                       csrValA,csrRowPtrA,csrColIndA,info); }
      //
      static inline cusparseStatus_t
      CSRSM_solve(cusparseHandle_t handle, cusparseOperation_t transA,
                  int m, int n, const double *alpha,
                  const cusparseMatDescr_t descrA,
                  const double *csrValA, const int *csrRowPtrA,
                  const int *csrColIndA, cusparseSolveAnalysisInfo_t info,
                  const double *X, int ldx,
                  double *Y, int ldy)
      { return cusparseDcsrsm_solve(handle,transA,m,n,
                                    alpha,descrA,csrValA,csrRowPtrA,csrColIndA,
                                    info,X,ldx,Y,ldy); }
    };
#endif

#ifdef HAVE_KOKKOS_CUDA_COMPLEX_FLOAT
    // complex<float> support
    template <>
    class CUSPARSETemplateAdaptors<std::complex<float> > {
      public:
      //
      static inline cusparseStatus_t
      CSRMM(cusparseHandle_t handle, cusparseOperation_t transA,
            int m, int n, int k, int nnz, const cuComplex *alpha,
            const cusparseMatDescr_t descrA, const cuComplex *csrValA,
            const int *csrRowPtrA, const int *csrColIndA,
            const cuComplex *B, int ldb,
            const cuComplex *beta, cuComplex *C, int ldc)
      { return cusparseCcsrmm(handle,transA,m,n,k,nnz,
                              alpha,descrA,csrValA,csrRowPtrA,csrColIndA,
                              B,ldb,beta,C,ldc); }
      //
      static inline cusparseStatus_t
      CSRSM_analysis(cusparseHandle_t handle, cusparseOperation_t transA,
                     int m, int nnz, const cusparseMatDescr_t descrA,
                     const cuComplex *csrValA, const int *csrRowPtrA,
                     const int *csrColIndA, cusparseSolveAnalysisInfo_t info)
      { return cusparseCcsrsm_analysis(handle,transA,m,nnz,descrA,
                                       csrValA,csrRowPtrA,csrColIndA,info); }
      //
      static inline cusparseStatus_t
      CSRSM_solve(cusparseHandle_t handle, cusparseOperation_t transA,
                  int m, int n, const cuComplex *alpha,
                  const cusparseMatDescr_t descrA,
                  const cuComplex *csrValA, const int *csrRowPtrA,
                  const int *csrColIndA, cusparseSolveAnalysisInfo_t info,
                  const cuComplex *X, int ldx,
                  cuComplex *Y, int ldy)
      { return cusparseCcsrsm_solve(handle,transA,m,n,
                                    alpha,descrA,csrValA,csrRowPtrA,csrColIndA,
                                    info,X,ldx,Y,ldy); }
    };
#endif

#ifdef HAVE_KOKKOS_CUDA_COMPLEX_DOUBLE
    // complex<double> support
    template <>
    class CUSPARSETemplateAdaptors<std::complex<double> > {
      public:
      //
      static inline cusparseStatus_t
      CSRMM(cusparseHandle_t handle, cusparseOperation_t transA,
            int m, int n, int k, int nnz, const cuDoubleComplex *alpha,
            const cusparseMatDescr_t descrA, const cuDoubleComplex *csrValA,
            const int *csrRowPtrA, const int *csrColIndA,
            const cuDoubleComplex *B, int ldb,
            const cuDoubleComplex *beta, cuDoubleComplex *C, int ldc)
      { return cusparseZcsrmm(handle,transA,m,n,k,nnz,
                              alpha,descrA,csrValA,csrRowPtrA,csrColIndA,
                              B,ldb,beta,C,ldc); }
      //
      static inline cusparseStatus_t
      CSRSM_analysis(cusparseHandle_t handle, cusparseOperation_t transA,
                     int m, int nnz, const cusparseMatDescr_t descrA,
                     const cuDoubleComplex *csrValA, const int *csrRowPtrA,
                     const int *csrColIndA, cusparseSolveAnalysisInfo_t info)
      { return cusparseZcsrsm_analysis(handle,transA,m,nnz,descrA,
                                       csrValA,csrRowPtrA,csrColIndA,info); }
      //
      static inline cusparseStatus_t
      CSRSM_solve(cusparseHandle_t handle, cusparseOperation_t transA,
                  int m, int n, const cuDoubleComplex *alpha,
                  const cusparseMatDescr_t descrA,
                  const cuDoubleComplex *csrValA, const int *csrRowPtrA,
                  const int *csrColIndA, cusparseSolveAnalysisInfo_t info,
                  const cuDoubleComplex *X, int ldx,
                  cuDoubleComplex *Y, int ldy)
      { return cusparseZcsrsm_solve(handle,transA,m,n,
                                    alpha,descrA,csrValA,csrRowPtrA,csrColIndA,
                                    info,X,ldx,Y,ldy); }
    };
#endif

  } // end of namespace CUSPARSEdetails

  //! \class CUSPARSECrsGraph
  /** \brief CRS sparse graph class supporting the CUSPARSE library.
  */
  template <class Node>
  class CUSPARSECrsGraph : public CrsGraphBase<int,Node>
  {
    public:
      CUSPARSECrsGraph(size_t numRows, const RCP<Node> &node, const RCP<ParameterList> &params);
      bool isEmpty() const;
      void setStructure(const ArrayRCP<const size_t>  &ptrs,
                        const ArrayRCP<const int> &inds);
      void setDeviceData(const ArrayRCP<const int> &devptrs, const ArrayRCP<const int> &devinds);
      inline ArrayRCP<const size_t> getPointers() const;
      inline ArrayRCP<const int>    getIndices() const;
      inline ArrayRCP<const int>    getDevPointers() const;
      inline ArrayRCP<const int>    getDevIndices() const;
      inline bool isInitialized() const;
      void setMatDesc(Teuchos::EUplo uplo, Teuchos::EDiag diag);
      void getMatDesc(Teuchos::EUplo &uplo, Teuchos::EDiag &diag) const;
      RCP<cusparseMatDescr_t> getMatDesc() const;
    private:
      bool isInitialized_;
      bool isEmpty_;
      // cusparse matrix description handle
      RCP<cusparseMatDescr_t> matdescr_;
      Teuchos::EUplo uplo_;
      Teuchos::EDiag diag_;
      // graph data
      ArrayRCP<const size_t>  host_rowptrs_;
      ArrayRCP<const int>     dev_rowptrs_;
      ArrayRCP<const int>     host_colinds_, dev_colinds_;
      // TODO: add CSC data, for efficient tranpose multiply
  };

  //! \class CUSPARSECrsMatrix
  /** \brief CRS sparse matrix class supporting the CUSPARSE library.
  */
  template <class Scalar,
            class Node>
  class CUSPARSECrsMatrix : public CrsMatrixBase<Scalar,int,Node>
  {
    public:
      CUSPARSECrsMatrix(const RCP<const CUSPARSECrsGraph<Node> > &graph, const RCP<ParameterList> &params);
      void setValues(const ArrayRCP<const Scalar> &vals);
      void setDeviceData(const ArrayRCP<const Scalar> &devvals);
      inline ArrayRCP<const Scalar> getValues() const;
      inline ArrayRCP<const Scalar> getDevValues() const;
      inline bool isInitialized() const;
      void setAnalyses(const RCP<cusparseSolveAnalysisInfo_t> &analysisNoTrans,
                       const RCP<cusparseSolveAnalysisInfo_t> &analysisTrans,
                       const RCP<cusparseSolveAnalysisInfo_t> &analysisConjTrans);
      void getAnalyses(RCP<cusparseSolveAnalysisInfo_t> &analysisNoTrans,
                       RCP<cusparseSolveAnalysisInfo_t> &analysisTrans,
                       RCP<cusparseSolveAnalysisInfo_t> &analysisConjTrans) const;
    private:
      bool isInitialized_;
      // cusparse analysis handles
      RCP<cusparseSolveAnalysisInfo_t> analysisNoTrans_, analysisConjTrans_, analysisTrans_;
      // matrix data
      ArrayRCP<const Scalar> vals_, dev_vals_;
      // TODO: add CSC data, for efficient transpose multiply
  };

  template <class Node>
  CUSPARSECrsGraph<Node>::CUSPARSECrsGraph(size_t numRows, const RCP<Node> &node, const RCP<ParameterList> &params)
  : CrsGraphBase<int,Node>(numRows,node,params)
  , isInitialized_(false)
  , isEmpty_(false)
  {
    CUSPARSEdetails::Session::init();
    // Make sure that users only specialize for Kokkos Node types that are CUDA Nodes
    Teuchos::CompileTimeAssert<Node::isCUDANode == false> cta; (void)cta;
  }

  template <class Node>
  bool CUSPARSECrsGraph<Node>::isEmpty() const
  { return isEmpty_; }

  template <class Node>
  void CUSPARSECrsGraph<Node>::setStructure(
                      const ArrayRCP<const size_t>  &ptrs,
                      const ArrayRCP<const int> &inds)
  {
    std::string tfecfFuncName("setStructure(ptrs,inds)");
    const size_t numrows = this->getNumRows();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        (size_t)ptrs.size() != numrows+1
        || ptrs[0] != 0
        || (size_t)inds.size() != ptrs[numrows],
        std::runtime_error, " graph data not coherent."
    )
    const size_t numEntries = ptrs[numrows];
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == true,
        std::runtime_error, " matrix has already been initialized"
    )
    if (numrows == 0 || numEntries == 0) isEmpty_ = true;
    host_rowptrs_ = ptrs;
    host_colinds_ = inds;
    isInitialized_ = true;
  }

  template <class Node>
  ArrayRCP<const size_t> CUSPARSECrsGraph<Node>::getPointers() const
  { return host_rowptrs_; }

  template <class Node>
  ArrayRCP<const int> CUSPARSECrsGraph<Node>::getIndices() const
  { return host_colinds_; }

  template <class Node>
  ArrayRCP<const int> CUSPARSECrsGraph<Node>::getDevPointers() const
  { return dev_rowptrs_; }

  template <class Node>
  ArrayRCP<const int> CUSPARSECrsGraph<Node>::getDevIndices() const
  { return dev_colinds_; }

  template <class Node>
  void CUSPARSECrsGraph<Node>::setDeviceData(const ArrayRCP<const int> &devptrs,
                                             const ArrayRCP<const int>    &devinds)
  { dev_rowptrs_ = devptrs; dev_colinds_ = devinds; }

  template <class Node>
  bool CUSPARSECrsGraph<Node>::isInitialized() const
  { return isInitialized_; }

  template <class Node>
  void CUSPARSECrsGraph<Node>::setMatDesc(Teuchos::EUplo uplo, Teuchos::EDiag diag)
  {
    std::string tfecfFuncName("setMatDesc()");
    //
    uplo_ = uplo;
    diag_ = diag;
    //
    matdescr_ = CUSPARSEdetails::createMatDescr();
    cusparseStatus_t status = cusparseSetMatIndexBase(*matdescr_, CUSPARSE_INDEX_BASE_ZERO);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(status != CUSPARSE_STATUS_SUCCESS, std::runtime_error, "error setting matrix descriptor (index base).")
    // upper or lower
    if (uplo == Teuchos::UPPER_TRI) {
      status = cusparseSetMatFillMode(*matdescr_, CUSPARSE_FILL_MODE_UPPER);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(status != CUSPARSE_STATUS_SUCCESS, std::runtime_error, "error setting matrix descriptor (upper).")
      status = cusparseSetMatType(*matdescr_, CUSPARSE_MATRIX_TYPE_TRIANGULAR);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(status != CUSPARSE_STATUS_SUCCESS, std::runtime_error, "error setting matrix descriptor (triangular).")
    }
    else if (uplo == Teuchos::LOWER_TRI) {
      status = cusparseSetMatFillMode(*matdescr_, CUSPARSE_FILL_MODE_LOWER);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(status != CUSPARSE_STATUS_SUCCESS, std::runtime_error, "error setting matrix descriptor (lower).")
      status = cusparseSetMatType(*matdescr_, CUSPARSE_MATRIX_TYPE_TRIANGULAR);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(status != CUSPARSE_STATUS_SUCCESS, std::runtime_error, "error setting matrix descriptor (triangular).")
    }
    else {
      status = cusparseSetMatType(*matdescr_, CUSPARSE_MATRIX_TYPE_GENERAL);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(status != CUSPARSE_STATUS_SUCCESS, std::runtime_error, "error setting matrix descriptor (general).")
    }
    if (diag == Teuchos::UNIT_DIAG) {
      status = cusparseSetMatDiagType(*matdescr_, CUSPARSE_DIAG_TYPE_UNIT);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(status != CUSPARSE_STATUS_SUCCESS, std::runtime_error, "error setting matrix descriptor (unit).")
    }
    else {
      status = cusparseSetMatDiagType(*matdescr_, CUSPARSE_DIAG_TYPE_NON_UNIT);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(status != CUSPARSE_STATUS_SUCCESS, std::runtime_error, "error setting matrix descriptor (non-unit).")
    }
  }

  template <class Node>
  void CUSPARSECrsGraph<Node>::getMatDesc(Teuchos::EUplo &uplo, Teuchos::EDiag &diag) const
  {
    uplo = uplo_;
    diag = diag_;
  }

  template <class Node>
  RCP<cusparseMatDescr_t> CUSPARSECrsGraph<Node>::getMatDesc() const
  {
    return matdescr_;
  }

  template <class Scalar, class Node>
  CUSPARSECrsMatrix<Scalar,Node>::CUSPARSECrsMatrix(const RCP<const CUSPARSECrsGraph<Node> > &graph, const RCP<ParameterList> &params)
  : CrsMatrixBase<Scalar,int,Node>(graph,params)
  , isInitialized_(false)
  {
    // Make sure that users only specialize for Kokkos Node types that are CUDA Nodes
    Teuchos::CompileTimeAssert<Node::isCUDANode == false> cta; (void)cta;
  }

  template <class Scalar, class Node>
  void CUSPARSECrsMatrix<Scalar,Node>::setValues(const ArrayRCP<const Scalar> &vals)
  {
    std::string tfecfFuncName("setValues(vals)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == true,
        std::runtime_error, " matrix is already initialized."
    )
    vals_ = vals;
    isInitialized_ = true;
  }

  template <class Scalar, class Node>
  ArrayRCP<const Scalar> CUSPARSECrsMatrix<Scalar,Node>::getValues() const
  { return vals_; }

  template <class Scalar, class Node>
  bool CUSPARSECrsMatrix<Scalar,Node>::isInitialized() const
  { return isInitialized_; }

  template <class Scalar, class Node>
  ArrayRCP<const Scalar> CUSPARSECrsMatrix<Scalar,Node>::getDevValues() const
  { return dev_vals_; }

  template <class Scalar, class Node>
  void CUSPARSECrsMatrix<Scalar,Node>::setDeviceData(const ArrayRCP<const Scalar> &devvals)
  { dev_vals_ = devvals; }

  template <class Scalar, class Node>
  void CUSPARSECrsMatrix<Scalar,Node>::setAnalyses(const RCP<cusparseSolveAnalysisInfo_t> &analysisNoTrans,
                                                   const RCP<cusparseSolveAnalysisInfo_t> &analysisTrans,
                                                   const RCP<cusparseSolveAnalysisInfo_t> &analysisConjTrans)
  {
    analysisNoTrans_   = analysisNoTrans;
    analysisTrans_     = analysisTrans;
    analysisConjTrans_ = analysisConjTrans;
  }

  template <class Scalar, class Node>
  void CUSPARSECrsMatrix<Scalar,Node>::getAnalyses(RCP<cusparseSolveAnalysisInfo_t> &analysisNoTrans,
                                                   RCP<cusparseSolveAnalysisInfo_t> &analysisTrans,
                                                   RCP<cusparseSolveAnalysisInfo_t> &analysisConjTrans) const
  {
    analysisNoTrans   = analysisNoTrans_;
    analysisTrans     = analysisTrans_;
    analysisConjTrans = analysisConjTrans_;
  }

  /// \class CUSPARSEOps
  /// \brief Default implementation of sparse matrix-vector multiply
  ///   and solve routines, for host-based Kokkos Node types.
  /// \ingroup kokkos_crs_ops
  ///
  /// \tparam Scalar The type of entries of the sparse matrix.
  /// \tparam Node The Kokkos Node type.
  template <class Scalar, class Node>
  class CUSPARSEOps {
  public:
    //@{
    //! @name Typedefs and structs

    //! The type of the individual entries of the sparse matrix.
    typedef Scalar  scalar_type;
    //! The type of the (local) indices describing the structure of the sparse matrix.
    typedef int     ordinal_type;
    //! The Kokkos Node type.
    typedef Node    node_type;
    //! The type of this object, the sparse operator object
    typedef CUSPARSEOps<Scalar,Node> sparse_ops_type;

    /** \brief Typedef for local graph class */
    template <class O, class N>
    struct graph {
      typedef CUSPARSECrsGraph<N> graph_type;
    };

    /** \brief Typedef for local matrix class */
    template <class S, class O, class N>
    struct matrix {
      typedef CUSPARSECrsMatrix<S,N> matrix_type;
    };

    /// \brief Sparse operations type for a different scalar type.
    ///
    /// The bind_scalar struct defines the type responsible for sparse
    /// operations for a scalar type S2, which may be different from
    /// \c Scalar.
    ///
    /// Use by Tpetra CrsMatrix to bind a potentially "void" scalar type to the appropriate scalar.
    ///
    /// This always specifies a specialization of \c
    /// CUSPARSEOps, regardless of the scalar type S2.
    ///
    /// \tparam S2 A scalar type possibly different from \c Scalar.
    template <class S2>
    struct bind_scalar {
      typedef CUSPARSEOps<S2,Node> other_type;
    };

    //@}
    //! @name Constructors/Destructor
    //@{

    //! Constructor accepting and retaining a node object.
    CUSPARSEOps(const RCP<Node> &node);

    //! Destructor
    ~CUSPARSEOps();

    //@}
    //! @name Accessor routines.
    //@{

    //! The Kokkos Node with which this object was instantiated.
    RCP<Node> getNode() const;

    //@}

    //! @name Initialization of graph and matrix
    //@{

    //! \brief Allocate and initialize the storage for the matrix values.
    static ArrayRCP<size_t> allocRowPtrs(const RCP<Node> &node, const ArrayView<const size_t> &rowPtrs);

    //! \brief Allocate and initialize the storage for a sparse graph.
    template <class T>
    static ArrayRCP<T> allocStorage(const RCP<Node> &node, const ArrayView<const size_t> &ptrs);

    //! Finalize a graph is null for CUSPARSE.
    static void finalizeGraph(Teuchos::EUplo uplo, Teuchos::EDiag diag, CUSPARSECrsGraph<Node> &graph, const RCP<ParameterList> &params);

    //! Finalize the matrix of an already-finalized graph.
    static void finalizeMatrix(const CUSPARSECrsGraph<Node> &graph, CUSPARSECrsMatrix<Scalar,Node> &matrix, const RCP<ParameterList> &params);

    //! Finalize a graph and a matrix.
    static void finalizeGraphAndMatrix(Teuchos::EUplo uplo, Teuchos::EDiag diag, CUSPARSECrsGraph<Node> &graph, CUSPARSECrsMatrix<Scalar,Node> &matrix, const RCP<ParameterList> &params);

    //! Initialize sparse operations with a graph and matrix
    void setGraphAndMatrix(const RCP<const CUSPARSECrsGraph<Node> > &graph,
                           const RCP<const CUSPARSECrsMatrix<Scalar,Node> > &matrix);

    //@}
    //! @name Computational methods
    //@{

    /// \brief Y := alpha * Op(A) * X.
    ///
    /// Apply the local sparse matrix A (or its transpose or conjugate
    /// transpose) to a multivector X, overwriting Y with the result.
    /// Op(A) means A, the transpose of A, or the conjugate transpose
    /// of A, depending on the \c trans argument.
    ///
    /// \tparam DomainScalar The type of entries in the input
    ///   multivector X.  This may differ from the type of entries in
    ///   A or in Y.
    ///
    /// \tparam RangeScalar The type of entries in the output
    ///   multivector Y.  This may differ from the type of entries in
    ///   A or in X.
    ///
    /// \param trans [in] Whether to apply the matrix, its transpose,
    ///   or its conjugate transpose (if applicable).
    ///
    /// \param alpha [in] Scalar constant \f$\alpha\f$ by which to
    ///   multiply the result: \f$Y := \alpha A X\f$.
    ///
    /// \param X [in] Input multivector.
    ///
    /// \param Y [out] Result multivector.
    template <class DomainScalar, class RangeScalar>
    void
    multiply (Teuchos::ETransp trans,
              RangeScalar alpha,
              const MultiVector<DomainScalar,Node> &X,
              MultiVector<RangeScalar,Node> &Y) const;

    /// \brief Y := Y + alpha * Op(A) * X.
    ///
    /// Apply the local sparse matrix A (or its transpose or conjugate
    /// transpose) to a multivector X, accumulating the result into Y.
    /// Op(A) means A, the transpose of A, or the conjugate transpose
    /// of A, depending on the \c trans argument.
    ///
    /// \tparam DomainScalar The type of entries in the input
    ///   multivector X.  This may differ from the type of entries in
    ///   A or in Y.
    ///
    /// \tparam RangeScalar The type of entries in the output
    ///   multivector Y.  This may differ from the type of entries in
    ///   A or in X.
    ///
    /// \param trans [in] Whether to apply the matrix, its transpose,
    ///   or its conjugate transpose (if applicable).
    ///
    /// \param alpha [in] Scalar constant \f$\alpha\f$ by which to
    ///   multiply the result: \f$Y := Y + \alpha A X\f$.
    ///
    /// \param X [in] Input multivector.
    ///
    /// \param Y [in/out] Result multivector.
    template <class DomainScalar, class RangeScalar>
    void
    multiply (Teuchos::ETransp trans,
              RangeScalar alpha,
              const MultiVector<DomainScalar,Node> &X,
              RangeScalar beta,
              MultiVector<RangeScalar,Node> &Y) const;

    /// \brief Solve Y = Op(A) X for X, where we assume A is triangular.
    ///
    /// Solve the (upper or lower) triangular system Y = Op(A) X.
    /// Op(A) means A, the transpose of A, or the conjugate transpose
    /// of A, depending on the \c trans argument.
    ///
    /// \tparam DomainScalar The type of entries in the input
    ///   multivector X.  This may differ from the type of entries in
    ///   A or in Y.
    ///
    /// \tparam RangeScalar The type of entries in the output
    ///   multivector Y.  This may differ from the type of entries in
    ///   A or in X.
    ///
    /// \param trans [in] Whether to solve with the matrix, its
    ///   transpose, or its conjugate transpose (if applicable).
    ///
    /// \param uplo [in] UPPER_TRI if the matrix is upper triangular,
    ///   else LOWER_TRI if the matrix is lower triangular.
    ///
    /// \param diag [in] UNIT_DIAG if the matrix has an implicit unit diagonal,
    ///   else NON_UNIT_DIAG (diagonal entries are explicitly stored in the matrix).
    ///
    /// \param Y [in] Input multivector.
    ///
    /// \param X [out] Result multivector.
    template <class DomainScalar, class RangeScalar>
    void
    solve (Teuchos::ETransp trans,
           const MultiVector<DomainScalar,Node> &Y,
           MultiVector<RangeScalar,Node> &X) const;

    //@}

  private:
    //! Copy constructor (protected and unimplemented)
    CUSPARSEOps(const CUSPARSEOps& source);

    //! The Kokkos Node instance given to this object's constructor.
    RCP<Node> node_;

    int numRows_, numNZ_;
    bool isInitialized_;

    ArrayRCP<const int> rowPtrs_, colInds_;
    ArrayRCP<const Scalar> rowVals_;
    RCP<cusparseMatDescr_t>          matdescr_;
    RCP<cusparseSolveAnalysisInfo_t> aiNoTrans_, aiTrans_, aiConjTrans_;
  };


  // ======= matrix finalization ===========
  template <class Scalar, class Node>
  void CUSPARSEOps<Scalar,Node>::finalizeGraph(Teuchos::EUplo uplo, Teuchos::EDiag diag,
                                               CUSPARSECrsGraph<Node> &graph, const RCP<ParameterList> &params)
  {
    const std::string prefix("finalizeGraph()");
    RCP<Node> node = graph.getNode();
    TEUCHOS_TEST_FOR_EXCEPTION(
        graph.isInitialized() == false,
        std::runtime_error, prefix << ": graph has not yet been initialized."
    )
    // diag: have to allocate and indicate
    ArrayRCP<int>    devinds, devptrs;
    ArrayRCP<const int>    hostinds = graph.getIndices();
    ArrayRCP<const size_t> hostptrs = graph.getPointers();
    const size_t numRows = graph.getNumRows();
    // set description
    graph.setMatDesc(uplo,diag);
    // allocate and initialize data
    if (diag == Teuchos::UNIT_DIAG) {
      // cuSPARSE, unfortunately, always assumes that the diagonal
      // entries are present in the storage.  Therefore, this flag
      // only specifies whether they are considered or not; they are
      // assumed to be present, and neglecting them will result in
      // incorrect behavior (causing the adjacent entry to be
      // neglected instead).  Therefore, because our API doesn't give
      // us explicit diagonal entries if diag == Teuchos::UNIT_DIAG,
      // we must allocate space for them.
      const size_t numnz = hostinds.size() + numRows;
      devptrs = node->template allocBuffer<int>( numRows+1 );
      if (numnz) devinds = node->template allocBuffer<int>( numnz );
      ArrayRCP<int> h_devptrs = node->viewBufferNonConst(WriteOnly, numRows+1, devptrs);
      ArrayRCP<int> h_devinds = node->viewBufferNonConst(WriteOnly, numnz,     devinds);
      for (size_t r=0; r < numRows; ++r) {
        h_devptrs[r] = (int)(hostptrs[r]+r);
        if (uplo == Teuchos::LOWER_TRI) {
          // leave one space at the end
          std::copy (hostinds.begin()+hostptrs[r],
                     hostinds.begin()+hostptrs[r+1],
                     h_devinds.begin()+h_devptrs[r]);
        }
        else {
          // leave one space at the beginning
          std::copy (hostinds.begin()+hostptrs[r],
                     hostinds.begin()+hostptrs[r+1],
                     h_devinds.begin()+h_devptrs[r]+1);
        }
      }
      h_devptrs[numRows] = (int)(hostptrs[numRows]+numRows);
      // copy back
      h_devptrs = null;
      h_devinds = null;
    }
    else {
      // our format == their format; just allocate and copy
      const size_t numnz = hostinds.size();
      devptrs = node->template allocBuffer<int>( numRows+1 );
      ArrayRCP<int> h_devptrs = node->viewBufferNonConst(WriteOnly, numRows+1, devptrs);
      // have to copy on host because of the cast from size_t to int
      std::copy(hostptrs.begin(), hostptrs.end(), h_devptrs.begin());
      h_devptrs = null;
      if (numnz) {
        devinds = node->template allocBuffer<int>( numnz );
        node->copyToBuffer( numnz,     hostinds(), devinds );
      }
    }
    // set the data
    graph.setDeviceData(devptrs,devinds);
  }

  // ======= matrix finalization ===========
  template <class Scalar, class Node>
  void CUSPARSEOps<Scalar,Node>::
  finalizeMatrix (const CUSPARSECrsGraph<Node> &graph,
                  CUSPARSECrsMatrix<Scalar,Node> &matrix,
                  const RCP<ParameterList> &params)
  {
    std::string FuncName("Kokkos::CUSPARSEOps::finalizeMatrix()");
    RCP<Node> node = graph.getNode();
    TEUCHOS_TEST_FOR_EXCEPTION(
        matrix.isInitialized() == false,
        std::runtime_error, FuncName << ": matrix has not yet been initialized."
    )
    // diag: have to allocate and indicate
    Teuchos::EUplo uplo;
    Teuchos::EDiag diag;
    graph.getMatDesc(uplo,diag);
    ArrayRCP<Scalar>       devvals;
    ArrayRCP<const size_t> hostptrs = graph.getPointers();
    ArrayRCP<const Scalar> hostvals = matrix.getValues();
    const size_t numRows = graph.getNumRows();
    if (diag == Teuchos::UNIT_DIAG) {
      // cuSPARSE, unfortunately, always assumes that the diagonal
      // entries are present in the storage.  Therefore, this flag
      // only specifies whether they are considered or not; they are
      // assumed to be present, and neglecting them will result in
      // incorrect behavior (causing the adjacent entry to be
      // neglected instead).  Therefore, because our API doesn't give
      // us diagonal entries if diag == Teuchos::UNIT_DIAG, we must
      // allocate space for them.
      const size_t numnz = hostptrs[numRows] + numRows;
      devvals = node->template allocBuffer<Scalar>( numnz );
      ArrayRCP<Scalar> h_devvals = node->viewBufferNonConst(WriteOnly, numnz, devvals);
      for (size_t r=0; r < numRows; ++r) {
        if (uplo == Teuchos::LOWER_TRI) {
          // leave one space at the end
          std::copy( hostvals.begin()+hostptrs[r], hostvals.begin()+hostptrs[r+1], h_devvals.begin()+hostptrs[r]+r );
        }
        else {
          // leave one space at the beginning
          std::copy( hostvals.begin()+hostptrs[r], hostvals.begin()+hostptrs[r+1], h_devvals.begin()+hostptrs[r]+r+1 );
        }
      }
      // copy back
      h_devvals = null;
    }
    else {
      const size_t numnz = hostptrs[numRows];
      if (numnz) {
        devvals = node->template allocBuffer<Scalar>( numnz );
        node->copyToBuffer( numnz,     hostvals(), devvals );
      }
    }
    matrix.setDeviceData(devvals);

    const int numnz = (int)devvals.size();
    RCP<cusparseMatDescr_t> descr = graph.getMatDesc();
    ArrayRCP<const int> devptrs = graph.getDevPointers(),
                        devinds = graph.getDevIndices();

    RCP<const cusparseHandle_t> hndl = CUSPARSEdetails::Session::getHandle();
    // look at the parameter list and do any analyses requested for solves
    RCP<cusparseSolveAnalysisInfo_t> ai_non, ai_trans, ai_conj;
    if (params != null && params->get("Prepare Solve",false)) {
      ai_non = CUSPARSEdetails::createSolveAnalysisInfo();
      cusparseStatus_t stat = CUSPARSEdetails::CUSPARSETemplateAdaptors<Scalar>::CSRSM_analysis(
          *hndl, CUSPARSE_OPERATION_NON_TRANSPOSE,numRows,
          numnz, *descr, devvals.getRawPtr(), devptrs.getRawPtr(), devinds.getRawPtr(),
          *ai_non);
      TEUCHOS_TEST_FOR_EXCEPTION(stat != CUSPARSE_STATUS_SUCCESS, std::runtime_error,
          FuncName << ": CSRSM_analysis(non-trans) returned error " << stat);
    }
    if (params != null && params->get("Prepare Tranpose Solve",false)) {
      ai_trans = CUSPARSEdetails::createSolveAnalysisInfo();
      cusparseStatus_t stat = CUSPARSEdetails::CUSPARSETemplateAdaptors<Scalar>::CSRSM_analysis(
          *hndl, CUSPARSE_OPERATION_TRANSPOSE,numRows,
          numnz, *descr, devvals.getRawPtr(), devptrs.getRawPtr(), devinds.getRawPtr(),
          *ai_trans);
      TEUCHOS_TEST_FOR_EXCEPTION(stat != CUSPARSE_STATUS_SUCCESS, std::runtime_error,
          FuncName << ": CSRSM_analysis(trans) returned error " << stat);
    }
    if (params != null && params->get("Prepare Conjugate Tranpose Solve",false)) {
      ai_conj = CUSPARSEdetails::createSolveAnalysisInfo();
      cusparseStatus_t stat = CUSPARSEdetails::CUSPARSETemplateAdaptors<Scalar>::CSRSM_analysis(
          *hndl, CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE,numRows,
          numnz, *descr, devvals.getRawPtr(), devptrs.getRawPtr(), devinds.getRawPtr(),
          *ai_conj);
      TEUCHOS_TEST_FOR_EXCEPTION(stat != CUSPARSE_STATUS_SUCCESS, std::runtime_error,
          FuncName << ": CSRSM_analysis(conj-trans) returned error " << stat);
    }
    matrix.setAnalyses( ai_non, ai_trans, ai_conj );
    //
    if (params != null && params->get("Prepare Transpose Multiply",false)) {
      // finish: compute CSC for transpose
      TEUCHOS_TEST_FOR_EXCEPT(true);
    }
  }

  // ======= graph and matrix finalization ===========
  template <class Scalar, class Node>
  void CUSPARSEOps<Scalar,Node>::finalizeGraphAndMatrix(Teuchos::EUplo uplo, Teuchos::EDiag diag,
                                                        CUSPARSECrsGraph<Node> &graph, CUSPARSECrsMatrix<Scalar,Node> &matrix,
                                                        const RCP<ParameterList> &params)
  {
    std::string FuncName("Kokkos::CUSPARSEOps::finalizeGraphAndMatrix(graph,matrix,params)");
    TEUCHOS_TEST_FOR_EXCEPTION(
        graph.isInitialized() == false,
        std::runtime_error, FuncName << ": graph has not yet been initialized."
    )
    TEUCHOS_TEST_FOR_EXCEPTION(
        matrix.isInitialized() == false,
        std::runtime_error, FuncName << ": matrix has not yet been initialized."
    )
    // no benefit to doing them together; do them separately
    finalizeGraph(uplo,diag,graph,params);
    finalizeMatrix(graph,matrix,params);
  }


  // ======= pointer allocation ===========
  template <class Scalar, class Node>
  ArrayRCP<size_t>
  CUSPARSEOps<Scalar,Node>::allocRowPtrs(const RCP<Node> &/*node*/, const ArrayView<const size_t> &numEntriesPerRow)
  {
    // alloc page-locked ("pinned") memory on the host, specially allocated and specially deallocated
    CUDANodeHostPinnedDeallocator<size_t> dealloc;
    ArrayRCP<size_t> ptrs = dealloc.alloc(numEntriesPerRow.size() + 1);
    ptrs[0] = 0;
    std::partial_sum( numEntriesPerRow.getRawPtr(), numEntriesPerRow.getRawPtr()+numEntriesPerRow.size(), ptrs.begin()+1 );
    return ptrs;
  }

  // ======= other allocation ===========
  template <class Scalar, class Node>
  template <class T>
  ArrayRCP<T>
  CUSPARSEOps<Scalar,Node>::allocStorage(const RCP<Node> &/*node*/, const ArrayView<const size_t> &rowPtrs)
  {
    // alloc page-locked ("pinned") memory on the host, specially allocated and specially deallocated
    const size_t totalNumEntries = *(rowPtrs.end()-1);
    CUDANodeHostPinnedDeallocator<T> dealloc;
    ArrayRCP<T> buf = dealloc.alloc(totalNumEntries);
    std::fill(buf.begin(), buf.end(), Teuchos::ScalarTraits<T>::zero() );
    return buf;
  }
  template<class Scalar, class Node>
  CUSPARSEOps<Scalar,Node>::CUSPARSEOps(const RCP<Node> &node)
  : node_(node)
  , numRows_(0)
  , numNZ_(0)
  , isInitialized_(false)
  {
    CUSPARSEdetails::Session::init();
    // Make sure that users only specialize for Kokkos Node types that are CUDA Nodes
    Teuchos::CompileTimeAssert<Node::isCUDANode == false> cta; (void)cta;
  }

  template<class Scalar, class Node>
  CUSPARSEOps<Scalar,Node>::~CUSPARSEOps()
  { }

  template <class Scalar, class Node>
  RCP<Node> CUSPARSEOps<Scalar,Node>::getNode() const {
    return node_;
  }

  template <class Scalar, class Node>
  void CUSPARSEOps<Scalar,Node>::setGraphAndMatrix(const RCP<const CUSPARSECrsGraph<Node> > &graph, const RCP<const CUSPARSECrsMatrix<Scalar,Node> > &matrix)
  {
    std::string tfecfFuncName("setGraphAndMatrix(graph,matrix)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == true,
        std::runtime_error, " operators already initialized.");
    // get cusparse data from the matrix
    numRows_ = graph->getNumRows();
    matdescr_ = graph->getMatDesc();
    rowPtrs_ = graph->getDevPointers();
    colInds_ = graph->getDevIndices();
    rowVals_ = matrix->getDevValues();
    numNZ_ = colInds_.size();
    matrix->getAnalyses( aiNoTrans_, aiTrans_, aiConjTrans_ );
    isInitialized_ = true;
  }

  template <class Scalar, class Node>
  template <class DomainScalar, class RangeScalar>
  void CUSPARSEOps<Scalar,Node>::solve(Teuchos::ETransp trans,
                                       const MultiVector<DomainScalar,Node> &Y,
                                             MultiVector< RangeScalar,Node> &X) const
  {
    // CUSPARSE doesn't support mixed precision; partial specialize, then nix the generic versions
    Teuchos::CompileTimeAssert<Teuchos::TypeTraits::is_same<DomainScalar,Scalar>::value == false ||
                               Teuchos::TypeTraits::is_same< RangeScalar,Scalar>::value == false > cta; (void)cta;
    //
    std::string tfecfFuncName("solve()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == false,
        std::runtime_error, "sparse operators have not been initialized with graph and matrix data; call setGraphAndMatrix() first.")
    // get pointers,stride from X and Y
    int stride_x = (int)X.getStride(),
        stride_y = (int)Y.getStride();
    const Scalar * data_y = Y.getValues().getRawPtr();
    Scalar * data_x = X.getValuesNonConst().getRawPtr();
    const int numMatRows = X.getNumRows();
    const int numRHS     = X.getNumCols();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        X.getNumRows() != Y.getNumRows(),
        std::runtime_error, "X and Y do not have the same number of row vectors.")
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        X.getNumCols() != Y.getNumCols(),
        std::runtime_error, "X and Y do not have the same number of column vectors.")
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        Y.getNumRows() != numRows_,
        std::runtime_error, "Y does not have the same number of rows as does the matrix.")
    //
    RCP<cusparseSolveAnalysisInfo_t> solveInfo;
    cusparseOperation_t                     op;
    if (trans == Teuchos::NO_TRANS) {
      solveInfo = aiNoTrans_;
      op = CUSPARSE_OPERATION_NON_TRANSPOSE;
    }
    else if (trans == Teuchos::TRANS) {
      solveInfo = aiTrans_;
      op = CUSPARSE_OPERATION_TRANSPOSE;
    }
    else if (trans == Teuchos::CONJ_TRANS) {
      solveInfo = aiConjTrans_;
      op = CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE;
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          true, std::runtime_error, "invalid transformation (none ofe of NO_TRANS, TRANS or CONJ_TRANS): " << trans);
    }
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      solveInfo == null, std::runtime_error,
      "solve info not computed at matrix finalization time for requested transformation: " << trans);
    RCP<const cusparseHandle_t> hndl = CUSPARSEdetails::Session::getHandle();
    const Scalar s_alpha = Teuchos::ScalarTraits<Scalar>::one();
    cusparseStatus_t stat = CUSPARSEdetails::CUSPARSETemplateAdaptors<Scalar>::CSRSM_solve(
        *hndl, op, numMatRows, numRHS, &s_alpha,
        *matdescr_, rowVals_.getRawPtr(), rowPtrs_.getRawPtr(), colInds_.getRawPtr(), *solveInfo,
        data_y, stride_y, data_x, stride_x);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        stat != CUSPARSE_STATUS_SUCCESS,
        std::runtime_error, ": CSRSM_solve returned error " << stat);
    return;
  }


  template <class Scalar, class Node>
  template <class DomainScalar, class RangeScalar>
  void CUSPARSEOps<Scalar,Node>::multiply(Teuchos::ETransp trans,
                                          RangeScalar alpha,
                                          const MultiVector<DomainScalar,Node> &X,
                                                MultiVector< RangeScalar,Node> &Y) const
  {
    // CUSPARSE doesn't support mixed precision; partial specialize, then nix the generic versions
    Teuchos::CompileTimeAssert<Teuchos::TypeTraits::is_same<DomainScalar,Scalar>::value == false ||
                               Teuchos::TypeTraits::is_same< RangeScalar,Scalar>::value == false > cta; (void)cta;
    //
    std::string tfecfFuncName("multiply(trans,alpha,X,Y)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == false,
        std::runtime_error, "sparse operators have not been initialized with graph and matrix data; call setGraphAndMatrix() first.")
    // get pointers,stride from X and Y
    int stride_x = (int)X.getStride(),
        stride_y = (int)Y.getStride();
    const Scalar * data_x = X.getValues().getRawPtr();
    Scalar * data_y = Y.getValuesNonConst().getRawPtr();
    const int numMatRows = Y.getNumRows();
    const int numMatCols = X.getNumRows();
    const int numRHS     = X.getNumCols();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        X.getNumCols() != Y.getNumCols(),
        std::runtime_error, "X and Y do not have the same number of column vectors.")
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        Y.getNumRows() != (size_t)numRows_,
        std::runtime_error, "Y does not have the same number of rows as does the matrix.")
    // call mat-vec
    cusparseOperation_t op;
    if      (trans == Teuchos::NO_TRANS)   op = CUSPARSE_OPERATION_NON_TRANSPOSE;
    else if (trans == Teuchos::TRANS)      op = CUSPARSE_OPERATION_TRANSPOSE;
    else if (trans == Teuchos::CONJ_TRANS) op = CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE;
    RCP<const cusparseHandle_t> sess = CUSPARSEdetails::Session::getHandle();
    const Scalar s_alpha = (Scalar)alpha,
                 s_beta  = Teuchos::ScalarTraits<Scalar>::zero();
#ifdef HAVE_KOKKOS_DEBUG
    cudaError_t err = cudaGetLastError();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( cudaSuccess != err, std::runtime_error,
        ": cudaGetLastError() returned error after function call:\n"
        << cudaGetErrorString(err) );
    err = cudaThreadSynchronize();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( cudaSuccess != err, std::runtime_error,
        ": cudaThreadSynchronize() returned error after function call:\n"
        << cudaGetErrorString(err) );
#endif
    cusparseStatus_t stat;
    stat = CUSPARSEdetails::CUSPARSETemplateAdaptors<Scalar>::CSRMM(
        *sess, op, numMatRows, numRHS, numMatCols, numNZ_, &s_alpha,
        *matdescr_, rowVals_.getRawPtr(), rowPtrs_.getRawPtr(), colInds_.getRawPtr(),
        data_x, stride_x, &s_beta, data_y, stride_y);
#ifdef HAVE_KOKKOS_DEBUG
    err = cudaGetLastError();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( cudaSuccess != err, std::runtime_error,
        ": cudaGetLastError() returned error after function call:\n"
        << cudaGetErrorString(err) );
    err = cudaThreadSynchronize();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( cudaSuccess != err, std::runtime_error,
        ": cudaThreadSynchronize() returned error after function call:\n"
        << cudaGetErrorString(err) );
#endif
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        stat != CUSPARSE_STATUS_SUCCESS,
        std::runtime_error, ": CSRMM returned error " << stat);
    return;
  }


  template <class Scalar, class Node>
  template <class DomainScalar, class RangeScalar>
  void CUSPARSEOps<Scalar,Node>::multiply(Teuchos::ETransp trans,
                                          RangeScalar alpha, const MultiVector<DomainScalar,Node> &X,
                                          RangeScalar beta, MultiVector<RangeScalar,Node> &Y) const
  {
    // CUSPARSE doesn't support mixed precision; partial specialize, then nix the generic versions
    Teuchos::CompileTimeAssert<Teuchos::TypeTraits::is_same<DomainScalar,Scalar>::value == false ||
                               Teuchos::TypeTraits::is_same< RangeScalar,Scalar>::value == false > cta; (void)cta;
    //
    std::string tfecfFuncName("multiply(trans,alpha,X,beta,Y)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == false,
        std::runtime_error, "sparse operators have not been initialized with graph and matrix data; call setGraphAndMatrix() first.")
    // get pointers,stride from X and Y
    int stride_x = (int)X.getStride(),
        stride_y = (int)Y.getStride();
    const Scalar * data_x = X.getValues().getRawPtr();
    Scalar * data_y = Y.getValuesNonConst().getRawPtr();
    const int numMatRows = Y.getNumRows();
    const int numMatCols = X.getNumRows();
    const int numRHS     = X.getNumCols();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        X.getNumCols() != Y.getNumCols(),
        std::runtime_error, "X and Y do not have the same number of column vectors.")
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        Y.getNumRows() != (size_t)numRows_,
        std::runtime_error, "Y does not have the same number of rows as does the matrix.")
    // call mat-vec
    cusparseOperation_t op;
    RCP<const cusparseHandle_t> hndl = CUSPARSEdetails::Session::getHandle();
    if      (trans == Teuchos::NO_TRANS)   op = CUSPARSE_OPERATION_NON_TRANSPOSE;
    else if (trans == Teuchos::TRANS)      op = CUSPARSE_OPERATION_TRANSPOSE;
    else if (trans == Teuchos::CONJ_TRANS) op = CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE;
    const Scalar s_alpha = (Scalar)alpha,
                 s_beta  = (Scalar)beta;
    cusparseStatus_t stat = CUSPARSEdetails::CUSPARSETemplateAdaptors<Scalar>::CSRMM(
        *hndl, op, numMatRows, numRHS, numMatCols, numNZ_, &s_alpha,
        *matdescr_, rowVals_.getRawPtr(), rowPtrs_.getRawPtr(), colInds_.getRawPtr(),
        data_x, stride_x, &s_beta, data_y, stride_y);
#ifdef HAVE_KOKKOS_DEBUG
    cudaError_t err = cudaGetLastError();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( cudaSuccess != err, std::runtime_error,
        ": cudaGetLastError() returned error after function call:\n"
        << cudaGetErrorString(err) );
    err = cudaThreadSynchronize();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( cudaSuccess != err, std::runtime_error,
        ": cudaThreadSynchronize() returned error after function call:\n"
        << cudaGetErrorString(err) );
#endif
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        stat != CUSPARSE_STATUS_SUCCESS,
        std::runtime_error, ": CSRMM returned error " << stat);
    return;
  }



} // namespace Kokkos

#endif /* KOKKOS_CUSPARSEOPS_HPP */
