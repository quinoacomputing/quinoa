/*
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
*/
#ifndef KOKKOS_BLAS1_MV_IMPL_NORMW_HPP_
#define KOKKOS_BLAS1_MV_IMPL_NORMW_HPP_

#include <TpetraKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>

namespace KokkosBlas {
namespace Impl {

/// \brief Weighted 2-norm functor for single vectors.
///
/// \tparam RV 1-D output View
/// \tparam XMV 2-D input View
/// \tparam SizeType Index type.  Use int (32 bits) if possible.
template<class RV, class XMV, class SizeType = typename XMV::size_type>
struct MV_Nrm2w_Functor
{
  typedef typename XMV::execution_space             execution_space;
  typedef SizeType                                        size_type;
  typedef typename XMV::non_const_value_type            xvalue_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<xvalue_type> IPT;
  typedef typename IPT::mag_type                           mag_type;
  typedef mag_type                                     value_type[];

  const size_type value_count;
  RV r_;
  XMV X_;
  XMV W_;

  MV_Nrm2w_Functor (const RV& r, const XMV& X, const XMV& W) :
    value_count (X.dimension_1 ()), r_ (r), X_ (X), W_ (W)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i, value_type sums) const {
    // NOTE (mfh 09 Apr 2015) Don't bother optimizing this with e.g.,
    // vector pragmas.  I only left this function in place for
    // backwards compatibility.  If you want weird functions like this
    // on MultiVector, and you want them to be super fast, you should
    // write your own Kokkos kernels for them.  Then you can load them
    // up with all the pragmas you like.
    for (size_type j = 0; j < value_count; ++j) {
      const mag_type tmp_X = IPT::norm (X_(i,j));
      const mag_type tmp_W = IPT::norm (W_(i,j));
      const mag_type quotient = tmp_X / tmp_W;

      sums[j] += quotient * quotient;
    }
  }

  KOKKOS_INLINE_FUNCTION void init (value_type update) const
  {
    for (size_type j = 0; j < value_count; ++j) {
      update[j] = Kokkos::Details::ArithTraits<mag_type>::zero ();
    }
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type update,
        const volatile value_type source) const
  {
    for (size_type j = 0; j < value_count; ++j) {
      update[j] += source[j];
    }
  }

  // On device, write the reduction result to the output View.
  KOKKOS_INLINE_FUNCTION void final (const value_type dst) const
  {
    for (size_type j = 0; j < value_count; ++j) {
      r_(j) = dst[j];
    }
  }
};

/// \brief Weighted 2-norm functor for single vectors.
///
/// "Weighted" means X(i) / W(i).
///
/// \tparam RV 0-D output View
/// \tparam XV 1-D input View
/// \tparam SizeType Index type.  Use int (32 bits) if possible.
template<class RV, class XV, class SizeType = typename XV::size_type>
struct V_Nrm2w_Functor
{
  typedef typename XV::execution_space              execution_space;
  typedef SizeType                                        size_type;
  typedef typename XV::non_const_value_type             xvalue_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<xvalue_type> IPT;
  typedef typename IPT::mag_type                           mag_type;
  typedef mag_type                                       value_type;

  RV r_;
  XV X_;
  XV W_;

  V_Nrm2w_Functor (const RV& r, const XV& X, const XV& W) :
    r_ (r), X_ (X), W_ (W)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i, value_type& sum) const {
    const mag_type tmp_X = IPT::norm (X_(i));
    const mag_type tmp_W = IPT::norm (W_(i));
    const mag_type quotient = tmp_X / tmp_W;

    sum += quotient * quotient;
  }

  KOKKOS_INLINE_FUNCTION void init (value_type& update) const
  {
    update = Kokkos::Details::ArithTraits<mag_type>::zero ();
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type& update,
        const volatile value_type& source) const
  {
    update += source;
  }

  // On device, write the reduction result to the output View.
  KOKKOS_INLINE_FUNCTION void final (const value_type& dst) const
  {
    r_() = dst;
  }
};

/// \brief Implementation of KokkosBlas::nrm2w_squared for
///   multivectors and single vectors.
template<class RV, class XMV, int rank = XMV::rank>
struct Nrm2w {};

//! Specialization for multivectors.
template<class RV, class XMV>
struct Nrm2w<RV, XMV, 2> {
  typedef typename XMV::execution_space execution_space;
  typedef typename XMV::size_type size_type;

  static void nrm2w_squared (const RV& r, const XMV& X, const XMV& W)
  {
    const size_type numRows = X.dimension_0 ();
    const size_type numCols = X.dimension_1 ();

    // int is generally faster than size_t, but check for overflow first.
    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      typedef MV_Nrm2w_Functor<RV, XMV, int> functor_type;
      Kokkos::RangePolicy<execution_space, int> policy (0, numRows);
      functor_type op (r, X, W);
      Kokkos::parallel_reduce (policy, op);
    }
    else {
      typedef MV_Nrm2w_Functor<RV, XMV, size_type> functor_type;
      Kokkos::RangePolicy<execution_space, size_type> policy (0, numRows);
      functor_type op (r, X, W);
      Kokkos::parallel_reduce (policy, op);
    }
  }
};

//! Specialization for single vectors.
template<class R, class XV>
struct Nrm2w<R, XV, 1> {
  typedef typename XV::execution_space execution_space;
  typedef typename XV::size_type size_type;

  static void nrm2w_squared (const R& r, const XV& X, const XV& W)
  {
    const size_type numRows = X.dimension_0 ();

    // int is generally faster than size_t, but check for overflow first.
    if (numRows < static_cast<size_type> (INT_MAX)) {
      typedef V_Nrm2w_Functor<R, XV, int> functor_type;
      Kokkos::RangePolicy<execution_space, int> policy (0, numRows);
      functor_type op (r, X, W);
      Kokkos::parallel_reduce (policy, op);
    }
    else {
      typedef V_Nrm2w_Functor<R, XV, size_type> functor_type;
      Kokkos::RangePolicy<execution_space, size_type> policy (0, numRows);
      functor_type op (r, X, W);
      Kokkos::parallel_reduce (policy, op);
    }
  }
};

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Nrm2w for rank == 2.  This is NOT for users!!!
// All the declarations of full specializations go in this header
// file.  We may spread out definitions (see _DEF macro below) across
// one or more .cpp files.
//

#define KOKKOSBLAS_IMPL_MV_NRM2W_RANK2_DECL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
template<> \
struct Nrm2w<Kokkos::View<Kokkos::Details::InnerProductSpaceTraits< SCALAR >::mag_type*, \
                          EXEC_SPACE::array_layout, \
                          Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
             Kokkos::View<const SCALAR**, \
                          LAYOUT, \
                          Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
             2 > \
{ \
  typedef Kokkos::View<Kokkos::Details::InnerProductSpaceTraits< SCALAR >::mag_type*, \
                       EXEC_SPACE::array_layout, \
                       Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const SCALAR**, \
                       LAYOUT, \
                       Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XMV; \
  typedef XMV::execution_space execution_space; \
  typedef XMV::size_type size_type; \
 \
  static void nrm2w_squared (const RV& r, const XMV& X, const XMV& W); \
};

//
// Declarations of full specializations of Impl::Nrm2w for rank == 2.
// Their definitions go in .cpp file(s) in this source directory.
//

#ifdef TPETRAKERNELS_BUILD_EXECUTION_SPACE_SERIAL

KOKKOSBLAS_IMPL_MV_NRM2W_RANK2_DECL( double, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace )

#endif // TPETRAKERNELS_BUILD_EXECUTION_SPACE_SERIAL

#ifdef TPETRAKERNELS_BUILD_EXECUTION_SPACE_OPENMP

KOKKOSBLAS_IMPL_MV_NRM2W_RANK2_DECL( double, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace )

#endif // TPETRAKERNELS_BUILD_EXECUTION_SPACE_OPENMP

#ifdef TPETRAKERNELS_BUILD_EXECUTION_SPACE_PTHREAD

KOKKOSBLAS_IMPL_MV_NRM2W_RANK2_DECL( double, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace )

#endif // TPETRAKERNELS_BUILD_EXECUTION_SPACE_PTHREAD

#ifdef TPETRAKERNELS_BUILD_EXECUTION_SPACE_CUDA

KOKKOSBLAS_IMPL_MV_NRM2W_RANK2_DECL( double, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace )

#endif // TPETRAKERNELS_BUILD_EXECUTION_SPACE_CUDA

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Nrm2w for rank == 2.  This is NOT for users!!!
//

#define KOKKOSBLAS_IMPL_MV_NRM2W_RANK2_DEF( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
void \
Nrm2w<Kokkos::View<Kokkos::Details::InnerProductSpaceTraits< SCALAR >::mag_type*, \
                   EXEC_SPACE::array_layout, \
                   Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<const SCALAR**, \
                   LAYOUT, \
                   Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      2 >:: \
nrm2w_squared (const RV& r, const XMV& X, const XMV& W) \
{ \
  const size_type numRows = X.dimension_0 (); \
  const size_type numCols = X.dimension_1 (); \
 \
  if (numRows < static_cast<size_type> (INT_MAX) && \
      numRows * numCols < static_cast<size_type> (INT_MAX)) { \
    typedef MV_Nrm2w_Functor<RV, XMV, int> functor_type; \
    Kokkos::RangePolicy<execution_space, int> policy (0, numRows); \
    functor_type op (r, X, W); \
    Kokkos::parallel_reduce (policy, op); \
  } \
  else { \
    typedef MV_Nrm2w_Functor<RV, XMV, size_type> functor_type; \
    Kokkos::RangePolicy<execution_space, size_type> policy (0, numRows); \
    functor_type op (r, X, W); \
    Kokkos::parallel_reduce (policy, op); \
  } \
}



} // namespace Impl
} // namespace KokkosBlas

#endif // KOKKOS_BLAS1_MV_IMPL_NORMW_HPP_
