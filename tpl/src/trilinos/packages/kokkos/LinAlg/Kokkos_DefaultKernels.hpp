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

#ifndef KOKKOS_DEFAULT_KERNELS_
#define KOKKOS_DEFAULT_KERNELS_

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultSparseOps.hpp"
#include "Kokkos_DefaultBlockSparseOps.hpp"
#include "Kokkos_DefaultRelaxation.hpp"
#ifdef HAVE_KOKKOS_CUSPARSE
#include "Kokkos_CUSPARSEOps.hpp"
#endif
// #include "Kokkos_FirstTouchSparseOps.hpp"

namespace Kokkos {

  /** \brief Traits class providing default kernel types for CRS, block CRS and relaxation kernels.
      \ingroup kokkos_crs_ops
   */
  template <class Scalar, class Ordinal, class Node>
  struct DefaultKernels {
    typedef DefaultHostSparseOps <void  ,Ordinal,Node>  SparseOps;
    typedef DefaultBlockSparseOps<Scalar,Ordinal,Node>  BlockSparseOps;
    typedef DefaultRelaxation    <Scalar,Ordinal,Node>  Relaxations;
  };

// #ifndef HAVE_KOKKOS_NO_FIRST_TOUCH_MATVEC_ALLOCATION
//   class TBBNode;
//   template <class Scalar, class Ordinal>
//   struct DefaultKernels<Scalar,Ordinal,TBBNode> {
//     typedef FirstTouchSparseOps  <void  ,Ordinal,TBBNode>  SparseOps;
//     typedef DefaultBlockSparseOps<Scalar,Ordinal,TBBNode>  BlockSparseOps;
//     typedef DefaultRelaxation    <Scalar,Ordinal,TBBNode>  Relaxations;
//   };
//   class TPINode;
//   template <class Scalar, class Ordinal>
//   struct DefaultKernels<Scalar,Ordinal,TPINode> {
//     typedef FirstTouchSparseOps  <void  ,Ordinal,TPINode>  SparseOps;
//     typedef DefaultBlockSparseOps<Scalar,Ordinal,TPINode>  BlockSparseOps;
//     typedef DefaultRelaxation    <Scalar,Ordinal,TPINode>  Relaxations;
//   };
// #endif

  /** \brief Traits class providing default kernel types for CRS, block CRS and relaxation kernels.
      \ingroup kokkos_crs_ops
    
      For ThrustGPUNode, defaults are the same as in general, except that the default sparse ops should be provided by 
      DefaultDeviceSparseOps.
   */
   class ThrustGPUNode;
   template <class Scalar, class Ordinal>
   struct DefaultKernels<Scalar,Ordinal,ThrustGPUNode> {
     // empty == fail
   };
#ifdef HAVE_KOKKOS_CUSPARSE
   template <class Scalar>
   struct DefaultKernels<Scalar,int,ThrustGPUNode> {
     typedef CUSPARSEOps<void,ThrustGPUNode> SparseOps;
   };
#endif

}

#endif
