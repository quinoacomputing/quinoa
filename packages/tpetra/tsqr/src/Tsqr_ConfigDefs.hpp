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

#ifndef TSQR_CONFIGDEFS_HPP
#define TSQR_CONFIGDEFS_HPP

/// \file Tsqr_ConfigDefs.hpp
/// \brief Include this header to get TSQR's configuration options.

// Users should not include TpetraTSQR_config.h directly.
// Include Tsqr_ConfigDefs.hpp instead.
#include <TpetraTSQR_config.h>

/// \namespace TSQR
/// \brief Implementation of the Tall Skinny QR (TSQR) factorization.
///
/// This namespace contains a full hybrid-parallel (MPI + Kokkos)
/// implementation of the Tall Skinny QR (TSQR) factorization.  The
/// following paper describes the implementation:
///
/// Mark Hoemmen.  "A communication-avoiding, hybrid-parallel,
/// rank-revealing orthogonalization method."  IEEE International
/// Parallel and Distributed Processing Symposium (IPDPS), April 2011.
///
/// For further details, see the following:
///
/// Marghoob Mohiyuddin, Mark Hoemmen, James Demmel, and Kathy Yelick.
/// "Minimizing Communication in Sparse Matrix Solvers."  In
/// Proceedings of Supercomputing 2009, November 2009.
///
/// James Demmel, Laura Grigori, Mark Hoemmen, and Julien Langou.
/// "Communication-optimal parallel and sequential QR and LU
/// factorizations."  SIAM Journal on Scientific Computing, Volume 34,
/// Issue 1, 2012.
namespace TSQR {
  //
  // We declare the TSQR namespace here so that Doxygen will find it
  // and pull in all its documentation.
  //

  /// \namespace TSQR::Test
  /// \brief Accuracy and performance tests for TSQR.
  ///
  /// The classes and routines here are not intended for consumers of
  /// TSQR, but may be helpful as examples.
  namespace Test {
    //
    // We declare the TSQR::Test namespace here so that Doxygen will
    // find it and pull in all its documentation.
    //
  } // namespace Test
} // namespace TSQR

#endif // TSQR_CONFIGDEFS_HPP
