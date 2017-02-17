/*
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
*/

#include "Tpetra_ConfigDefs.hpp"

#if defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION) && defined(HAVE_TPETRA_QD)

#include <qd/dd_real.h>

// definitions
#include "Tpetra_Map_decl.hpp"
#include "Tpetra_Map_def.hpp"
#include "Tpetra_MultiVector_decl.hpp"
#include "Tpetra_MultiVector_def.hpp"
#include "Tpetra_Vector_decl.hpp"
#include "Tpetra_Vector_def.hpp"

// nodes
#include <Kokkos_SerialNode.hpp>
#if defined(HAVE_KOKKOS_TBB)
#  include <Kokkos_TBBNode.hpp>
#endif
#if defined(HAVE_KOKKOS_THREADPOOL)
#  include <Kokkos_TPINode.hpp>
#endif
#if defined(HAVE_KOKKOS_OPENMP)
#  include <Kokkos_OpenMPNode.hpp>
#endif

/* these examples require some explicit instantiations that is not enabled in the build of the library 
*/

#define INSTANT_ALL( NODE ) \
  TPETRA_MULTIVECTOR_INSTANT(dd_real,int,int,NODE) \
  TPETRA_VECTOR_INSTANT(     dd_real,int,int,NODE)

namespace Tpetra {

  INSTANT_ALL(Kokkos::SerialNode)
#if defined(HAVE_KOKKOS_TBB)
  INSTANT_ALL(Kokkos::TBBNode)
#endif
#if defined(HAVE_KOKKOS_THREADPOOL)
  INSTANT_ALL(Kokkos::TPINode)
#endif
#if defined(HAVE_KOKKOS_OPENMP)
  INSTANT_ALL(Kokkos::OpenMPNode)
#endif

}

#endif

