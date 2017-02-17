// @HEADER
// ***********************************************************************
//
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov)
//
// ***********************************************************************
// @HEADER

#include "Thyra_SpmdVectorSpaceUtilities.hpp"
#include "Teuchos_CommHelpers.hpp"


namespace Thyra {


Ordinal SpmdVectorSpaceUtilities::computeMapCode(
  const Teuchos::Comm<Ordinal> &comm, const Ordinal localSubDim
  )
{
  using Teuchos::outArg;
  using Teuchos::REDUCE_SUM;
  using Teuchos::reduceAll;
  //
  // Here we will make a map code out of just the local sub-dimension on each
  // processor.  If each processor has the same number of local elements, then
  // the map codes will be the same and this is all you need for RTOp
  // compatibility.
  //
  const int procRank = comm.getSize ();
  Ordinal mapCode = -1;
  Ordinal localCode = localSubDim % (procRank+1) + localSubDim;
  reduceAll<Ordinal, Ordinal> (comm, REDUCE_SUM, localCode, outArg (mapCode));
  return mapCode;
}


Ordinal SpmdVectorSpaceUtilities::computeLocalOffset(
  const Teuchos::Comm<Ordinal> &comm, const Ordinal localSubDim
  )
{
  using Teuchos::outArg;
  using Teuchos::REDUCE_SUM;
  using Teuchos::scan;

  Ordinal localOffset;
  const Ordinal _localOffset = localSubDim;
  scan<Ordinal, Ordinal> (comm, REDUCE_SUM, _localOffset, outArg (localOffset));
  localOffset -= localSubDim;
  return localOffset;
}


Ordinal SpmdVectorSpaceUtilities::computeGlobalDim(
  const Teuchos::Comm<Ordinal> &comm, const Ordinal localSubDim
  )
{
  using Teuchos::outArg;
  using Teuchos::REDUCE_SUM;
  using Teuchos::reduceAll;

  Ordinal globalDim = -1;
  reduceAll<Ordinal, Ordinal> (comm, REDUCE_SUM, localSubDim, outArg (globalDim));
  return globalDim;
}


} // namespace Thyra
