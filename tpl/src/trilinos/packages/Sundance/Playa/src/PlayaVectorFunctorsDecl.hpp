/* @HEADER@ */
// ************************************************************************
// 
//                 Playa: Programmable Linear Algebra
//                 Copyright 2012 Sandia Corporation
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
// Questions? Contact Kevin Long (kevin.long@ttu.edu)
// 

/* @HEADER@ */

#ifndef PLAYA_VECTORFUNCTORDECL_HPP
#define PLAYA_VECTORFUNCTORDECL_HPP

#include "PlayaDefs.hpp"
#include "PlayaMPIComm.hpp"
#include "Teuchos_RCP.hpp"
#include "PlayaGeneralizedIndex.hpp"

namespace PlayaFunctors
{

using Playa::GeneralizedIndex;
using Teuchos::RCP;
using Playa::MPIComm;
using Playa::MPIOp;
using Playa::MPIDataType;

/**
 * \brief This traits class specifies the return type of a reduction functor. 
 * If not specialized, the default return type will be a Scalar.
 *
 * @author Kevin Long (kevin.long@ttu.edu)
 */
template <class Scalar, class FunctorType>
class VectorFunctorTraits
{
public:
  typedef Scalar ReturnType;
};



/** 
 * \brief Base class for reduction functors
 *
 * @author Kevin Long (kevin.long@ttu.edu)
*/
template <class Scalar>
class ReductionFunctorBase
{
public:
  /** Construct with a communicator */
  ReductionFunctorBase(
    const MPIComm& comm
    )
    : comm_(comm) {}

  /** Callback for any postprocessing step (for example, MPI all-reduce) */
  virtual void postProc() const = 0 ;

protected:

  /** Return the MPI communicator */
  const MPIComm& comm() const {return comm_;}

private:
  MPIComm comm_;
};

  
/**
 * \brief IndexedValue is the return type for reduction operations such
 * as MinLoc that return a location and a value. 
 */
template <class Scalar>
struct IndexedValue
{
  /** Value */
  Scalar what;
  /** Index */
  int where;
};
  
}

#endif
