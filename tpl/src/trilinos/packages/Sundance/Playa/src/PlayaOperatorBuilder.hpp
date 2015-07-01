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


#ifndef PLAYA_OPERATORBUILDER_HPP
#define PLAYA_OPERATORBUILDER_HPP

#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaLinearCombinationDecl.hpp"
#include "PlayaVectorType.hpp"
#include "PlayaVectorSpaceDecl.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_GlobalMPISession.hpp"


using namespace Playa;
using namespace Teuchos;


namespace Playa
{
  /** Base class for building test operators */
  template <class Scalar>
  class OperatorBuilder
  {
  public:
    /** */
    OperatorBuilder(int nLocal, const VectorType<Scalar>& vecType);
    /** */
    OperatorBuilder(int nLocalDomain, int nLocalRange,
                    const VectorType<Scalar>& vecType);
    /** */
    OperatorBuilder(const VectorSpace<Scalar>& domain,
                    const VectorSpace<Scalar>& range,
                    const VectorType<Scalar>& vecType);
    /** */
    virtual ~OperatorBuilder(){;}

    /** */
    const VectorType<Scalar>& vecType() const {return vecType_;}

    /** */
    const VectorSpace<Scalar>& domain() const {return domain_;}

    /** */
    const VectorSpace<Scalar>& range() const {return range_;}

    /** */
    virtual LinearOperator<Scalar> getOp() const = 0 ; 

  protected:

  private:
    VectorType<Scalar> vecType_;

    VectorSpace<Scalar> domain_;

    VectorSpace<Scalar> range_;
  };

  template <class Scalar> 
  inline OperatorBuilder<Scalar>
  ::OperatorBuilder(int nLocalRows, const VectorType<Scalar>& vecType)
    : vecType_(vecType), domain_(), range_()
  {
    range_ = vecType_.createEvenlyPartitionedSpace(MPIComm::world(), nLocalRows);
    domain_ = range_;
  }

  template <class Scalar> 
  inline OperatorBuilder<Scalar>
  ::OperatorBuilder(int nLocalDomain, int nLocalRange,
                    const VectorType<Scalar>& vecType)
    : vecType_(vecType), domain_(), range_()
  {
    range_ = vecType_.createEvenlyPartitionedSpace(MPIComm::world(), nLocalRange);
    domain_ = vecType_.createEvenlyPartitionedSpace(MPIComm::world(), nLocalDomain);
  }

  

  template <class Scalar> 
  inline OperatorBuilder<Scalar>
  ::OperatorBuilder(const VectorSpace<Scalar>& domain,
                    const VectorSpace<Scalar>& range,
                    const VectorType<Scalar>& vecType)
    : vecType_(vecType), domain_(domain), range_(range)
  {}

}

#endif
