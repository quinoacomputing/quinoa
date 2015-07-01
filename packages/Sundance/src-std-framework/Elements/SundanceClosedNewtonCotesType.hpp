/* @HEADER@ */
// ************************************************************************
// 
//                             Sundance
//                 Copyright 2011 Sandia Corporation
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

#ifndef SUNDANCE_CLOSEDNEWTONCOTESTYPE_H
#define SUNDANCE_CLOSEDNEWTONCOTESTYPE_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceQuadratureFamily.hpp"
#include "SundanceQuadratureTypeBase.hpp"

namespace Sundance
{
using namespace Teuchos;
  

/** 
 * Generates closed Newton-Cotes quadrature rules
 */
class ClosedNewtonCotesType : public QuadratureTypeBase
{
public:
  /** */
  ClosedNewtonCotesType();

  /** */
  virtual ~ClosedNewtonCotesType(){;}

  /** */
  virtual XMLObject toXML() const ;

  /** Indicate whether the given cell type is supported at any order */
  virtual bool supportsCellType(const CellType& cellType) const ;
    
  /** Indicate whether the given cell type is supported at the
   * specified order */
  virtual bool supports(const CellType& cellType, int order) const ;
    
  /** Return the max quadrature order available on the given cell type */
  virtual int maxOrder(const CellType& cellType) const ;

  /** Create a quadrature family of the specified order */
  virtual QuadratureFamily createQuadFamily(int order) const ;

  /** Indicate whether there is a maximum order for quadrature rules
   * available on the given cell type. */
  virtual bool hasLimitedOrder(const CellType& cellType) const {return true;}

  /* handleable boilerplate */
  GET_RCP(QuadratureTypeBase);

};
}


#endif
