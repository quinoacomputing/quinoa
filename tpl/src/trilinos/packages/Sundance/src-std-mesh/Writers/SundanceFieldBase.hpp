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

#ifndef SUNDANCE_FIELDBASE_H
#define SUNDANCE_FIELDBASE_H

#include "SundanceDefs.hpp"
#include "PlayaHandleable.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceCellSet.hpp"



namespace Sundance
{
class CellFilter;
}


namespace Sundance
{
class Mesh;
using Sundance::CellFilter;

using Sundance::Map;
using Teuchos::Array;
using Teuchos::RefCountPtr;
/**
 *
 */
class FieldBase : public Playa::Handleable<FieldBase>
{
public:
  /** */
  FieldBase(){;}

  /** virtual dtor */
  virtual ~FieldBase(){;}

  /** */
  virtual int numElems() const {return 1;}

  /** */
  virtual double getData(int cellDim, int cellID, int elem) const = 0 ;

  /** */
  virtual bool isDefined(int cellDim, int cellID, int elem) const = 0 ;

  /** */
  virtual bool isPointData() const = 0 ;

  /** */
  virtual bool isCellData() const {return !isPointData();}

  /** Get a batch of data. 
   * \param batch Output array of data values. This is a 2D array packed
   * into a 1D vector with function index as the faster running index.
   */
  virtual void getDataBatch(int cellDim, const Array<int>& cellID,
    const Array<int>& funcElem, Array<double>& batch) const ;

  /**
   * Return the cell filter on which this field is defined 
   */
  virtual const CellFilter& domain() const ;

};

using Sundance::CellSet;
using Sundance::CellFilter;

CellSet connectedNodeSet(const CellFilter& f, const Mesh& mesh);
RCP<Array<int> > cellSetToLIDArray(const CellSet& cs);
}

#endif
