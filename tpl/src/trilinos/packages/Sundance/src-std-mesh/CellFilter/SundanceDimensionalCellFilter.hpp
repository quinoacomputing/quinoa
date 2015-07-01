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

#ifndef SUNDANCE_DIMENSIONALCELLFILTER_H
#define SUNDANCE_DIMENSIONALCELLFILTER_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceCellFilterBase.hpp"


namespace Sundance
{
using namespace Teuchos;
  


/**
 * DimensionalCellFilter is a filter that identifies all cells of a 
 * specified dimension. 
 *
 * <h4> Example: </h4> get all faces in a 3D mesh
 *
 * \code
 * Mesh myMesh = myReader.getMesh();
 * CellFilter faceFilter = new DimensionalCellFilter(2);
 * CellSet faces = faceFilter.getCells(myMesh);
 * \endcode
 */
class DimensionalCellFilter : public CellFilterBase 
{
public:
  /** */
  DimensionalCellFilter(int dim);

  /** */
  virtual ~DimensionalCellFilter(){;}

  /** */
  virtual int dimension(const Mesh& mesh) const {return dim_;}

  /** */
  virtual XMLObject toXML() const ;

  /** */
  virtual std::string typeName() const {return "DimensionalCellFilter";}

  /** */
  virtual std::string description() const 
    {return "Cells(d=" + Teuchos::toString(dim_) + ")";}

  /** */
  virtual bool lessThan(const CellFilterStub* other) const ;

  /* */
  GET_RCP(CellFilterStub);

protected:
  /** get the cells */
  virtual CellSet internalGetCells(const Mesh& mesh) const ;
    
private:
    
  int dim_;


};

}



#endif
