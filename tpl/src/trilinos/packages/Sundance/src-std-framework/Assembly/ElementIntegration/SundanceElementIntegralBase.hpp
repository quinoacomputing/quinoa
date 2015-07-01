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

#ifndef SUNDANCE_ELEMENTINTEGRALBASE_H
#define SUNDANCE_ELEMENTINTEGRALBASE_H

#include "SundanceDefs.hpp"
#include "SundanceCellJacobianBatch.hpp"
#include "SundanceQuadratureFamily.hpp"
#include "SundanceBasisFamily.hpp"
#include "Teuchos_Array.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;
  
  namespace Internal
  {
    using namespace Teuchos;
    
    /** 
     * ElementIntegralFunctional provides an abstract base class for storing
     * information about the geometry/topology over which the form is defined
     * and the quadrature rule.
     * Computational routines are provided by subclasses; this is for storage
     * and information only, providing no particular interface to integration.
    */
    class ElementIntegralBase
      : public TSFExtended::ParameterControlledObjectWithVerbosity<ElementIntegralBase>,
	public TSFExtended::Printable
    {
    public:
      /** Constructor */
      ElementIntegralBase( int spatialDim,
			   const CellType& maxCellType,
			   int dim, 
			   const CellType& cellType,
			   const QuadratureFamily& quad,
			   const ParameterList& verbParams 
			   = *ElementIntegralBase::defaultVerbParams()):
	ParameterControlledObjectWithVerbosity<ElementIntegralBase>("Integration",verbParams),
	spatialDim_( spatialDim ), 
	maxCellType_( maxCellType ) ,
	dim_( dim ),
	cellType_( cellType ),
	quad_( quad ),
	verbParams_( verbParams ) {;}
	
      
      /** Destructor */
      virtual ~ElementIntegralBase() {;}
      
      virtual int spatialDim() const { return spatialDim_; }
      virtual int dim() const { return dim_; }
      virtual const CellType & maxCellType() const { return maxCellType_; }
      virtual const CellType & cellType() const { return cellType_; }
      virtual const QuadratureFamily & quad() const { return quad_; }
      virtual const ParameterList &verbParams() const { return verbParams_; }
 
      static RefCountPtr<ParameterList> defaultVerbParams()
      {
	static RefCountPtr<ParameterList> rtn = rcp(new ParameterList("Integration"));
	static int first = true;
	if (first)
	  {
	    rtn->set<int>("setup", 0);
	    rtn->set<int>("transformation", 0);
	    rtn->set<int>("integration", 0);
	    rtn->set<int>("extract weak form", 0);
	    rtn->set<int>("find groups", 0);
	    first = false;
	  }
	return rtn;
      }
      
    private:
      const int spatialDim_;
      const CellType &maxCellType_;
      const int dim_;
      const CellType &cellType_;
      const QuadratureFamily & quad_;
      const ParameterList & verbParams_;
      
    };
  }
}
#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
