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

#ifndef SUNDANCE_IQI_LF_VN_FACET_H
#define SUNDANCE_IQI_LF_VN_FACET_H

#include "SundanceDefs.hpp"
#include "SundanceCellJacobianBatch.hpp"
#include "SundanceQuadratureFamily.hpp"
#include "SundanceBasisFamily.hpp"
#include "Teuchos_Array.hpp"
#include "SundanceElementIntegralLinearFormFacet.hpp"
#include "Intrepid_FieldContainer.hpp"

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
     * 
     IQI: IntrepidQuadratureIntegral
     HdivLF: LinearForm where test function lies in H(div)
     VN: the test function has its normal component taken
     Facet: this only makes sense on codimension 1 facets.

     This class evaluates integrals of the form
     \int_e f v.n
     where e is a facet of co-dimension 1, f is tabulated
     at quadrature points, and v is a test function in an H(div) space.
    */
    class IQI_HdivLF_VN_Facet : public ElementIntegralLinearFormFacet
    {
    public:
      /** Constructor */
      // this is only going to take maxCellType
      // as an argument, then use SundanceStdMesh::facetType(maxCellType,spatialDim-1,0)
      // to pass up the chain since we can only work on facets of codimension 1
      IQI_HdivLF_VN_Facet( int spatialDim,
			   const CellType& maxCellType,
			   const BasisFamily& testBasis,
			   const QuadratureFamily& quad,
			   const ParameterList& verbParams 
			   = *ElementIntegralBase::defaultVerbParams());

      /** Destructor */
      virtual ~IQI_HdivLF_VN_Facet() {;}
      
      /** evaluates integral of coeff against divergence of each basis function */
      virtual void evaluate( CellJacobianBatch& JTrans,
			     CellJacobianBatch& JVol ,
			     const Array<int>& facetIndex ,
			     const double* const coeff,
			     RefCountPtr<Array<double> >& A) const;


    };
  }
}
#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
