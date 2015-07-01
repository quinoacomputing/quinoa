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




#include "SundanceIQI_HdivLF_DivV_Cell.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_HDIV_TET_I1_FEM.hpp"

using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace Teuchos;

 
IQI_HdivLF_DivV_Cell::IQI_HdivLF_DivV_Cell( int spatialDim,
					    const CellType& maxCellType,
					    const BasisFamily& testBasis,
					    const QuadratureFamily& quad,
					    const ParameterList& verbParams):
  ElementIntegralLinearFormCell( spatialDim ,
				 maxCellType ,
				 testBasis ,
				 quad ,
				 verbParams ),
  DivV_( testBasis.nReferenceDOFs(maxCellType,maxCellType) , quad.getNumPoints(maxCellType) ),
  QP_( quad.getNumPoints( maxCellType ) , spatialDim ),
  QW_( quad.getNumPoints( maxCellType ) )
{
  // bypass testBasis.refEval and use Intrepid to fill DivV
  // warning: only works on tets right now.
  Intrepid::Basis_HDIV_TET_I1_FEM<double,Intrepid::FieldContainer<double> > myBasis;

  // move quadrature points into a field container
  Array<Point> qpSundance;
  Array<double> qwSundance;
  quad.getPoints( maxCellType , qpSundance , qwSundance );

  for (int i=0;i<(int)qpSundance.size();i++) {
    for (int j=0;j<3;j++) {
      QP_(i,j) = qpSundance[i][j];
    }
    QW_(i)=qwSundance[i];
  }

  // now tabulate the divergences 
  myBasis.getValues( DivV_ , QP_ , Intrepid::OPERATOR_DIV );

}

void IQI_HdivLF_DivV_Cell::evaluate( CellJacobianBatch& JTrans,
				     const double* const coeff,
				     RefCountPtr<Array<double> >& A) const
{
  const int nqp = quad().getNumPoints( maxCellType() );
  const int ncell = JTrans.numCells();
  const int nbf = testBasis().nReferenceDOFs(maxCellType(),maxCellType());

  // wrap A into a rank 2 field container
  Teuchos::Array<int> Aindices(2);
  Aindices[0] = nqp;
  Aindices[1] = nbf;
  Intrepid::FieldContainer<double> AFC(Aindices,*A);

  // wrap coeff into another field container.
  // by way of a Teuchos array
  // by way of discarding the const

  /* this surprisingly doesn't depend on the Jacobian ! */

  for (int c=0;c<ncell;c++) {
    for (int bf=0;bf<nbf;bf++) {
      AFC(c,bf) = 0.0;
      for (int q=0;q<nqp;q++) {
	AFC(c,bf) += QW_(q) * coeff[c*nqp+q] * DivV_(bf,q );
      }
    }
  }

  return;
}
