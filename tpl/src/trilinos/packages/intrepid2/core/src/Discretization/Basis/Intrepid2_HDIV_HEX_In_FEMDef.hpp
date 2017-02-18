// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_HDIV_HEX_In_FEMDef.hpp
    \brief  Definition file for FEM basis functions of degree n for H(div) functions on HEX cells.
    \author Created by R. Kirby, P. Bochev, D. Ridzal and K. Peterson.
*/

namespace Intrepid2 {

  template<class Scalar, class ArrayScalar>
  Basis_HDIV_HEX_In_FEM<Scalar,ArrayScalar>::Basis_HDIV_HEX_In_FEM( int order ,
								      const ArrayScalar & ptsClosed ,
								      const ArrayScalar & ptsOpen):
    closedBasis_( order , ptsClosed ),
    openBasis_( order-1 , ptsOpen ),
    closedPts_( ptsClosed ),
    openPts_( ptsOpen )
  {
    this -> basisDegree_       = order;
    this -> basisCardinality_  = 3 * closedBasis_.getCardinality() * openBasis_.getCardinality() * openBasis_.getCardinality();
    this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >() );
    this -> basisType_         = BASIS_FEM_FIAT;
    this -> basisCoordinates_  = COORDINATES_CARTESIAN;
    this -> basisTagsAreSet_   = false;

    Array<Array<RCP<Basis<Scalar,ArrayScalar > > > > bases(3);
    bases[0].resize(3); bases[1].resize(3); bases[2].resize(3);
    bases[0][0] = rcp( &closedBasis_ , false );
    bases[0][1] = rcp( &openBasis_ , false );
    bases[0][2] = rcp( &openBasis_ , false );
    bases[1][0] = rcp( &openBasis_ , false );
    bases[1][1] = rcp( &closedBasis_ , false );
    bases[1][2] = rcp( &openBasis_ , false );
    bases[2][0] = rcp( &openBasis_ , false );
    bases[2][1] = rcp( &openBasis_ , false );
    bases[2][2] = rcp( &closedBasis_ , false );
    this->setBases( bases );

    initializeTags();
    this->basisTagsAreSet_ = true;
  }

  template<class Scalar, class ArrayScalar>
  Basis_HDIV_HEX_In_FEM<Scalar,ArrayScalar>::Basis_HDIV_HEX_In_FEM( int order , const EPointType &pointType ):
    closedBasis_( order , pointType==POINTTYPE_SPECTRAL?POINTTYPE_SPECTRAL:POINTTYPE_EQUISPACED ),
    openBasis_( order-1 , pointType==POINTTYPE_SPECTRAL?POINTTYPE_SPECTRAL_OPEN:POINTTYPE_EQUISPACED ),
    closedPts_( order+1 , 1 ),
    openPts_( order , 1 )
  {
    this -> basisDegree_       = order;
    this -> basisCardinality_  = 3 * closedBasis_.getCardinality() * openBasis_.getCardinality() * openBasis_.getCardinality();
    this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >() );
    this -> basisType_         = BASIS_FEM_FIAT;
    this -> basisCoordinates_  = COORDINATES_CARTESIAN;
    this -> basisTagsAreSet_   = false;

    PointTools::getLattice<Scalar,ArrayScalar >( closedPts_ ,
                                                            shards::CellTopology(shards::getCellTopologyData<shards::Line<2> >()) ,
                                                            order ,
                                                            0 ,
                                                            pointType==POINTTYPE_SPECTRAL?POINTTYPE_WARPBLEND:POINTTYPE_EQUISPACED );

    if (pointType == POINTTYPE_SPECTRAL)
      {
	PointTools::getGaussPoints<Scalar,ArrayScalar >( openPts_ ,
								    order - 1 );
      }
    else
      {
	PointTools::getLattice<Scalar,ArrayScalar >( openPts_ ,
								shards::CellTopology(shards::getCellTopologyData<shards::Line<2> >()) ,
								order - 1,
								0 ,
								POINTTYPE_EQUISPACED );

      }


    Array<Array<RCP<Basis<Scalar,ArrayScalar > > > > bases(3);
    bases[0].resize(3); bases[1].resize(3); bases[2].resize(3);
    bases[0][0] = rcp( &closedBasis_ , false );
    bases[0][1] = rcp( &openBasis_ , false );
    bases[0][2] = rcp( &openBasis_ , false );
    bases[1][0] = rcp( &openBasis_ , false );
    bases[1][1] = rcp( &closedBasis_ , false );
    bases[1][2] = rcp( &openBasis_ , false );
    bases[2][0] = rcp( &openBasis_ , false );
    bases[2][1] = rcp( &openBasis_ , false );
    bases[2][2] = rcp( &closedBasis_ , false );
    this->setBases( bases );

    initializeTags();
    this->basisTagsAreSet_ = true;
  }
  

  
  template<class Scalar, class ArrayScalar>
  void Basis_HDIV_HEX_In_FEM<Scalar, ArrayScalar>::initializeTags() {
    
    // Basis-dependent intializations
    int tagSize  = 4;        // size of DoF tag
    int posScDim = 0;        // position in the tag, counting from 0, of the subcell dim 
    int posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
    int posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell
    
    std::vector<int> tags( tagSize * this->getCardinality() );
    
    const std::vector<std::vector<int> >& closedDofTags = closedBasis_.getAllDofTags();
    const std::vector<std::vector<int> >& openDofTags = openBasis_.getAllDofTags();

    std::map<int,std::map<int,int> > total_dof_per_entity;
    std::map<int,std::map<int,int> > current_dof_per_entity;

    // four vertices
    for (int i=0;i<4;i++) {
      total_dof_per_entity[0][i] = 0;
      current_dof_per_entity[0][i] = 0;
    }
    // twelve edges
    for (int i=0;i<12;i++) {
      total_dof_per_entity[1][i] = 0;
      current_dof_per_entity[1][i] = 0;
    }
    // six faces
    for (int i=0;i<6;i++) {
      total_dof_per_entity[2][i] = 0;
      current_dof_per_entity[2][i] = 0;
    }
    total_dof_per_entity[3][0] = 0;
    current_dof_per_entity[3][0] = 0;

    // tally dof on each facet.  none on vertex or edge
    for (int i=0;i<6;i++) {
      total_dof_per_entity[2][i] = openBasis_.getCardinality() * openBasis_.getCardinality();
    }

    total_dof_per_entity[2][0] = this->getCardinality() - 6 * openBasis_.getCardinality()* openBasis_.getCardinality();

    int tagcur = 0;
    // loop over the x-component basis functions, which are (psi(x)phi(y)phi(z),0,0)
    // for psi in the closed basis and phi in the open
    for (int k=0;k<openBasis_.getCardinality();k++) {
      const int odimk = openDofTags[k][0];
      const int oentk = openDofTags[k][1];
      for (int j=0;j<openBasis_.getCardinality();j++) {
	const int odimj = openDofTags[j][0];
	const int oentj = openDofTags[j][1];
	for (int i=0;i<closedBasis_.getCardinality();i++) {
	  const int cdim = closedDofTags[i][0];
	  const int cent = closedDofTags[i][1];
	  int dofdim;
	  int dofent;
	  ProductTopology::lineProduct3d(cdim,cent,odimj,oentj,odimk,oentk,dofdim,dofent);
	  tags[4*tagcur] = dofdim;
	  tags[4*tagcur+1] = dofent;
	  tags[4*tagcur+2] = current_dof_per_entity[dofdim][dofent];
	  current_dof_per_entity[dofdim][dofent]++;
	  tags[4*tagcur+3] = total_dof_per_entity[dofdim][dofent];
	  tagcur++;
	}
      }
    }

    // now we have to do it for the y-component basis functions, which are
    // (0,phi(x)psi(y)phi(z),0) for psi in the closed basis and phi in the open
    for (int k=0;k<openBasis_.getCardinality();k++) {
      const int odimk = openDofTags[k][0];
      const int oentk = openDofTags[k][1];
      for (int j=0;j<closedBasis_.getCardinality();j++) {
	const int cdim = closedDofTags[j][0];
	const int cent = closedDofTags[j][1];
	for (int i=0;i<openBasis_.getCardinality();i++) {
	  const int odimi = openDofTags[i][0];
	  const int oenti = openDofTags[i][1];
	  int dofdim;
	  int dofent;
	  ProductTopology::lineProduct3d(odimi,oenti,cdim,cent,odimk,oentk,dofdim,dofent);
	  tags[4*tagcur] = dofdim;
	  tags[4*tagcur+1] = dofent;
	  tags[4*tagcur+2] = current_dof_per_entity[dofdim][dofent];
	  current_dof_per_entity[dofdim][dofent]++;
	  tags[4*tagcur+3] = total_dof_per_entity[dofdim][dofent];
	  tagcur++;
	}
      }
    }

    // now we have to do it for the z-component basis functions, which are
    // (0,0,phi(x)phi(y)psi(z)) for psi in the closed basis and phi in the open
    for (int k=0;k<closedBasis_.getCardinality();k++) {
      const int cdim = closedDofTags[k][0];
      const int cent = closedDofTags[k][1];
      for (int j=0;j<openBasis_.getCardinality();j++) {
	const int odimj = openDofTags[j][0];
	const int oentj = openDofTags[j][1];
	for (int i=0;i<openBasis_.getCardinality();i++) {
	  const int odimi = openDofTags[i][0];
	  const int oenti = openDofTags[i][1];
	  int dofdim;
	  int dofent;
	  ProductTopology::lineProduct3d(odimi,oenti,odimj,oentj,cdim,cent,dofdim,dofent);
	  tags[4*tagcur] = dofdim;
	  tags[4*tagcur+1] = dofent;
	  tags[4*tagcur+2] = current_dof_per_entity[dofdim][dofent];
	  current_dof_per_entity[dofdim][dofent]++;
	  tags[4*tagcur+3] = total_dof_per_entity[dofdim][dofent];
	  tagcur++;
	}
      }
    }
  
//     for (int i=0;i<this->getCardinality();i++) {
//       for (int j=0;j<4;j++) {
// 	std::cout << tags[4*i+j] << " ";
//       }
//       std::cout << std::endl;
//     }
  
    // Basis-independent function sets tag and enum data in tagToOrdinal_ and ordinalToTag_ arrays:
    Intrepid2::setOrdinalTagData(this -> tagToOrdinal_,
				this -> ordinalToTag_,
				&(tags[0]),
				this -> basisCardinality_,
				tagSize,
				posScDim,
				posScOrd,
				posDfOrd);
  }


  template<class Scalar, class ArrayScalar>
  void Basis_HDIV_HEX_In_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar &        outputValues,
							      const ArrayScalar &  inputPoints,
							      const EOperator      operatorType) const {
    
    // Verify arguments
#ifdef HAVE_INTREPID2_DEBUG
    Intrepid2::getValues_HDIV_Args<Scalar, ArrayScalar>(outputValues,
						       inputPoints,
						       operatorType,
						       this -> getBaseCellTopology(),
						       this -> getCardinality() );
#endif
    
    // Number of evaluation points = dim 0 of inputPoints
    int dim0 = inputPoints.dimension(0);
    
    // separate out points
    ArrayScalar xPoints(dim0,1);
    ArrayScalar yPoints(dim0,1);
    ArrayScalar zPoints(dim0,1);

    
    for (int i=0;i<dim0;i++) {
      xPoints(i,0) = inputPoints(i,0);
      yPoints(i,0) = inputPoints(i,1);
      zPoints(i,0) = inputPoints(i,2);
    }
    
    switch (operatorType) {
    case OPERATOR_VALUE:
      {
	ArrayScalar closedBasisValsXPts( closedBasis_.getCardinality() , dim0 );
	ArrayScalar closedBasisValsYPts( closedBasis_.getCardinality() , dim0 );
	ArrayScalar closedBasisValsZPts( closedBasis_.getCardinality() , dim0 );
	ArrayScalar openBasisValsXPts( openBasis_.getCardinality() , dim0 );
	ArrayScalar openBasisValsYPts( openBasis_.getCardinality() , dim0 );
	ArrayScalar openBasisValsZPts( openBasis_.getCardinality() , dim0 );
	
	closedBasis_.getValues( closedBasisValsXPts , xPoints , OPERATOR_VALUE );
	closedBasis_.getValues( closedBasisValsYPts , yPoints , OPERATOR_VALUE );
	closedBasis_.getValues( closedBasisValsZPts , zPoints , OPERATOR_VALUE );
	openBasis_.getValues( openBasisValsXPts , xPoints , OPERATOR_VALUE );
	openBasis_.getValues( openBasisValsYPts , yPoints , OPERATOR_VALUE );
	openBasis_.getValues( openBasisValsZPts , zPoints , OPERATOR_VALUE );
	
	int bfcur = 0;
	// x component bfs are (closed(x) open(y) open(z),0,0)
	for (int k=0;k<openBasis_.getCardinality();k++) {
	  for (int j=0;j<openBasis_.getCardinality();j++) {
	    for (int i=0;i<closedBasis_.getCardinality();i++) {
	      for (int l=0;l<dim0;l++) {
		outputValues(bfcur,l,0) = closedBasisValsXPts(i,l) * openBasisValsYPts(j,l) * openBasisValsZPts(k,l);
		outputValues(bfcur,l,1) = 0.0;
		outputValues(bfcur,l,2) = 0.0;
	      }
	      bfcur++;
	    }
	  }
	}
	
	// y component bfs are (0,open(x) closed(y) open(z),0)
	for (int k=0;k<openBasis_.getCardinality();k++) {
	  for (int j=0;j<closedBasis_.getCardinality();j++) {
	    for (int i=0;i<openBasis_.getCardinality();i++) {
	      for (int l=0;l<dim0;l++) {
		outputValues(bfcur,l,0) = 0.0;
		outputValues(bfcur,l,1) = openBasisValsXPts(i,l) * closedBasisValsYPts(j,l) * openBasisValsZPts(k,l);
		outputValues(bfcur,l,2) = 0.0;
	      }
	      bfcur++;
	    }
	  }
	}

	// z component bfs are (0,0,open(x) open(y) closed(z))
	for (int k=0;k<closedBasis_.getCardinality();k++) {
	  for (int j=0;j<openBasis_.getCardinality();j++) {
	    for (int i=0;i<openBasis_.getCardinality();i++) {
	      for (int l=0;l<dim0;l++) {
		outputValues(bfcur,l,0) = 0.0;
		outputValues(bfcur,l,1) = 0.0;
		outputValues(bfcur,l,2) = openBasisValsXPts(i,l) * openBasisValsYPts(j,l) * closedBasisValsZPts(k,l);
	      }
	      bfcur++;
	    }
	  }
	}


      }
      break;
    case OPERATOR_DIV:
      {
	ArrayScalar closedBasisDerivsXPts( closedBasis_.getCardinality() , dim0 , 1 );
	ArrayScalar closedBasisDerivsYPts( closedBasis_.getCardinality() , dim0 , 1 );
	ArrayScalar closedBasisDerivsZPts( closedBasis_.getCardinality() , dim0 , 1 );
	ArrayScalar openBasisValsXPts( openBasis_.getCardinality() , dim0 );
	ArrayScalar openBasisValsYPts( openBasis_.getCardinality() , dim0 );
	ArrayScalar openBasisValsZPts( openBasis_.getCardinality() , dim0 );
	
	closedBasis_.getValues( closedBasisDerivsXPts , xPoints , OPERATOR_D1 );
	closedBasis_.getValues( closedBasisDerivsYPts , yPoints , OPERATOR_D1 );
	closedBasis_.getValues( closedBasisDerivsZPts , zPoints , OPERATOR_D1 );
	openBasis_.getValues( openBasisValsXPts , xPoints , OPERATOR_VALUE );
	openBasis_.getValues( openBasisValsYPts , yPoints , OPERATOR_VALUE );
	openBasis_.getValues( openBasisValsZPts , zPoints , OPERATOR_VALUE );
	
	int bfcur = 0;
	
	// x component basis functions first
	for (int k=0;k<openBasis_.getCardinality();k++) {
	  for (int j=0;j<openBasis_.getCardinality();j++) {
	    for (int i=0;i<closedBasis_.getCardinality();i++) {
	      for (int l=0;l<dim0;l++) {
		outputValues(bfcur,l) = closedBasisDerivsXPts(i,l,0) * openBasisValsYPts(j,l) * openBasisValsZPts(k,l);
	      }
	      bfcur++;
	    }
	  }
	}
	
	// now y component basis functions
	for (int k=0;k<openBasis_.getCardinality();k++) {
	  for (int j=0;j<closedBasis_.getCardinality();j++) {
	    for (int i=0;i<openBasis_.getCardinality();i++) {
	      for (int l=0;l<dim0;l++) {
		outputValues(bfcur,l) = openBasisValsXPts(i,l) * closedBasisDerivsYPts(j,l,0) * openBasisValsZPts(k,l);
	      }
	      bfcur++;
	    }
	  }
	}

	// now z component basis functions
	for (int k=0;k<closedBasis_.getCardinality();k++) {
	  for (int j=0;j<openBasis_.getCardinality();j++) {
	    for (int i=0;i<openBasis_.getCardinality();i++) {
	      for (int l=0;l<dim0;l++) {
		outputValues(bfcur,l) = openBasisValsXPts(i,l) * openBasisValsYPts(j,l) * closedBasisDerivsZPts(k,l,0);
	      }
	      bfcur++;
	    }
	  }
	}
      }
      break;
    case OPERATOR_CURL:
      TEUCHOS_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_CURL), std::invalid_argument,
			  ">>> ERROR (Basis_HDIV_HEX_In_FEM): CURL is invalid operator for HDIV Basis Functions");
      break;
      
    case OPERATOR_GRAD:
      TEUCHOS_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_GRAD), std::invalid_argument,
			  ">>> ERROR (Basis_HDIV_HEX_In_FEM): GRAD is invalid operator for HDIV Basis Functions");
      break;
      
    case OPERATOR_D1:
    case OPERATOR_D2:
    case OPERATOR_D3:
    case OPERATOR_D4:
    case OPERATOR_D5:
    case OPERATOR_D6:
    case OPERATOR_D7:
    case OPERATOR_D8:
    case OPERATOR_D9:
    case OPERATOR_D10:
      TEUCHOS_TEST_FOR_EXCEPTION( ( (operatorType == OPERATOR_D1)    ||
			    (operatorType == OPERATOR_D2)    ||
			    (operatorType == OPERATOR_D3)    ||
			    (operatorType == OPERATOR_D4)    ||
			    (operatorType == OPERATOR_D5)    ||
			    (operatorType == OPERATOR_D6)    ||
			    (operatorType == OPERATOR_D7)    ||
			    (operatorType == OPERATOR_D8)    ||
			    (operatorType == OPERATOR_D9)    ||
			    (operatorType == OPERATOR_D10) ),
			  std::invalid_argument,
			  ">>> ERROR (Basis_HDIV_HEX_In_FEM): Invalid operator type");
      break;
      
    default:
      TEUCHOS_TEST_FOR_EXCEPTION( ( (operatorType != OPERATOR_VALUE) &&
			    (operatorType != OPERATOR_GRAD)  &&
			    (operatorType != OPERATOR_CURL)  &&
			    (operatorType != OPERATOR_DIV)   &&
			    (operatorType != OPERATOR_D1)    &&
			    (operatorType != OPERATOR_D2)    &&
			    (operatorType != OPERATOR_D3)    &&
			    (operatorType != OPERATOR_D4)    &&
			    (operatorType != OPERATOR_D5)    &&
			    (operatorType != OPERATOR_D6)    &&
			    (operatorType != OPERATOR_D7)    &&
			    (operatorType != OPERATOR_D8)    &&
			    (operatorType != OPERATOR_D9)    &&
			    (operatorType != OPERATOR_D10) ),
			  std::invalid_argument,
			  ">>> ERROR (Basis_HDIV_HEX_In_FEM): Invalid operator type");
    }
  }
  
  
  
  template<class Scalar, class ArrayScalar>
  void Basis_HDIV_HEX_In_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar&           outputValues,
							      const ArrayScalar &    inputPoints,
							      const ArrayScalar &    cellVertices,
							      const EOperator        operatorType) const {
    TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
			">>> ERROR (Basis_HDIV_HEX_In_FEM): FEM Basis calling an FVD member function");
  }

  template<class Scalar, class ArrayScalar>
  void Basis_HDIV_HEX_In_FEM<Scalar,ArrayScalar>::getDofCoords(ArrayScalar & DofCoords) const
  {
    // x-component basis functions
    int cur = 0;

    for (int k=0;k<openPts_.dimension(0);k++)
      {
	for (int j=0;j<openPts_.dimension(0);j++)
	  {
	    for (int i=0;i<closedPts_.dimension(0);i++)
	      {
		DofCoords(cur,0) = closedPts_(i,0);
		DofCoords(cur,1) = openPts_(j,0);
		DofCoords(cur,2) = openPts_(k,0);
		cur++;
	      }
	  }
      }

    // y-component basis functions
    for (int k=0;k<openPts_.dimension(0);k++)
      {
	for (int j=0;j<closedPts_.dimension(0);j++)
	  {
	    for (int i=0;i<openPts_.dimension(0);i++)
	      {
		DofCoords(cur,0) = openPts_(i,0);
		DofCoords(cur,1) = closedPts_(j,0);
		DofCoords(cur,2) = openPts_(k,0);
		cur++;
	      }
	  }
      }

    // z-component basis functions
    for (int k=0;k<closedPts_.dimension(0);k++)
      {
	for (int j=0;j<openPts_.dimension(0);j++)
	  {
	    for (int i=0;i<openPts_.dimension(0);i++)
	      {
		DofCoords(cur,0) = openPts_(i,0);
		DofCoords(cur,1) = openPts_(j,0);
		DofCoords(cur,2) = closedPts_(k,0);
		cur++;
	      }
	  }
      }



    return;
  }
  
}// namespace Intrepid2
