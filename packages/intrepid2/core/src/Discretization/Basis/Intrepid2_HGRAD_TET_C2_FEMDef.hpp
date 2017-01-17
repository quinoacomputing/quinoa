#ifndef INTREPID2_HGRAD_TET_C2_FEMDEF_HPP
#define INTREPID2_HGRAD_TET_C2_FEMDEF_HPP
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

/** \file   Intrepid_HGRAD_TET_C2_FEMDef.hpp
    \brief  Definition file for the default H(grad)-compatible FEM basis of degree 2 on Tetrahedron cell
    \author Created by P. Bochev and D. Ridzal.
*/

namespace Intrepid2 {

  template<class Scalar, class ArrayScalar>
  Basis_HGRAD_TET_C2_FEM<Scalar, ArrayScalar>::Basis_HGRAD_TET_C2_FEM()
  {
    this -> basisCardinality_  = 10;
    this -> basisDegree_       = 2;    
    this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<4> >() );
    this -> basisType_         = BASIS_FEM_DEFAULT;
    this -> basisCoordinates_  = COORDINATES_CARTESIAN;
    this -> basisTagsAreSet_   = false;

    initializeTags();
    this->basisTagsAreSet_ = true;
  }
  
  
template<class Scalar, class ArrayScalar>
void Basis_HGRAD_TET_C2_FEM<Scalar, ArrayScalar>::initializeTags() {
  
  // Basis-dependent intializations
  int tagSize  = 4;        // size of DoF tag
  int posScDim = 0;        // position in the tag, counting from 0, of the subcell dim 
  int posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
  int posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell

  // An array with local DoF tags assigned to basis functions, in the order of their local enumeration 
  int tags[]  = { 0, 0, 0, 1,
                  0, 1, 0, 1,
                  0, 2, 0, 1,
                  0, 3, 0, 1,
                  1, 0, 0, 1,
                  1, 1, 0, 1,
                  1, 2, 0, 1,
                  1, 3, 0, 1,
                  1, 4, 0, 1,
                  1, 5, 0, 1,
  };
  
  // Basis-independent function sets tag and enum data in tagToOrdinal_ and ordinalToTag_ arrays:
  Intrepid2::setOrdinalTagData(this -> tagToOrdinal_,
                              this -> ordinalToTag_,
                              tags,
                              this -> basisCardinality_,
                              tagSize,
                              posScDim,
                              posScOrd,
                              posDfOrd);
}



template<class Scalar, class ArrayScalar>
void Basis_HGRAD_TET_C2_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar &        outputValues,
                                                             const ArrayScalar &  inputPoints,
                                                             const EOperator      operatorType) const {
  
  // Verify arguments
#ifdef HAVE_INTREPID2_DEBUG
  Intrepid2::getValues_HGRAD_Args<Scalar, ArrayScalar>(outputValues,
                                                      inputPoints,
                                                      operatorType,
                                                      this -> getBaseCellTopology(),
                                                      this -> getCardinality() );
#endif
  
  // Number of evaluation points = dim 0 of inputPoints
  int dim0 = inputPoints.dimension(0);  
  
  // Temporaries: (x,y,z) coordinates of the evaluation point
  Scalar x = 0.0;                                    
  Scalar y = 0.0;   
  Scalar z = 0.0;
  
  switch (operatorType) {
    
    case OPERATOR_VALUE:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0, 0);
        y = inputPoints(i0, 1);
        z = inputPoints(i0, 2);
        
        // outputValues is a rank-2 array with dimensions (basisCardinality_, dim0)
        outputValues(0, i0) = (-1. + x + y + z)*(-1. + 2.*x + 2.*y + 2.*z);
        outputValues(1, i0) = x*(-1. + 2.*x);
        outputValues(2, i0) = y*(-1. + 2.*y);
        outputValues(3, i0) = z*(-1. + 2.*z);

        outputValues(4, i0) = -4.*x*(-1. + x + y + z);
        outputValues(5, i0) =  4.*x*y;
        outputValues(6, i0) = -4.*y*(-1. + x + y + z);
        outputValues(7, i0) = -4.*z*(-1. + x + y + z);
        outputValues(8, i0) =  4.*x*z;
        outputValues(9, i0) =  4.*y*z;

      }
      break;
      
    case OPERATOR_GRAD:
    case OPERATOR_D1:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0, 0);
        y = inputPoints(i0, 1);
        z = inputPoints(i0, 2);
        
        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, spaceDim)
        outputValues(0, i0, 0) = -3.+ 4.*x + 4.*y + 4.*z;
        outputValues(0, i0, 1) = -3.+ 4.*x + 4.*y + 4.*z;
        outputValues(0, i0, 2) = -3.+ 4.*x + 4.*y + 4.*z; 
        
        outputValues(1, i0, 0) = -1.+ 4.*x; 
        outputValues(1, i0, 1) =  0.;
        outputValues(1, i0, 2) =  0.;
        
        outputValues(2, i0, 0) =  0.;        
        outputValues(2, i0, 1) = -1.+ 4.*y;
        outputValues(2, i0, 2) =  0.;
  
        outputValues(3, i0, 0) =  0.;         
        outputValues(3, i0, 1) =  0.;
        outputValues(3, i0, 2) = -1.+ 4.*z;
 
        
        outputValues(4, i0, 0) = -4.*(-1.+ 2*x + y + z);         
        outputValues(4, i0, 1) = -4.*x;
        outputValues(4, i0, 2) = -4.*x;
        
        outputValues(5, i0, 0) =  4.*y;       
        outputValues(5, i0, 1) =  4.*x; 
        outputValues(5, i0, 2) =  0.;

        outputValues(6, i0, 0) = -4.*y;          
        outputValues(6, i0, 1) = -4.*(-1.+ x + 2*y + z);
        outputValues(6, i0, 2) = -4.*y; 

        outputValues(7, i0, 0) = -4.*z;          
        outputValues(7, i0, 1) = -4.*z;
        outputValues(7, i0, 2) = -4.*(-1.+ x + y + 2*z);

        outputValues(8, i0, 0) =  4.*z;     
        outputValues(8, i0, 1) =  0.;
        outputValues(8, i0, 2) =  4.*x;

        outputValues(9, i0, 0) =  0.;         
        outputValues(9, i0, 1) =  4.*z;
        outputValues(9, i0, 2) =  4.*y;
        
      }
      break;
      
    case OPERATOR_CURL:
      TEUCHOS_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_CURL), std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_TET_C2_FEM): CURL is invalid operator for rank-0 (scalar) functions in 3D");
      break;
      
    case OPERATOR_DIV:
      TEUCHOS_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_DIV), std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_TET_C2_FEM): DIV is invalid operator for rank-0 (scalar) functions in 3D");
      break;
      
    case OPERATOR_D2:
      for (int i0 = 0; i0 < dim0; i0++) {
          
        outputValues(0, i0, 0) =  4.;
        outputValues(0, i0, 1) =  4.;
        outputValues(0, i0, 2) =  4.;
        outputValues(0, i0, 3) =  4.;
        outputValues(0, i0, 4) =  4.;
        outputValues(0, i0, 5) =  4.;
        
        outputValues(1, i0, 0) =  4.;
        outputValues(1, i0, 1) =  0.;
        outputValues(1, i0, 2) =  0.;
        outputValues(1, i0, 3) =  0.;
        outputValues(1, i0, 4) =  0.;
        outputValues(1, i0, 5) =  0.;

        outputValues(2, i0, 0) =  0.;
        outputValues(2, i0, 1) =  0.;
        outputValues(2, i0, 2) =  0.;
        outputValues(2, i0, 3) =  4.;
        outputValues(2, i0, 4) =  0.;
        outputValues(2, i0, 5) =  0.;

        outputValues(3, i0, 0) =  0.;
        outputValues(3, i0, 1) =  0.;
        outputValues(3, i0, 2) =  0.;
        outputValues(3, i0, 3) =  0.;
        outputValues(3, i0, 4) =  0.;
        outputValues(3, i0, 5) =  4.;

        outputValues(4, i0, 0) = -8.;
        outputValues(4, i0, 1) = -4.;
        outputValues(4, i0, 2) = -4.;
        outputValues(4, i0, 3) =  0.;
        outputValues(4, i0, 4) =  0.;
        outputValues(4, i0, 5) =  0.;

        outputValues(5, i0, 0) =  0.;
        outputValues(5, i0, 1) =  4.;
        outputValues(5, i0, 2) =  0.;
        outputValues(5, i0, 3) =  0.;
        outputValues(5, i0, 4) =  0.;
        outputValues(5, i0, 5) =  0.;

        outputValues(6, i0, 0) =  0.;
        outputValues(6, i0, 1) = -4.;
        outputValues(6, i0, 2) =  0.;
        outputValues(6, i0, 3) = -8.;
        outputValues(6, i0, 4) = -4.;
        outputValues(6, i0, 5) =  0;

        outputValues(7, i0, 0) =  0.;
        outputValues(7, i0, 1) =  0.;
        outputValues(7, i0, 2) = -4.;
        outputValues(7, i0, 3) =  0.;
        outputValues(7, i0, 4) = -4.;
        outputValues(7, i0, 5) = -8.;

        outputValues(8, i0, 0) =  0.;
        outputValues(8, i0, 1) =  0.;
        outputValues(8, i0, 2) =  4.;
        outputValues(8, i0, 3) =  0.;
        outputValues(8, i0, 4) =  0.;
        outputValues(8, i0, 5) =  0.;

        outputValues(9, i0, 0) =  0.;
        outputValues(9, i0, 1) =  0.;
        outputValues(9, i0, 2) =  0.;
        outputValues(9, i0, 3) =  0.;
        outputValues(9, i0, 4) =  4.;
        outputValues(9, i0, 5) =  0.;
      }
      break;
      
    case OPERATOR_D3:
    case OPERATOR_D4:
    case OPERATOR_D5:
    case OPERATOR_D6:
    case OPERATOR_D7:
    case OPERATOR_D8:
    case OPERATOR_D9:
    case OPERATOR_D10:
      {
        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, DkCardinality)
        int DkCardinality = Intrepid2::getDkCardinality(operatorType, 
                                                       this -> basisCellTopology_.getDimension() );
        for(int dofOrd = 0; dofOrd < this -> basisCardinality_; dofOrd++) {
          for (int i0 = 0; i0 < dim0; i0++) {
            for(int dkOrd = 0; dkOrd < DkCardinality; dkOrd++){
              outputValues(dofOrd, i0, dkOrd) = 0.0;
            }
          }
        }
      }
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION( !( Intrepid2::isValidOperator(operatorType) ), std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_TET_C2_FEM): Invalid operator type");
  }
}


  
template<class Scalar, class ArrayScalar>
void Basis_HGRAD_TET_C2_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar&           outputValues,
                                                             const ArrayScalar &    inputPoints,
                                                             const ArrayScalar &    cellVertices,
                                                             const EOperator        operatorType) const {
  TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                      ">>> ERROR (Basis_HGRAD_TET_C2_FEM): FEM Basis calling an FVD member function");
}
}// namespace Intrepid2
#endif
