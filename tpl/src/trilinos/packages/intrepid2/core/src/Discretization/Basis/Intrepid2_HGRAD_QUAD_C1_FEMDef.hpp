#ifndef INTREPID2_HGRAD_QUAD_C1_FEMDEF_HPP
#define INTREPID2_HGRAD_QUAD_C1_FEMDEF_HPP
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

/** \file   Intrepid_HGRAD_QUAD_C1_FEMDef.hpp
    \brief  Definition file for bi-linear FEM basis functions for H(grad) functions on QUAD cells.
    \author Created by P. Bochev and D. Ridzal.
*/

namespace Intrepid2 {

  template<class Scalar, class ArrayScalar>
  Basis_HGRAD_QUAD_C1_FEM<Scalar, ArrayScalar>::Basis_HGRAD_QUAD_C1_FEM()
  {
    this -> basisCardinality_  = 4;
    this -> basisDegree_       = 1;    
    this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >() );
    this -> basisType_         = BASIS_FEM_DEFAULT;
    this -> basisCoordinates_  = COORDINATES_CARTESIAN;
    this -> basisTagsAreSet_   = false;

    initializeTags();
    this->basisTagsAreSet_ = true;
  }
    
  
template<class Scalar, class ArrayScalar>
void Basis_HGRAD_QUAD_C1_FEM<Scalar, ArrayScalar>::initializeTags() {
  
  // Basis-dependent intializations
  int tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
  int posScDim = 0;        // position in the tag, counting from 0, of the subcell dim 
  int posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
  int posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell

  // An array with local DoF tags assigned to basis functions, in the order of their local enumeration 
  int tags[]  = { 0, 0, 0, 1,
                  0, 1, 0, 1,
                  0, 2, 0, 1,
                  0, 3, 0, 1};
  
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
void Basis_HGRAD_QUAD_C1_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar &        outputValues,
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
  
  // Temporaries: (x,y) coordinates of the evaluation point
  Scalar x = 0.0;                                    
  Scalar y = 0.0;                                    
  
  switch (operatorType) {
    
    case OPERATOR_VALUE:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0, 0);
        y = inputPoints(i0, 1);
        
        // outputValues is a rank-2 array with dimensions (basisCardinality_, dim0)
        outputValues(0, i0) = (1.0 - x)*(1.0 - y)/4.0;
        outputValues(1, i0) = (1.0 + x)*(1.0 - y)/4.0;
        outputValues(2, i0) = (1.0 + x)*(1.0 + y)/4.0;
        outputValues(3, i0) = (1.0 - x)*(1.0 + y)/4.0;
      }
      break;
      
    case OPERATOR_GRAD:
    case OPERATOR_D1:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0,0);
        y = inputPoints(i0,1);
        
        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, spaceDim)
        outputValues(0, i0, 0) = -(1.0 - y)/4.0;
        outputValues(0, i0, 1) = -(1.0 - x)/4.0;
        
        outputValues(1, i0, 0) =  (1.0 - y)/4.0;
        outputValues(1, i0, 1) = -(1.0 + x)/4.0;
        
        outputValues(2, i0, 0) =  (1.0 + y)/4.0;
        outputValues(2, i0, 1) =  (1.0 + x)/4.0;
 
        outputValues(3, i0, 0) = -(1.0 + y)/4.0;
        outputValues(3, i0, 1) =  (1.0 - x)/4.0;
      }
      break;
      
    case OPERATOR_CURL:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0,0);
        y = inputPoints(i0,1);
        
        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, spaceDim)
        outputValues(0, i0, 0) = -(1.0 - x)/4.0;
        outputValues(0, i0, 1) =  (1.0 - y)/4.0;
        
        outputValues(1, i0, 0) = -(1.0 + x)/4.0;
        outputValues(1, i0, 1) = -(1.0 - y)/4.0;
        
        outputValues(2, i0, 0) =  (1.0 + x)/4.0;
        outputValues(2, i0, 1) = -(1.0 + y)/4.0;
        
        outputValues(3, i0, 0) =  (1.0 - x)/4.0;
        outputValues(3, i0, 1) =  (1.0 + y)/4.0;
      }
      break;
      
    case OPERATOR_DIV:
      TEUCHOS_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_DIV), std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_QUAD_C1_FEM): DIV is invalid operator for rank-0 (scalar) functions in 2D");
      break;
      
    case OPERATOR_D2:
      for (int i0 = 0; i0 < dim0; i0++) {
        
        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, D2Cardinality=3) 
        outputValues(0, i0, 0) =  0.0;
        outputValues(0, i0, 1) =  0.25;
        outputValues(0, i0, 2) =  0.0;

        outputValues(1, i0, 0) =  0.0;
        outputValues(1, i0, 1) = -0.25;
        outputValues(1, i0, 2) =  0.0;
        
        outputValues(2, i0, 0) =  0.0;
        outputValues(2, i0, 1) =  0.25;
        outputValues(2, i0, 2) =  0.0;
        
        outputValues(3, i0, 0) =  0.0;
        outputValues(3, i0, 1) = -0.25;
        outputValues(3, i0, 2) =  0.0;
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
                          ">>> ERROR (Basis_HGRAD_QUAD_C1_FEM): Invalid operator type");
  }
}


  
template<class Scalar, class ArrayScalar>
void Basis_HGRAD_QUAD_C1_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar&           outputValues,
                                                             const ArrayScalar &    inputPoints,
                                                             const ArrayScalar &    cellVertices,
                                                             const EOperator        operatorType) const {
  TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                      ">>> ERROR (Basis_HGRAD_QUAD_C1_FEM): FEM Basis calling an FVD member function");
}



template<class Scalar, class ArrayScalar>
void Basis_HGRAD_QUAD_C1_FEM<Scalar, ArrayScalar>::getDofCoords(ArrayScalar & DofCoords) const {
#ifdef HAVE_INTREPID2_DEBUG
  // Verify rank of output array.
  TEUCHOS_TEST_FOR_EXCEPTION( !(DofCoords.rank() == 2), std::invalid_argument,
                      ">>> ERROR: (Intrepid2::Basis_HGRAD_QUAD_C1_FEM::getDofCoords) rank = 2 required for DofCoords array");
  // Verify 0th dimension of output array.
  TEUCHOS_TEST_FOR_EXCEPTION( !( static_cast<index_type>(DofCoords.dimension(0)) == static_cast<index_type>(this -> basisCardinality_) ), std::invalid_argument,
                      ">>> ERROR: (Intrepid2::Basis_HGRAD_QUAD_C1_FEM::getDofCoords) mismatch in number of DoF and 0th dimension of DofCoords array");
  // Verify 1st dimension of output array.
  TEUCHOS_TEST_FOR_EXCEPTION( !( static_cast<index_type>(DofCoords.dimension(1)) == static_cast<index_type>(this -> basisCellTopology_.getDimension()) ), std::invalid_argument,
                      ">>> ERROR: (Intrepid2::Basis_HGRAD_QUAD_C1_FEM::getDofCoords) incorrect reference cell (1st) dimension in DofCoords array");
#endif

  DofCoords(0,0) = -1.0;   DofCoords(0,1) = -1.0;
  DofCoords(1,0) =  1.0;   DofCoords(1,1) = -1.0;
  DofCoords(2,0) =  1.0;   DofCoords(2,1) =  1.0;
  DofCoords(3,0) = -1.0;   DofCoords(3,1) =  1.0;
}

}// namespace Intrepid2
#endif
