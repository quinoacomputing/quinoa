#ifndef INTREPID_HGRAD_PYR_I2_FEMDEF_HPP
#define INTREPID_HGRAD_PYR_I2_FEMDEF_HPP
// @HEADER
// ************************************************************************
//
//                           Intrepid Package
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
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_HGRAD_PYR_I2_FEMDef.hpp
    \brief  Definition file for an incomplete bi-quadratic FEM basis function for H(grad) functions on PYRAMID cells. 
    \author Created by P. Bochev and D. Ridzal and K. Peterson.
*/

namespace Intrepid {

  template<class Scalar, class ArrayScalar>
  Basis_HGRAD_PYR_I2_FEM<Scalar, ArrayScalar>::Basis_HGRAD_PYR_I2_FEM()
  {
    this -> basisCardinality_  = 13;
    this -> basisDegree_       = 2;    
    this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Pyramid<5> >() );
    this -> basisType_         = BASIS_FEM_DEFAULT;
    this -> basisCoordinates_  = COORDINATES_CARTESIAN;
    this -> basisTagsAreSet_   = false;
  }
  
  
template<class Scalar, class ArrayScalar>
void Basis_HGRAD_PYR_I2_FEM<Scalar, ArrayScalar>::initializeTags() {
  
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
                  0, 4, 0, 1,
                  1, 0, 0, 1,
                  1, 1, 0, 1,
                  1, 2, 0, 1,
                  1, 3, 0, 1,
                  1, 4, 0, 1,
                  1, 5, 0, 1,
                  1, 6, 0, 1,
                  1, 7, 0, 1,
  };
  
  // Basis-independent function sets tag and enum data in tagToOrdinal_ and ordinalToTag_ arrays:
  Intrepid::setOrdinalTagData(this -> tagToOrdinal_,
                              this -> ordinalToTag_,
                              tags,
                              this -> basisCardinality_,
                              tagSize,
                              posScDim,
                              posScOrd,
                              posDfOrd);
}



template<class Scalar, class ArrayScalar>
void Basis_HGRAD_PYR_I2_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar &    outputValues,
                                                             const ArrayScalar &  inputPoints,
                                                             const EOperator      operatorType) const {
  
  // Verify arguments
#ifdef HAVE_INTREPID_DEBUG
  Intrepid::getValues_HGRAD_Args<Scalar, ArrayScalar>(outputValues,
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
  const Scalar eps = std::numeric_limits<Scalar>::epsilon( );
  
  switch (operatorType) {
    
    case OPERATOR_VALUE:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0, 0);
        y = inputPoints(i0, 1);
        z = inputPoints(i0, 2);

        //be sure that the basis functions are defined when z is very close to 1.

        if(fabs(z-1.0) < eps) {
          if(z <= 1.0) z = 1.0-eps;
          else  z = 1.0+eps;
          //z = 1.0-eps;
        }

        Scalar w = 1.0/(1.0 - z);
        
        // outputValues is a rank-2 array with dimensions (basisCardinality_, dim0)
        outputValues(0, i0) = 0.25 * (-x - y - 1.0)*((1.0-x)*(1.0-y) - z + x*y*z*w);
        outputValues(1, i0) = 0.25 * ( x - y - 1.0)*((1.0+x)*(1.0-y) - z - x*y*z*w);
        outputValues(2, i0) = 0.25 * ( x + y - 1.0)*((1.0+x)*(1.0+y) - z + x*y*z*w);
        outputValues(3, i0) = 0.25 * (-x + y - 1.0)*((1.0-x)*(1.0+y) - z - x*y*z*w);

        outputValues(4, i0) =  z * (2.0*z - 1.0);

        outputValues(5, i0) = 0.5 * (1.0 + x - z)*(1.0 - x - z)*(1.0 - y - z)*w;
        outputValues(6, i0) = 0.5 * (1.0 + y - z)*(1.0 - y - z)*(1.0 + x - z)*w;
        outputValues(7, i0) = 0.5 * (1.0 + x - z)*(1.0 - x - z)*(1.0 + y - z)*w;
        outputValues(8, i0) = 0.5 * (1.0 + y - z)*(1.0 - y - z)*(1.0 - x - z)*w;

        outputValues(9, i0) = z*(1.0 - x - z)*(1.0 - y - z)*w;
        outputValues(10,i0) = z*(1.0 + x - z)*(1.0 - y - z)*w;
        outputValues(11,i0) = z*(1.0 + x - z)*(1.0 + y - z)*w;
        outputValues(12,i0) = z*(1.0 - x - z)*(1.0 + y - z)*w;
      }
      break;
      
    case OPERATOR_GRAD:
    case OPERATOR_D1:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0,0);
        y = inputPoints(i0,1);
        z = inputPoints(i0,2);

        //be sure that the basis functions are defined when z is very close to 1.

        if(fabs(z-1.0) < eps) {
          if(z <= 1.0) z = 1.0-eps;
          else  z = 1.0+eps;
        }

        Scalar w = 1.0/(1.0 - z);
        
        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, spaceDim)
        outputValues(0, i0, 0) =  0.25*(-1.0-x-y)*(-1.0+y + y*z*w) - 0.25*((1.0-x)*(1.0-y)-z + x*y*z*w);
        outputValues(0, i0, 1) =  0.25*(-1.0-x-y)*(-1.0+x + x*z*w) - 0.25*((1.0-x)*(1.0-y)-z + x*y*z*w);
        outputValues(0, i0, 2) =  0.25*(-1.0-x-y)*(-1.0 + x*y*w + x*y*z*w*w);
        
        outputValues(1, i0, 0) =  0.25*(-1.0+x-y)*( 1.0-y - y*z*w) + 0.25*((1.0+x)*(1.0-y)-z - x*y*z*w);
        outputValues(1, i0, 1) =  0.25*(-1.0+x-y)*(-1.0-x - x*z*w) - 0.25*((1.0+x)*(1.0-y)-z - x*y*z*w);
        outputValues(1, i0, 2) =  0.25*(-1.0+x-y)*(-1.0 - x*y*w - x*y*z*w*w);
        
        outputValues(2, i0, 0) =  0.25*(-1.0+x+y)*(1.0+y + y*z*w) + 0.25*((1.0+x)*(1.0+y)-z + x*y*z*w);
        outputValues(2, i0, 1) =  0.25*(-1.0+x+y)*(1.0+x + x*z*w) + 0.25*((1.0+x)*(1.0+y)-z + x*y*z*w);
        outputValues(2, i0, 2) =  0.25*(-1.0+x+y)*(-1.0 + x*y*w + x*y*z*w*w);
        
        outputValues(3, i0, 0) =  0.25*(-1.0-x+y)*(-1.0-y - y*z*w) - 0.25*((1.0-x)*(1.0+y)-z - x*y*z*w);
        outputValues(3, i0, 1) =  0.25*(-1.0-x+y)*( 1.0-x - x*z*w) + 0.25*((1.0-x)*(1.0+y)-z - x*y*z*w);
        outputValues(3, i0, 2) =  0.25*(-1.0-x+y)*(-1.0 - x*y*w - x*y*z*w*w);
        
        outputValues(4, i0, 0) =  0.0;
        outputValues(4, i0, 1) =  0.0;
        outputValues(4, i0, 2) =  -1.0 + 4.0*z;
        
        outputValues(5, i0, 0) = -x*w*(1.0-y-z);
        outputValues(5, i0, 1) = -0.5*(1.0-x-z)*(1.0+x-z)*w;
        outputValues(5, i0, 2) =  0.5*y*x*x*w*w + 0.5*y - 1.0+z;
        
        outputValues(6, i0, 0) =  0.5*(1.0-y-z)*(1.0+y-z)*w;
        outputValues(6, i0, 1) = -y*w*(1.0+x-z);
        outputValues(6, i0, 2) = -0.5*x*y*y*w*w - 0.5*x - 1.0+z; 
        
        outputValues(7, i0, 0) = -x*w*(1.0+y-z);
        outputValues(7, i0, 1) =  0.5*(1.0-x-z)*(1.0+x-z)*w;
        outputValues(7, i0, 2) = -0.5*y*x*x*w*w - 0.5*y - 1.0+z;
        
        outputValues(8, i0, 0) = -0.5*(1.0-y-z)*(1.0+y-z)*w;
        outputValues(8, i0, 1) = -y*w*(1.0-x-z);
        outputValues(8, i0, 2) =  0.5*x*y*y*w*w + 0.5*x - 1.0+z;
        
        outputValues(9, i0, 0) = -(1.0-y-z)*z*w;
        outputValues(9, i0, 1) = -(1.0-x-z)*z*w;
        outputValues(9, i0, 2) =  x*y*w*w + 1.0 - x - y - 2.0*z;
        
        outputValues(10,i0, 0) =  (1.0-y-z)*z*w;
        outputValues(10,i0, 1) = -(1.0+x-z)*z*w;
        outputValues(10,i0, 2) = -x*y*w*w + 1.0 + x - y - 2.0*z;
        
        outputValues(11,i0, 0) =  (1.0+y-z)*z*w;
        outputValues(11,i0, 1) =  (1.0+x-z)*z*w;
        outputValues(11,i0, 2) =  x*y*w*w + 1.0 + x + y - 2.0*z;
        
        outputValues(12,i0, 0) = -(1.0+y-z)*z*w;
        outputValues(12,i0, 1) =  (1.0-x-z)*z*w;
        outputValues(12,i0, 2) = -x*y*w*w + 1.0 - x + y - 2.0*z;
        
      }
      break;
      
    case OPERATOR_CURL:
      TEUCHOS_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_CURL), std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_PYR_I2_FEM): CURL is invalid operator for rank-0 (scalar) functions in 3D");
      break;
      
    case OPERATOR_DIV:
      TEUCHOS_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_DIV), std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_PYR_I2_FEM): DIV is invalid operator for rank-0 (scalar) functions in 3D");
      break;
      
    case OPERATOR_D2:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0,0);
        y = inputPoints(i0,1);
        z = inputPoints(i0,2);

        if(fabs(z-1.0) < eps) {
          if(z <= 1.0) z = 1.0-eps;
          else  z = 1.0+eps;
        }
        Scalar w = 1.0/(1.0 - z);
        
        outputValues(0, i0, 0) = -0.5*(-1.0+y+y*z*w);
        outputValues(0, i0, 1) = -(-0.25 + 0.5*x + 0.5*y + 0.5*z)*w;
        outputValues(0, i0, 2) =  0.25 + (-0.25*y-0.5*x*y-0.25*y*y)*w*w; 
        outputValues(0, i0, 3) = -0.5*(-1.0+x+x*z*w);
        outputValues(0, i0, 4) =  0.25 + (-0.25*x*x-0.25*x-0.5*x*y)*w*w; 
        outputValues(0, i0, 5) =  0.5*x*y*(-1.0-x-y)*w*w*w;
        
        outputValues(1, i0, 0) =  0.5*(1.0-y-y*z*w); 
        outputValues(1, i0, 1) =-(0.25 + 0.5*x - 0.5*y - 0.5*z)*w;
        outputValues(1, i0, 2) = -0.25 + (0.25*y-0.5*x*y+0.25*y*y)*w*w;
        outputValues(1, i0, 3) = -0.5*(-1.0-x-x*z*w); 
        outputValues(1, i0, 4) =  0.25 + (-0.25*x*x + 0.25*x + 0.5*x*y)*w*w; 
        outputValues(1, i0, 5) = -0.5*x*y*(-1.0+x-y)*w*w*w; 
        
        outputValues(2, i0, 0) =  0.5*(1.0+y+y*z*w);
        outputValues(2, i0, 1) =-(-0.25 - 0.5*x - 0.5*y + 0.5*z)*w; 
        outputValues(2, i0, 2) = -0.25 + (-0.25*y+0.5*x*y+0.25*y*y)*w*w;
        outputValues(2, i0, 3) =  0.5*(1.0+x+x*z*w); 
        outputValues(2, i0, 4) = -0.25 + (0.25*x*x -0.25*x + 0.5*x*y)*w*w;  
        outputValues(2, i0, 5) =  0.5*x*y*(-1.0+x+y)*w*w*w;
        
        outputValues(3, i0, 0) = -0.5*(-1.0-y-y*z*w);
        outputValues(3, i0, 1) =-(0.25 - 0.5*x + 0.5*y - 0.5*z)*w; 
        outputValues(3, i0, 2) =  0.25 + (0.25*y+0.5*x*y-0.25*y*y)*w*w; 
        outputValues(3, i0, 3) =  0.5*(1.0-x-x*z*w);
        outputValues(3, i0, 4) = -0.25 + (0.25*x + 0.25*x*x - 0.5*x*y)*w*w;
        outputValues(3, i0, 5) = -0.5*x*y*(-1.0-x+y)*w*w*w;
        
        outputValues(4, i0, 0) =  0.0;
        outputValues(4, i0, 1) =  0.0;
        outputValues(4, i0, 2) =  0.0;
        outputValues(4, i0, 3) =  0.0; 
        outputValues(4, i0, 4) =  0.0;
        outputValues(4, i0, 5) =  4.0; 
        
        outputValues(5, i0, 0) = -(1.0-y-z)*w;
        outputValues(5, i0, 1) =  x*w; 
        outputValues(5, i0, 2) =  x*y*w*w;
        outputValues(5, i0, 3) =  0.0; 
        outputValues(5, i0, 4) =  0.5*x*x*w*w + 0.5; 
        outputValues(5, i0, 5) =  x*x*y*w*w*w + 1.0;
        
        outputValues(6, i0, 0) =  0.0;
        outputValues(6, i0, 1) = -y*w;
        outputValues(6, i0, 2) = -0.5*y*y*w*w - 0.5;
        outputValues(6, i0, 3) =-(1.0+x-z)*w; 
        outputValues(6, i0, 4) = -x*y*w*w; 
        outputValues(6, i0, 5) = -x*y*y*w*w*w + 1.0; 
        
        outputValues(7, i0, 0) = -(1.0+y-z)*w;
        outputValues(7, i0, 1) = -x*w;
        outputValues(7, i0, 2) = -x*y*w*w; 
        outputValues(7, i0, 3) =  0.0; 
        outputValues(7, i0, 4) = -0.5*x*x*w*w - 0.5;
        outputValues(7, i0, 5) = -x*x*y*w*w*w + 1.0; 
        
        outputValues(8, i0, 0) =  0.0;
        outputValues(8, i0, 1) =  y*w;
        outputValues(8, i0, 2) =  0.5*y*y*w*w + 0.5; 
        outputValues(8, i0, 3) = -(1.0-x-z)*w; 
        outputValues(8, i0, 4) =  x*y*w*w;
        outputValues(8, i0, 5) =  x*y*y*w*w*w + 1.0;

        outputValues(9, i0, 0) =  0.0;
        outputValues(9, i0, 1) =  z*w; 
        outputValues(9, i0, 2) =  y*w*w - 1.0;
        outputValues(9, i0, 3) =  0.0; 
        outputValues(9, i0, 4) =  x*w*w - 1.0; 
        outputValues(9, i0, 5) =  2.0*x*y*w*w*w - 2.0; 
         
        outputValues(10,i0, 0) =  0.0;
        outputValues(10,i0, 1) = -z*w;
        outputValues(10,i0, 2) = -y*w*w + 1.0;
        outputValues(10,i0, 3) =  0.0;
        outputValues(10,i0, 4) = -x*w*w - 1.0;
        outputValues(10,i0, 5) = -2.0*x*y*w*w*w - 2.0;
        
        outputValues(11,i0, 0) =  0.0;
        outputValues(11,i0, 1) =  z*w; 
        outputValues(11,i0, 2) =  y*w*w + 1.0;
        outputValues(11,i0, 3) =  0.0;
        outputValues(11,i0, 4) =  x*w*w + 1.0; 
        outputValues(11,i0, 5) =  2.0*x*y*w*w*w - 2.0;      
        
        outputValues(12,i0, 0) =  0.0;
        outputValues(12,i0, 1) = -z*w; 
        outputValues(12,i0, 2) = -y*w*w - 1.0; 
        outputValues(12,i0, 3) =  0.0;
        outputValues(12,i0, 4) = -x*w*w + 1.0;    
        outputValues(12,i0, 5) = -2.0*x*y*w*w*w - 2.0; 
        
      }
      break;
      
    case OPERATOR_D3:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0,0);
        y = inputPoints(i0,1);
        z = inputPoints(i0,2);

        if(fabs(z-1.0) < eps) {
          if(z <= 1.0) z = 1.0-eps;
          else  z = 1.0+eps;
        }

        Scalar w = 1.0/(1.0 - z);

        outputValues(0, i0, 0) =  0.0;
        outputValues(0, i0, 1) = -0.5*w;
        outputValues(0, i0, 2) = -0.5*y*w*w;
        outputValues(0, i0, 3) = -0.5*w;
        outputValues(0, i0, 4) =(-0.25 - 0.5*x - 0.5*y)*w*w;
        outputValues(0, i0, 5) =-(0.5*y + x*y + 0.5*y*y)*w*w*w;
        outputValues(0, i0, 6) =  0.0;
        outputValues(0, i0, 7) = -0.5*x*w*w;
        outputValues(0, i0, 8) =-(0.5*x + x*y + 0.5*x*x)*w*w*w;
        outputValues(0, i0, 9) =  1.5*x*y*(-1.0 - x - y)*w*w*w*w;  
          
        outputValues(1, i0, 0) =  0.0;
        outputValues(1, i0, 1) = -0.5*w;
        outputValues(1, i0, 2) = -0.5*y*w*w;
        outputValues(1, i0, 3) =  0.5*w;
        outputValues(1, i0, 4) = (0.25 - 0.5*x - 0.5*y)*w*w;
        outputValues(1, i0, 5) =-(-0.5*y + x*y - 0.5*y*y)*w*w*w;
        outputValues(1, i0, 6) =  0.0;
        outputValues(1, i0, 7) =  0.5*x*w*w;
        outputValues(1, i0, 8) =-(-0.5*x - x*y + 0.5*x*x)*w*w*w;
        outputValues(1, i0, 9) = -1.5*x*y*(-1.0 + x - y)*w*w*w*w;  
        
        outputValues(2, i0, 0) =  0.0;
        outputValues(2, i0, 1) =  0.5*w;
        outputValues(2, i0, 2) =  0.5*y*w*w;
        outputValues(2, i0, 3) =  0.5*w;
        outputValues(2, i0, 4) = (-0.25 + 0.5*x + 0.5*y)*w*w;
        outputValues(2, i0, 5) =-(0.5*y - x*y - 0.5*y*y)*w*w*w;
        outputValues(2, i0, 6) =  0.0;
        outputValues(2, i0, 7) =  0.5*x*w*w;
        outputValues(2, i0, 8) =-(0.5*x - x*y - 0.5*x*x)*w*w*w;
        outputValues(2, i0, 9) =  1.5*x*y*(-1.0 + x + y)*w*w*w*w;  
        
        outputValues(3, i0, 0) =  0.0;
        outputValues(3, i0, 1) =  0.5*w;
        outputValues(3, i0, 2) =  0.5*y*w*w;
        outputValues(3, i0, 3) = -0.5*w;
        outputValues(3, i0, 4) = (0.25 + 0.5*x - 0.5*y)*w*w;
        outputValues(3, i0, 5) =-(-0.5*y - x*y + 0.5*y*y)*w*w*w;
        outputValues(3, i0, 6) =  0.0;
        outputValues(3, i0, 7) = -0.5*x*w*w;
        outputValues(3, i0, 8) =-(-0.5*x + x*y - 0.5*x*x)*w*w*w;
        outputValues(3, i0, 9) = -1.5*x*y*(-1.0 - x + y)*w*w*w*w;  
        
        outputValues(4, i0, 0) =  0.0;
        outputValues(4, i0, 1) =  0.0;
        outputValues(4, i0, 2) =  0.0;
        outputValues(4, i0, 3) =  0.0;
        outputValues(4, i0, 4) =  0.0;
        outputValues(4, i0, 5) =  0.0;
        outputValues(4, i0, 6) =  0.0;
        outputValues(4, i0, 7) =  0.0;
        outputValues(4, i0, 8) =  0.0;
        outputValues(4, i0, 9) =  0.0;  
        
        outputValues(5, i0, 0) =  0.0;
        outputValues(5, i0, 1) =  w;
        outputValues(5, i0, 2) =  y*w*w;
        outputValues(5, i0, 3) =  0.0;
        outputValues(5, i0, 4) =  x*w*w;
        outputValues(5, i0, 5) =  2.0*y*x*w*w*w;
        outputValues(5, i0, 6) =  0.0;
        outputValues(5, i0, 7) =  0.0;
        outputValues(5, i0, 8) =  x*x*w*w*w;
        outputValues(5, i0, 9) =  3.0*x*x*y*w*w*w*w;  
        
        outputValues(6, i0, 0) =  0.0;
        outputValues(6, i0, 1) =  0.0;
        outputValues(6, i0, 2) =  0.0;
        outputValues(6, i0, 3) = -w;
        outputValues(6, i0, 4) = -y*w*w;
        outputValues(6, i0, 5) = -y*y*w*w*w;
        outputValues(6, i0, 6) =  0.0;
        outputValues(6, i0, 7) = -x*w*w ;
        outputValues(6, i0, 8) = -2.0*x*y*w*w*w;
        outputValues(6, i0, 9) = -3.0*x*y*y*w*w*w*w;  
        
        outputValues(7, i0, 0) =  0.0;
        outputValues(7, i0, 1) = -w;
        outputValues(7, i0, 2) = -y*w*w;
        outputValues(7, i0, 3) =  0.0;
        outputValues(7, i0, 4) = -x*w*w;
        outputValues(7, i0, 5) = -2.0*x*y*w*w*w;
        outputValues(7, i0, 6) =  0.0;
        outputValues(7, i0, 7) =  0.0;
        outputValues(7, i0, 8) = -x*x*w*w*w;
        outputValues(7, i0, 9) = -3.0*x*x*y*w*w*w*w;  
        
        outputValues(8, i0, 0) =  0.0;
        outputValues(8, i0, 1) =  0.0;
        outputValues(8, i0, 2) =  0.0;
        outputValues(8, i0, 3) =  w;
        outputValues(8, i0, 4) =  y*w*w;
        outputValues(8, i0, 5) =  y*y*w*w*w;
        outputValues(8, i0, 6) =  0.0;
        outputValues(8, i0, 7) =  x*w*w;
        outputValues(8, i0, 8) =  2.0*x*y*w*w*w;
        outputValues(8, i0, 9) =  3.0*x*y*y*w*w*w*w ;  
        
        outputValues(9, i0, 0) =  0.0;
        outputValues(9, i0, 1) =  0.0;
        outputValues(9, i0, 2) =  0.0;
        outputValues(9, i0, 3) =  0.0;
        outputValues(9, i0, 4) =  w*w;
        outputValues(9, i0, 5) =  2.0*y*w*w*w;
        outputValues(9, i0, 6) =  0.0;
        outputValues(9, i0, 7) =  0.0;
        outputValues(9, i0, 8) =  2.0*x*w*w*w;
        outputValues(9, i0, 9) =  6.0*x*y*w*w*w*w;  
        
        outputValues(10,i0, 0) =  0.0;
        outputValues(10,i0, 1) =  0.0;
        outputValues(10,i0, 2) =  0.0;
        outputValues(10,i0, 3) =  0.0;
        outputValues(10,i0, 4) = -w*w;
        outputValues(10,i0, 5) = -2.0*y*w*w*w;
        outputValues(10,i0, 6) =  0.0;
        outputValues(10,i0, 7) =  0.0;
        outputValues(10,i0, 8) = -2.0*x*w*w*w;
        outputValues(10,i0, 9) = -6.0*x*y*w*w*w*w;  
        
        outputValues(11,i0, 0) =  0.0;
        outputValues(11,i0, 1) =  0.0;
        outputValues(11,i0, 2) =  0.0;
        outputValues(11,i0, 3) =  0.0;
        outputValues(11,i0, 4) =  w*w;
        outputValues(11,i0, 5) =  2.0*y*w*w*w;
        outputValues(11,i0, 6) =  0.0;
        outputValues(11,i0, 7) =  0.0;
        outputValues(11,i0, 8) =  2.0*x*w*w*w;
        outputValues(11,i0, 9) =  6.0*x*y*w*w*w*w;  
        
        outputValues(12,i0, 0) =  0.0;
        outputValues(12,i0, 1) =  0.0;
        outputValues(12,i0, 2) =  0.0;
        outputValues(12,i0, 3) =  0.0;
        outputValues(12,i0, 4) = -w*w;
        outputValues(12,i0, 5) = -2.0*y*w*w*w;
        outputValues(12,i0, 6) =  0.0;
        outputValues(12,i0, 7) =  0.0;
        outputValues(12,i0, 8) = -2.0*x*w*w*w;
        outputValues(12,i0, 9) = -6.0*x*y*w*w*w*w;  
        
      }
      break;
      
    case OPERATOR_D4:
    case OPERATOR_D5:
    case OPERATOR_D6:
    case OPERATOR_D7:
    case OPERATOR_D8:
    case OPERATOR_D9:
    case OPERATOR_D10:
      {
        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, DkCardinality)
        int DkCardinality = Intrepid::getDkCardinality(operatorType, 
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
      TEUCHOS_TEST_FOR_EXCEPTION( !( Intrepid::isValidOperator(operatorType) ), std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_PYR_I2_FEM): Invalid operator type");
  }
}


  
template<class Scalar, class ArrayScalar>
void Basis_HGRAD_PYR_I2_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar&       outputValues,
                                                             const ArrayScalar &    inputPoints,
                                                             const ArrayScalar &    cellVertices,
                                                             const EOperator        operatorType) const {
  TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                      ">>> ERROR (Basis_HGRAD_PYR_I2_FEM): FEM Basis calling an FVD member function");
}
}// namespace Intrepid
#endif
