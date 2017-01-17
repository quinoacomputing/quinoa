#ifndef INTREPID2_HGRAD_TRI_CN_FEM_ORTHDEF_HPP
#define INTREPID2_HGRAD_TRI_CN_FEM_ORTHDEF_HPP
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

/** \file   Intrepid_HGRAD_TRI_Cn_FEM_ORTHDef.hpp
    \brief  Definition file for FEM orthogonal basis functions of arbitrary degree 
            for H(grad) functions on TRI.
    \author Created by R. Kirby
*/

namespace Intrepid2 {
  
template<class Scalar, class ArrayScalar>
Basis_HGRAD_TRI_Cn_FEM_ORTH<Scalar,ArrayScalar>::Basis_HGRAD_TRI_Cn_FEM_ORTH( int degree )
{
    this -> basisCardinality_  = (degree+1)*(degree+2)/2;
    this -> basisDegree_       = degree;
    this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Triangle<3> >() );
    this -> basisType_         = BASIS_FEM_HIERARCHICAL;
    this -> basisCoordinates_  = COORDINATES_CARTESIAN;
    this -> basisTagsAreSet_   = false;

    initializeTags();
    this->basisTagsAreSet_ = true;
}
  
  
  
template<class Scalar, class ArrayScalar>
void Basis_HGRAD_TRI_Cn_FEM_ORTH<Scalar, ArrayScalar>::initializeTags() {
  
  // Basis-dependent initializations
  int tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
  int posScDim = 0;        // position in the tag, counting from 0, of the subcell dim 
  int posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
  int posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell
  
  // An array with local DoF tags assigned to the basis functions, in the order of their local enumeration 
  int *tags = new int[tagSize * this->getCardinality()];
  for (int i=0;i<this->getCardinality();i++) {
    tags[4*i] = 2;
    tags[4*i+1] = 0;
    tags[4*i+2] = i;
    tags[4*i+3] = this->getCardinality();
  }
  
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
void Basis_HGRAD_TRI_Cn_FEM_ORTH<Scalar, ArrayScalar>::getValues(ArrayScalar &        outputValues,
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
  const int deg = this->getDegree();

  // add more here and put in appropriate extra case statements below to enable higher derivatives.
  void (*tabulators[])(ArrayScalar &, const int, const ArrayScalar &)
    = { TabulatorTri<Scalar,ArrayScalar,0>::tabulate ,
        TabulatorTri<Scalar,ArrayScalar,1>::tabulate ,
        TabulatorTri<Scalar,ArrayScalar,2>::tabulate };


  switch (operatorType) {
  case OPERATOR_VALUE:
    tabulators[0]( outputValues , deg , inputPoints );
    break;
  case OPERATOR_GRAD:
    tabulators[1]( outputValues , deg , inputPoints );
    break;
  case OPERATOR_D1:
  case OPERATOR_D2:
    // add more case OPEATOR_Dn statements if you've added more items to the
    // array above.
    tabulators[operatorType-OPERATOR_D1+1]( outputValues , deg , inputPoints );
    break;
  default:
    TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument,
                        ">>> ERROR (Basis_HGRAD_TRI_Cn_FEM_ORTH): invalid or unsupported operator" );

  }

  return;
}



template<class Scalar, class ArrayScalar>
void Basis_HGRAD_TRI_Cn_FEM_ORTH<Scalar, ArrayScalar>::getValues(ArrayScalar&           outputValues,
                                                                 const ArrayScalar &    inputPoints,
                                                                 const ArrayScalar &    cellVertices,
                                                                 const EOperator        operatorType) const {
  TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                      ">>> ERROR (Basis_HGRAD_TRI_Cn_FEM): FEM Basis calling an FVD member function");
}



template<typename Scalar, typename ArrayScalar>
void TabulatorTri<Scalar,ArrayScalar,0>::tabulate(ArrayScalar &outputValues ,
                                                  const int deg ,
                                                  const ArrayScalar &z )
{
  const int np = z.dimension( 0 );
  
  // each point needs to be transformed from Pavel's element
  // z(i,0) --> (2.0 * z(i,0) - 1.0)
  // z(i,1) --> (2.0 * z(i,1) - 1.0)
  
  // set up constant term
  int idx_cur = TabulatorTri<Scalar,ArrayScalar,0>::idx(0,0);
  int idx_curp1,idx_curm1;
  
  // set D^{0,0} = 1.0
  for (int i=0;i<np;i++) {
    outputValues(idx_cur,i) = Scalar( 1.0 ) + z(i,0) - z(i,0) + z(i,1) - z(i,1);
  }
  

  if (deg > 0) {
    Teuchos::Array<Scalar> f1(np),f2(np),f3(np);
    
    for (int i=0;i<np;i++) {
      f1[i] = 0.5 * (1.0+2.0*(2.0*z(i,0)-1.0)+(2.0*z(i,1)-1.0));
      f2[i] = 0.5 * (1.0-(2.0*z(i,1)-1.0));
      f3[i] = f2[i] * f2[i];
    }
    
    // set D^{1,0} = f1
    idx_cur = TabulatorTri<Scalar,ArrayScalar,0>::idx(1,0);
    for (int i=0;i<np;i++) {
      outputValues(idx_cur,i) = f1[i];
    }
    
    // recurrence in p
    for (int p=1;p<deg;p++) {
      idx_cur = TabulatorTri<Scalar,ArrayScalar,0>::idx(p,0);
      idx_curp1 = TabulatorTri<Scalar,ArrayScalar,0>::idx(p+1,0);
      idx_curm1 = TabulatorTri<Scalar,ArrayScalar,0>::idx(p-1,0);
      Scalar a = (2.0*p+1.0)/(1.0+p);
      Scalar b = p / (p+1.0);
      
      for (int i=0;i<np;i++) {
        outputValues(idx_curp1,i) = a * f1[i] * outputValues(idx_cur,i)
          - b * f3[i] * outputValues(idx_curm1,i);
      }
    }
    
    // D^{p,1}
    for (int p=0;p<deg;p++) {
      int idxp0 = TabulatorTri<Scalar,ArrayScalar,0>::idx(p,0);
      int idxp1 = TabulatorTri<Scalar,ArrayScalar,0>::idx(p,1);
      for (int i=0;i<np;i++) {
        outputValues(idxp1,i) = outputValues(idxp0,i)
          *0.5*(1.0+2.0*p+(3.0+2.0*p)*(2.0*z(i,1)-1.0));
      }
    }
    
    
    // recurrence in q
    for (int p=0;p<deg-1;p++) {
      for (int q=1;q<deg-p;q++) {
        int idxpqp1=TabulatorTri<Scalar,ArrayScalar,0>::idx(p,q+1);
        int idxpq=TabulatorTri<Scalar,ArrayScalar,0>::idx(p,q);
        int idxpqm1=TabulatorTri<Scalar,ArrayScalar,0>::idx(p,q-1);
        Scalar a,b,c;
        TabulatorTri<Scalar,ArrayScalar,0>::jrc((Scalar)(2*p+1),(Scalar)0,q,a,b,c);
        for (int i=0;i<np;i++) {
          outputValues(idxpqp1,i)
            = (a*(2.0*z(i,1)-1.0)+b)*outputValues(idxpq,i)
            - c*outputValues(idxpqm1,i);
        }
      }
    }
  }    
  
  // orthogonalize
  for (int p=0;p<=deg;p++) {
    for (int q=0;q<=deg-p;q++) {
      for (int i=0;i<np;i++) {
        outputValues(TabulatorTri<Scalar,ArrayScalar,0>::idx(p,q),i) *= sqrt( (p+0.5)*(p+q+1.0));
      }
    }
  }
  
  return;
}
 

 
template<typename Scalar, typename ArrayScalar>
void TabulatorTri<Scalar,ArrayScalar,1>::tabulate(ArrayScalar &outputValues ,
                                                  const int deg ,
                                                  const ArrayScalar &z ) 
{
  const int np = z.dimension(0);
  const int card = outputValues.dimension(0);
  FieldContainer<Sacado::Fad::SFad<Scalar,2> > dZ( z.dimension(0) , z.dimension(1) );
  for (int i=0;i<np;i++) {
    for (int j=0;j<2;j++) {
      dZ(i,j) = Sacado::Fad::SFad<Scalar,2>( z(i,j) );
      dZ(i,j).diff(j,2);
    }
  }
  FieldContainer<Sacado::Fad::SFad<Scalar,2> > dResult(card,np);

  TabulatorTri<Sacado::Fad::SFad<Scalar,2>,FieldContainer<Sacado::Fad::SFad<Scalar,2> >,0>::tabulate( dResult ,
                                                                                                  deg ,
                                                                                                  dZ );

  for (int i=0;i<card;i++) {
    for (int j=0;j<np;j++) {
      for (int k=0;k<2;k++) {
        outputValues(i,j,k) = dResult(i,j).dx(k);
      }
    }
  }

  return;

}



template<typename Scalar, typename ArrayScalar, unsigned derivOrder>
void TabulatorTri<Scalar,ArrayScalar,derivOrder>::tabulate( ArrayScalar &outputValues ,
                                                            const int deg ,
                                                            const ArrayScalar &z ) 
{
  const int np = z.dimension(0);
  const int card = outputValues.dimension(0);
  FieldContainer<Sacado::Fad::SFad<Scalar,2> > dZ( z.dimension(0) , z.dimension(1) );
  for (int i=0;i<np;i++) {
    for (int j=0;j<2;j++) {
      dZ(i,j) = Sacado::Fad::SFad<Scalar,2>( z(i,j) );
      dZ(i,j).diff(j,2);
    }
  }
  FieldContainer<Sacado::Fad::SFad<Scalar,2> > dResult(card,np,derivOrder+1);

  TabulatorTri<Sacado::Fad::SFad<Scalar,2>,FieldContainer<Sacado::Fad::SFad<Scalar,2> >,derivOrder-1>::tabulate(dResult ,
                                                                                                            deg ,
                                                                                                            dZ );

  for (int i=0;i<card;i++) {
    for (int j=0;j<np;j++) {
      outputValues(i,j,0) = dResult(i,j,0).dx(0);
      for (unsigned k=0;k<derivOrder;k++) {
        outputValues(i,j,k+1) = dResult(i,j,k).dx(1);
      }
    }
  }

  return;


}


}// namespace Intrepid2
#endif
