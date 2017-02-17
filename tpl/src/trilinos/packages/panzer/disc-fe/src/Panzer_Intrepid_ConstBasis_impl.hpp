#ifndef __Panzer_Intrepid_ConstBasis_impl_hpp__
#define __Panzer_Intrepid_ConstBasis_impl_hpp__

namespace panzer {
  
  
template<class Scalar, class ArrayScalar>
Basis_Constant<Scalar,ArrayScalar>::Basis_Constant(const shards::CellTopology & ct)
  {
    this -> basisCardinality_  = 1;
    this -> basisDegree_       = 0;    
    this -> basisCellTopology_ = ct;
    this -> basisType_         = Intrepid2::BASIS_FEM_DEFAULT;
    this -> basisCoordinates_  = Intrepid2::COORDINATES_CARTESIAN;
    this -> basisTagsAreSet_   = false;

    initializeTags();
    this->basisTagsAreSet_ = true;
  }
  
  
  
template<class Scalar, class ArrayScalar>
void Basis_Constant<Scalar, ArrayScalar>::initializeTags() {
  
  // Basis-dependent initializations
  int tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
  int posScDim = 0;        // position in the tag, counting from 0, of the subcell dim 
  int posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
  int posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell

  int dim = this -> basisCellTopology_.getDimension();
  
  // An array with local DoF tags assigned to the basis functions, in the order of their local enumeration 
  int tags[]  = { dim, 0, 0, 1 };
  
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
void Basis_Constant<Scalar, ArrayScalar>::getValues(ArrayScalar &        outputValues,
                                                            const ArrayScalar &  inputPoints,
                                                            const Intrepid2::EOperator      operatorType) const {
  
  // Verify arguments
#ifdef HAVE_INTREPID_DEBUG
  Intrepid2::getValues_HGRAD_Args<Scalar, ArrayScalar>(outputValues,
                                                      inputPoints,
                                                      operatorType,
                                                      this -> getBaseCellTopology(),
                                                      this -> getCardinality() );
#endif
  
  // Number of evaluation points = dim 0 of inputPoints
  int dim0 = inputPoints.dimension(0);  
  
  switch (operatorType) {
    
    case Intrepid2::OPERATOR_VALUE:
      for (int i0 = 0; i0 < dim0; i0++) {
        // outputValues is a rank-2 array with dimensions (basisCardinality_, dim0)
        outputValues(0, i0) = 1.0;
      }
      break;

    case Intrepid2::OPERATOR_GRAD:
    case Intrepid2::OPERATOR_D1:
    case Intrepid2::OPERATOR_CURL:
    case Intrepid2::OPERATOR_DIV:
    case Intrepid2::OPERATOR_D2:
    case Intrepid2::OPERATOR_D3:
    case Intrepid2::OPERATOR_D4:
    case Intrepid2::OPERATOR_D5:
    case Intrepid2::OPERATOR_D6:
    case Intrepid2::OPERATOR_D7:
    case Intrepid2::OPERATOR_D8:
    case Intrepid2::OPERATOR_D9:
    case Intrepid2::OPERATOR_D10:
    default:
      TEUCHOS_TEST_FOR_EXCEPTION( !( Intrepid2::isValidOperator(operatorType) ), std::invalid_argument,
                          ">>> ERROR (Basis_Constant): Invalid operator type");
  }
}
  

  
template<class Scalar, class ArrayScalar>
void Basis_Constant<Scalar, ArrayScalar>::getValues(ArrayScalar&           outputValues,
                                                            const ArrayScalar &    inputPoints,
                                                            const ArrayScalar &    cellVertices,
                                                            const Intrepid2::EOperator        operatorType) const {
  TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                      ">>> ERROR (Basis_Constant): FEM Basis calling an FVD member function");
}

}// namespace panzer
#endif
