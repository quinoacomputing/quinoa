// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef ABSTRACT_LIN_ALG_PACK_TYPES_H
#define ABSTRACT_LIN_ALG_PACK_TYPES_H

#include <memory>
#include <stdexcept>

#include "RTOp.h"
#include "RTOpPack_OldTypes.hpp"
#include "DenseLinAlgPack_Types.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_RCP.hpp"

namespace AbstractLinAlgPack {

#include "DenseLinAlgPack_PublicTypes.ud"

typedef RTOp_index_type  size_type;
typedef RTOp_value_type  value_type;
typedef RTOp_index_type  index_type;
using Teuchos::RCP;

#ifdef DOXYGEN_COMPILE // Doxygen needs a little help finding these links
/** \brief . */
typedef DenseLinAlgPack::VectorTmpl<value_type>             DVector;
/** \brief . */
typedef DenseLinAlgPack::VectorSliceTmpl<value_type>        DVectorSlice;
/** \brief . */
typedef DenseLinAlgPack::DMatrix                            DMatrix;
/** \brief . */
typedef DenseLinAlgPack::DMatrixSlice                       DMatrixSlice;
/** \brief . */
typedef DenseLinAlgPack::DMatrixSliceTriEle                 DMatrixSliceTriEle;
/** \brief . */
typedef DenseLinAlgPack::DMatrixSliceTri                    DMatrixSliceTri;
/** \brief . */
typedef DenseLinAlgPack::DMatrixSliceSym                    DMatrixSliceSym;
/** \brief . */
typedef RangePack::Range1D Range1D;
#endif


/** @name Exception classes */
//@{

/// Base class for precondition exceptions
class PreConditionException : public std::logic_error
{public: PreConditionException(const std::string& what_arg) : std::logic_error(what_arg) {}};

/// Base class for postcondition exceptions
class PostConditionException : public std::runtime_error
{public: PostConditionException(const std::string& what_arg) : std::runtime_error(what_arg) {}};

/// Base class for input exceptions (Preconditions).
class InputException : public PreConditionException
{public: InputException(const std::string& what_arg) : PreConditionException(what_arg) {}};

/// Base class for invalid setup for a class object when an exception is thrown
class SetupException : public PreConditionException
{public: SetupException(const std::string& what_arg) : PreConditionException(what_arg) {}};


//@}

/** @name Main interface library */
//@{

class InnerProduct;

class VectorSpaceFactory;
class VectorSpace;
class Vector;
class VectorMutable;

class MatrixBase;
class MatrixOp;
class MatrixNonsing;
class MatrixOpNonsing;
class MatrixSymOp;
class MatrixSymNonsing;
class MatrixSymOpNonsing;
class MatrixSymDiag;

class MultiVector;
class MultiVectorMutable;

class MatrixSymSecant;

class BasisSystem;
class BasisSystemPerm;
class BasisSystemFactory;

class Permutation;

// template classes

template <class T_Indice, class T_Value>	class SparseElement;
template <class T_Element, class T_Alloc>	class SparseVector;
template <class T_Element>					class SparseVectorSlice;

// concrete classes

class EtaVector;
class GenPermMatrixSlice;
typedef SparseVector<
  SparseElement<index_type,value_type>
  , std::allocator<
    SparseElement<index_type,value_type>
    >
  >												SpVector;
typedef SparseVectorSlice<
  SparseElement<index_type,value_type> >		SpVectorSlice;

//@}

/** @name Standard tools library */
//@{

// pure abstract classes

class PermVector;

class MatrixSymInitDiag;
class MatrixSymDiag;

// concrete subclasses

class BasisSystemComposite;
class VectorSpaceBlocked;
class VectorMutableBlocked;
class MatrixOpSubView;
class MatrixComposite;
class MatrixSymIdent;
class MatrixSymDiagStd;
class MatrixZero;
class MatrixPermAggr;
class MatrixOpNonsingAggr;

// testing classes

class VectorSpaceTester;
class VectorSpaceTesterSetOptions;
class MatrixOpNonsingTester;
class BasisSystemTester;
class BasisSystemTesterSetOptions;

//@}

/** @name Serial interface library */
//@{

// pure abstract classes

class MatrixOpSerial;
class MatrixNonsingSerial;
class MatrixSymOpSerial;
class MatrixSymNonsingSerial;
class MatrixOpNonsingSerial;
class MatrixSymOpNonsingSerial;
class MatrixSymDenseInitialize;
class MatrixSymDiagSparse;
class MatrixLoadSparseElements;
class MatrixConvertToSparse;
class MatrixExtractSparseElements;
class MatrixExtractInvCholFactor;
class MatrixSymOpGetGMSSymMutable;
class MatrixSymOpGetGMSSym;
class MatrixSymAddDelUpdateable;

//@}

/** @name Serial implementations library */
//@{

class VectorDenseEncap;
class VectorDenseMutableEncap;
class MatrixDenseEncap;
class MatrixDenseMutableEncap;
class MatrixDenseSymEncap;
class MatrixDenseSymMutableEncap;
class MatrixDenseTriEncap;

class PermutationSerial;
class VectorSpaceSerial;
class VectorMutableDense;
class VectorSparse;
class MatrixSparseCOORSerial;
class MatrixSymPosDefCholFactor;
class MatrixConvertToSparseEncap;
class MultiVectorMutableDense;

class MatrixSymDiagSparseStd;

//@}

/** @name Serial solvers library */
//@{

// Matrix scaling classes

class MatrixScaling_Strategy;

// Sparse linear solver classes

class DirectSparseSolver;        // Abstract interface
class DirectSparseSolverImp;     // Node implementation classs
class DirectSparseSolverMA28;    // Concrete subclass
class DirectSparseSolverMA48;    // ""
class DirectSparseSolverSuperLU; // ""

// BasisSystem classes

class BasisSystemPermDirectSparse;
class BasisSystemFactoryStd;

//@}

} // end namespace AbstractLinAlgPack

#endif // ABSTRACT_LIN_ALG_PACK_TYPES_H
