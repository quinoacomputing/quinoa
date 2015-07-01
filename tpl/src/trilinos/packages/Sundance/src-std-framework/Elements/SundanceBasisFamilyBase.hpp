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

#ifndef SUNDANCE_BASISFAMILYBASE_H
#define SUNDANCE_BASISFAMILYBASE_H

#include "SundanceCellJacobianBatch.hpp"
#include "SundanceDefs.hpp"
#include "SundanceBasisDOFTopologyBase.hpp"
#include "SundanceTensorBasisBase.hpp"
#include "SundanceBasisReferenceEvaluationBase.hpp"
#include "PlayaHandleable.hpp"
#include "SundanceMesh.hpp"
#include "PlayaPrintable.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "SundanceTypeUtils.hpp"
#include "Teuchos_XMLObject.hpp"

namespace Sundance {

using namespace Teuchos;


/** 
 *
 */
class BasisFamilyBase
  : public Playa::Handleable<BasisFamilyBase>,
    public Playa::Printable,
    public ObjectWithClassVerbosity<BasisFamilyBase>,
    public BasisDOFTopologyBase,
    public TensorBasisBase,
    public BasisReferenceEvaluationBase
{
public:

  /** */
  virtual int order() const = 0 ;

  /** */
  virtual bool lessThan(const BasisDOFTopologyBase* other) const ;

  /** \brief Indicates whether mapping the basis requires an additional
      correction */
  virtual bool requiresBasisTransformation() const { return false; }
  /** \brief Default transformation does nothing */
  virtual void preApplyTransformation( const CellType &maxCellType ,
				       const Mesh &mesh,
				       const Array<int> &cellLIDs,
				       const CellJacobianBatch& JVol,
				       RCP<Array<double> >& A
				       ) const {;}
  /** \brief Default transformation does nothing */
  virtual void postApplyTransformation( const CellType &maxCellType ,
				       const Mesh &mesh,
				       const Array<int> &cellLIDs,
					const CellJacobianBatch& JVol,
					RCP<Array<double> >& A
					) const {;}
  /** \brief Default transformation does nothing */
  virtual void preApplyTransformationTranspose( const CellType &maxCellType ,
						const Mesh &mesh,
						const Array<int> &cellLIDs,
						const CellJacobianBatch& JVol,
						Array<double>& A ) const {;}

};


/* ----- Subtypes of BasisFamilyBase that specify transformation type --- */

/**     
 * Base class for scalar-valued basis families. Bases for scalar
 * fields living in, e.g., H1 or L2, should derive from this class. 
 */
class ScalarBasis 
  : public BasisFamilyBase
{
public:
  /** Inform caller that my tensor order is zero */
  int tensorOrder() const {return 0;}

  /** Inform caller that I am a scalar basis */
  bool isScalarBasis() const {return true;}
  
   /** Return the dimension of the members of a scalar-valued basis  */
  int dim() const {return 1;}
};

/** */
class VectorBasis
  : public BasisFamilyBase
{
public:
  /** */
  VectorBasis(int dim) : dim_(dim) {}
  /** Inform caller that my tensor order is one */
  int tensorOrder() const {return 1;}

   /** Return the dimension of the members of a scalar-valued basis  */
  int dim() const {return dim_;}

private:
  int dim_;
};


/** 
 * Base class for bases living in H(div) 
 */
class HDivVectorBasis
  : public VectorBasis
{
public:
  /** */
  HDivVectorBasis(int dim) : VectorBasis(dim) {}
  
  /** Inform caller that I am an H(div) basis */
  bool isHDivBasis() const {return true;}

};

/** 
 * Base class for bases living in H(curl) 
 */
class HCurlVectorBasis
  : public VectorBasis
{
public:
  /** */
  HCurlVectorBasis(int dim) : VectorBasis(dim) {}

  /** Inform caller that I am an H(curl) basis */
  bool isHCurlBasis() const {return true;}

};




} // namespace Sundance


#endif
