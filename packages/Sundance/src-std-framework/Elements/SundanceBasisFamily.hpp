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

#ifndef SUNDANCE_BASISFAMILY_H
#define SUNDANCE_BASISFAMILY_H

#include "SundanceDefs.hpp"
#include "SundanceBasisFamilyBase.hpp"
#include "SundanceOrderedHandle.hpp"
#include "Teuchos_Array.hpp"

namespace Sundance
{
using Teuchos::XMLObject;
using Teuchos::tuple;
using Teuchos::Array;

class CommonFuncDataStub;

/** 
 * BasisFamily is the user-level handle class for specifying the basis with
 * which a test, unknown, or discrete function is represented.  Basis
 * functions can be vector-valued, as is the case with, for example, the
 * Nedelec basis in electromagnetics; the dim() method returns the spatial
 * dimension of the basis functions. Scalar-valued bases naturally have
 * dim()=1.

*/
class BasisFamily : public OrderedHandle<BasisFamilyBase>
{
public:
  /* handle ctor boilerplate */
  ORDERED_HANDLE_CTORS(BasisFamily, BasisFamilyBase);

  /** write to XML */
  XMLObject toXML() const ;

  /** 
   * \brief Return the polynomial order of the basis functions, for use in
   * determining the quadrature rule required to integrate a product of
   * basis functions exactly. The polynomial order will be the smallest
   * integer for which all mixed partial derivatives vanish exactly.
   *  
   * Note: in H(div) and H(curl) spaces the order of accuracy is not
   * always an integer, and the relationship between the order of accuracy
   * and the return value of the order() method is not necessarily simple
   * (for instance, it can depend on things such as the convexity of the
   * boundary). Thus it is better to think of this method as specifying
   * the required order of quadrature, and not the order of accuracy of
   * approximation nor the order to which the space is complete.
   */
  int order() const ;

  /** \brief
   * Return the number of DOFs for this basis on the given 
   * reference cell type, including its facets.
   */
  int nReferenceDOFsWithFacets(const CellType& maximalCellType,
    const CellType& cellType) const ;

  /** \brief
   * Return the number of DOFs for this basis on the given 
   * reference cell type, not including those on facets.
   */
  int nReferenceDOFsWithoutFacets(const CellType& maximalCellType,
    const CellType& cellType) const ;

  /** 
   * \brief Return the dimension of the members of 
   * a vector-valued basis. Return 1 if the basis
   * is scalar-valued. Otherwise, return the spatial dimension.
   */
  int dim() const ;

  /** 
   * \brief Return the tensor order of the members of 
   * a basis.
   */
  int tensorOrder() const {return ptr()->tensorOrder();}

  /** */
  bool operator==(const BasisFamily& other) const ;

  

  /** \brief Inform caller as to whether I am a scalar basis. */
  bool isScalarBasis() const {return ptr()->isScalarBasis();}

  /** \brief Inform caller as to whether I am an H(div) basis. */
  bool isHDivBasis() const {return ptr()->isHDivBasis();}

  /** \brief Inform caller as to whether I am an H(curl) basis. */
  bool isHCurlBasis() const {return ptr()->isHCurlBasis();}

  /** Sum up the dim() values for array of bases. */
  static int size(const Array<BasisFamily>& b) ;

  /** Extract the basis from an expression */
  static BasisFamily getBasis(const RCP<const CommonFuncDataStub>& funcData);

  /** Extract the basis from an expression */
  static RCP<BasisDOFTopologyBase> getBasisTopology(const RCP<const CommonFuncDataStub>& funcData);

  /** */
  void refEval(
    const CellType& cellType,
    const Array<Point>& pts,
    const SpatialDerivSpecifier& deriv,
    Array<Array<Array<double> > >& result,
    int verbosity) const ;

  /**  */
  void getConstrainsForHNDoF( const int indexInParent,
			      const int maxCellDim,
			      const int maxNrChild,
			      const int facetDim,
			      const int facetIndex,
			      const int nodeIndex,
			      Array<int>& localDoFs,
			      Array<int>& parentFacetDim,
			      Array<int>& parentFacetIndex,
			      Array<int>& parentFacetNode,
			      Array<double>& coefs
			      ) const;



  /** \brief Inform caller whether basis requires an element transformation */
  bool requiresBasisTransformation() const { return ptr()->requiresBasisTransformation(); }


  /** */
  virtual void preApplyTransformation( const CellType &maxCellType ,
				       const Mesh &mesh,
				       const Array<int> &cellLIDs,
				       const CellJacobianBatch& JVol,
				       RCP<Array<double> >& A
				       ) const 
  { 
    ptr()->preApplyTransformation( maxCellType ,
				   mesh,
				   cellLIDs,
				   JVol,
				   A);
  }
  /**  */
  virtual void postApplyTransformation( const CellType &maxCellType ,
					const Mesh &mesh,
					const Array<int> &cellLIDs,
					const CellJacobianBatch& JVol,
					RCP<Array<double> >& A ) const
  {
    ptr()->postApplyTransformation( maxCellType ,
				    mesh ,
				    cellLIDs ,
				    JVol,
				    A);
  }
  
  virtual void preApplyTransformationTranspose( const CellType &maxCellType ,
						const Mesh &mesh,
						const Array<int> &cellLIDs,
						const CellJacobianBatch& JVol,
						Array<double>& A ) const
  {
    ptr()->preApplyTransformationTranspose( maxCellType ,
					    mesh ,
					    cellLIDs ,
					    JVol ,
					    A );
  }

};

/** \relates BasisFamily */
inline Array<BasisFamily> List(const BasisFamily& a, const BasisFamily& b)
{
  return Array<BasisFamily>(tuple(a,b));
}

/** \relates BasisFamily */
inline Array<BasisFamily> List(const BasisFamily& a, const BasisFamily& b,
  const BasisFamily& c)
{
  return tuple(a,b,c);
}

/** \relates BasisFamily */
inline Array<BasisFamily> List(const BasisFamily& a, const BasisFamily& b,
  const BasisFamily& c, const BasisFamily& d)
{
  return tuple(a,b,c,d);
}

/** \relates BasisFamily */
inline Array<BasisFamily> List(const BasisFamily& a, const BasisFamily& b,
  const BasisFamily& c, const BasisFamily& d, 
  const BasisFamily& e)
{
  return tuple(a,b,c,d,e);
}


/** \relates BasisFamily */
inline Array<BasisFamily> List(const BasisFamily& a, const BasisFamily& b,
  const BasisFamily& c, const BasisFamily& d, 
  const BasisFamily& e, const BasisFamily& f)
{
  return tuple(a,b,c,d,e,f);
}

/** \relates BasisFamily */
inline Array<BasisFamily> List(const BasisFamily& a, const BasisFamily& b,
  const BasisFamily& c, const BasisFamily& d, 
  const BasisFamily& e, const BasisFamily& f, 
  const BasisFamily& g)
{
  return tuple(a,b,c,d,e,f,g);
}

/** \relates BasisFamily */
inline Array<BasisFamily> List(const BasisFamily& a, const BasisFamily& b,
  const BasisFamily& c, const BasisFamily& d, 
  const BasisFamily& e, const BasisFamily& f, 
  const BasisFamily& g, const BasisFamily& h)
{
  return tuple(a,b,c,d,e,f,g,h);
}

/** \relates BasisFamily */
inline Array<BasisFamily> List(const BasisFamily& a, const BasisFamily& b,
  const BasisFamily& c, const BasisFamily& d, 
  const BasisFamily& e, const BasisFamily& f, 
  const BasisFamily& g, const BasisFamily& h, 
  const BasisFamily& i)
{
  return tuple(a,b,c,d,e,f,g,h,i);
}

/** \relates BasisFamily */
inline Array<BasisFamily> List(const BasisFamily& a, const BasisFamily& b,
  const BasisFamily& c, const BasisFamily& d, 
  const BasisFamily& e, const BasisFamily& f, 
  const BasisFamily& g, const BasisFamily& h, 
  const BasisFamily& i, const BasisFamily& j)
{
  return tuple(a,b,c,d,e,f,g,h,i,j);
}

/** */
inline Array<BasisFamily> replicate(const BasisFamily& b, int n)
{
  Array<BasisFamily> rtn(n);
  for (int i=0; i<n; i++) rtn[i] = b;
  return rtn;
}


/** */
inline Array<BasisFamily> replicate(const Array<BasisFamily>& b, int n)
{
  Array<BasisFamily> rtn(n*b.size());
  for (int i=0; i<n*b.size(); i++) rtn[i] = b[0];
  return rtn;
}

class BasisArray : public Array<BasisFamily>
{
public:
  BasisArray() : Array<BasisFamily>() {;}

  BasisArray(int n) : Array<BasisFamily>(n) {;}

  BasisArray(const Array<BasisFamily>& a) 
    : Array<BasisFamily>(a) 
    {;}

};


/** \relates BasisFamily */
Array<std::pair<int, int> > vectorDimStructure(const Array<BasisFamily>& basis);


/** \relates BasisFamily */
Array<std::pair<int, int> > vectorDimStructure(const BasisFamily& basis);


/** \relates BasisFamily 
 * Indicate whether members of a basis have support on a boundary only
 * if their associated dofs live on the boundary. This will return true
 * for Lagrange, false for P1NC, RT, Nedelec, and most other bases. 
 *
 * This is used to simplify boundary integrals.
 */
bool basisRestrictableToBoundary(const BasisFamily& b);


}

#endif
