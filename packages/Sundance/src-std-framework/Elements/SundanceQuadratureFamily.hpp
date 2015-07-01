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

#ifndef SUNDANCE_QUADRATUREFAMILY_H
#define SUNDANCE_QUADRATUREFAMILY_H

#include "SundanceDefs.hpp"
#include "SundanceQuadratureFamilyBase.hpp"
#include "PlayaHandle.hpp"

namespace Sundance
{

/** 
 * QuadratureFamily is a geometry-independent specification of
 * a method by which quadrature is to be carried out. For example,
 * a GaussianQuadrature family will generate Gaussian
 * quadrature points on any cell type.
 */
class QuadratureFamily : public Playa::Handle<QuadratureFamilyStub>
{
public:
  /* */
  HANDLE_CTORS(QuadratureFamily, QuadratureFamilyStub);
  /** */
  XMLObject toXML() const ;

  /** */
  int order() const ;

  /** Returns the number of points in a rule of the given cell type 
      WARNING: this is slow.  Call it once and store the result. 
      TODO: make it pure virtual and override with queries in
      the derived classes, making them supply the information.  */
  int getNumPoints( const CellType& cellType ) const;


  /** Get the quadrature points and weights for the given cell type */
  void getPoints(const CellType& cellType, 
    Array<Point>& quadPoints,
    Array<double>& quadWeights) const ;

  /** Get quadrature points and weights for integration on a facet of a cell */
  void getFacetPoints(const CellType& cellType, 
    int facetDim,
    int facetIndex,
    Array<Point>& quadPoints,
    Array<double>& quadWeights) const ;

  /** Get the quadrature points and weights for the given cell type ,
   * which might be cut by a curve in the case of, see class for more info */
  void getAdaptedWeights(const CellType& cellType ,
  		         int cellDim,
	             int celLID,
	             int facetIndex ,
                 const Mesh& mesh ,
                 const ParametrizedCurve& globalCurve ,
                 Array<Point>& quadPoints ,
                 Array<double>& quadWeights ,
                 bool& isCut) const ;

private:
  /** Get quad points for a facet of a line */
  void getLineFacetQuad(int facetDim,
    int facetIndex,
    Array<Point>& quadPoints,
    Array<double>& quadWeights) const ;
  /** Get quad points for a facet of a triangle */
  void getTriangleFacetQuad(int facetDim,
    int facetIndex,
    Array<Point>& quadPoints,
    Array<double>& quadWeights) const ;

  /** Get quad points for a facet of a quadlateral */
  void getQuadFacetQuad(int facetDim,
    int facetIndex,
    Array<Point>& quadPoints,
    Array<double>& quadWeights) const ;

  /** Get quad points for a facet of a tet */
  void getTetFacetQuad(int facetDim,
    int facetIndex,
    Array<Point>& quadPoints,
    Array<double>& quadWeights) const ;

  /** Get quad points for a facet of a Brick cell */
  void getBrickFacetQuad(int facetDim,
    int facetIndex,
    Array<Point>& quadPoints,
    Array<double>& quadWeights) const ;
};


/** \relates QuadratureFamily */
void printQuad(std::ostream& os, 
  const Array<Point>& pts, const Array<double>& wgts);


}

#endif
