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

#ifndef SUNDANCE_ELEMENTINTEGRAL_H
#define SUNDANCE_ELEMENTINTEGRAL_H

#include "SundanceDefs.hpp"
#include "SundanceCellJacobianBatch.hpp"
#include "SundanceQuadratureFamily.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceParametrizedCurve.hpp"
#include "SundanceMesh.hpp"
#include "Teuchos_Array.hpp"


namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;

/** 
 * ElementIntegral encapsulates the common data needed for the
 * integration of groups of related zero-forms, one-forms, and two-forms.
 */
class ElementIntegral
{
public:
  /** Construct a zero-form */
  ElementIntegral(int spatialDim,
    const CellType& maxCellType,
    int dim, 
    const CellType& cellType,
    bool isInternalBdry,
    const ParametrizedCurve& globalCurve,
    const Mesh& mesh,
    int verb);

  /** Construct a one-form */
  ElementIntegral(int spatialDim,
    const CellType& maxCellType,
    int dim, 
    const CellType& cellType,
    const BasisFamily& testBasis,
    int alpha,
    int testDerivOrder,
    bool isInternalBdry,
    const ParametrizedCurve& globalCurve,
    const Mesh& mesh,
    int verb);

  /** Construct a two-form */
  ElementIntegral(int spatialDim,
    const CellType& maxCellType,
    int dim,
    const CellType& cellType,
    const BasisFamily& testBasis,
    int alpha,
    int testDerivOrder,
    const BasisFamily& unkBasis,
    int beta,
    int unkDerivOrder,
    bool isInternalBdry,
    const ParametrizedCurve& globalCurve,
    const Mesh& mesh,
    int verb);

  /** virtual dtor */
  virtual ~ElementIntegral(){;}

  /** Indicate whether this element integral is a 
   *  zero, one, or two form */
  int order() const {return order_;}
      
  /** Return the number of nodes associated with the test function */
  int nNodesTest() const {return nNodesTest_;}
      
  /** Return the number of nodes associated with the test function */
  int nNodesUnk() const {return nNodesUnk_;}

  /** Return the total number of elements in this local stiffness
   * matrix */
  int nNodes() const {return nNodes_;}

  /** Return the number of different facets for which integrals 
   * must be tabulated in the cases where an integral must be done 
   * by referring back to a maximal cell */
  int nFacetCases() const {return nFacetCases_;}

  /** Whether this is an integral on an internal boundary */
  bool isInternalBdry() const {return isInternalBdry_;}

  /** */
  static bool& alwaysUseCofacets();

  /** */
  void setVerb(
    int integrationVerb,
    int transformVerb);

  /** */
  int setupVerb() const {return setupVerb_;}
    
  /** */
  int integrationVerb() const {return integrationVerb_;}
    
  /** */
  int transformVerb() const {return transformVerb_;}

  /** */
  void describe(std::ostream& os) const ;

  /** */
  static int& transformationMatrixIsValid(int alpha);


  /** */
  static int& transformationMatrixIsValid(int alpha, int beta);

  /** */
  static void invalidateTransformationMatrices();

  /** */
  static double& totalFlops() {static double rtn = 0; return rtn;}

  BasisFamily &getTestBasis() { return testBasis_; }
  BasisFamily &getUnknownBasis() { return unkBasis_; }

protected:

  /** */
  void assertBilinearForm() const ;

  /** */
  void assertLinearForm() const ;

  /** */
  static void addFlops(const double& flops) {totalFlops() += flops;}

  /** The dimension of the cell being integrated */
  int dim() const {return dim_;}

  /** The dimension of the space in which the cell is embedded */
  int spatialDim() const {return spatialDim_;}
      
  /** Number of test function derivatives wrt reference coordinates that
   * are needed to evaluate this integral. Will always be equal to 
   * ipow(element dimension, differentiation order). */
  int nRefDerivTest() const {return nRefDerivTest_;}
      
  /** Number of unknown function derivatives wrt reference coordinates that
   * are needed to evaluate this integral. Will always be equal to 
   * ipow(element dimension, differentiation order). */
  int nRefDerivUnk() const {return nRefDerivUnk_;}

  /** The order to which the test function is differentiated in this
   * integral. */
  int testDerivOrder() const {return testDerivOrder_;}

  /** The order to which the unknown function is differentiated in this
   * integral. */
  int unkDerivOrder() const {return unkDerivOrder_;}

  /** */
  int alpha() const {return alpha_;}

  /** */
  int beta() const {return beta_;}

  /** */
  const CellType& cellType() const {return cellType_;}

  /** */
  const CellType& maxCellType() const {return maxCellType_;}

  /** */
  const CellType& evalCellType() const {return evalCellType_;}

  /** */
  const BasisFamily& testBasis() const {return testBasis_;}

  /** */
  const BasisFamily& unkBasis() const {return unkBasis_;}

  /** Workspace for element transformations involving one derivative*/
  static Array<double>& G(int gamma) ;

  /** Workspace for element transformations involving two derivatives */
  static Array<double>& G(int gamma, int delta) ;

  /** return base to the given power */
  static int ipow(int base, int power);

  /** The value below which chop() sets numbers to zero */
  static double chopVal() {static double rtn=1.0e-14; return rtn;}

  /** Chop a number to zero if it is smaller in magnitude than
   * the value chopVal() */
  static double chop(const double& x) 
    {
      if (::fabs(x) > chopVal()) return x;
      else return 0.0;
    }

  /** */
  void getQuad(const QuadratureFamily& quad, int evalCase,
    Array<Point>& quadPts, Array<double>& quadWeights) const ;

  /** */
  void createTwoFormTransformationMatrix(const CellJacobianBatch& JTrans,
    const CellJacobianBatch& JVol) const;
  /** */
  void createOneFormTransformationMatrix(const CellJacobianBatch& JTrans,
    const CellJacobianBatch& JVol) const;

  /** */
  const Mesh& mesh() const {return mesh_;}

  /** */
  const ParametrizedCurve& globalCurve() const {return globalCurve_;}

private:

  int setupVerb_;
  int integrationVerb_;
  int transformVerb_;

  int spatialDim_;

  int dim_;

  bool isInternalBdry_;

  int nFacetCases_;

  int testDerivOrder_;

  int nRefDerivTest_;

  int nNodesTest_;

  int unkDerivOrder_;

  int nRefDerivUnk_;

  int nNodesUnk_;

  int nNodes_;

  int order_;

  int alpha_;

  int beta_;

  CellType cellType_;

  CellType maxCellType_;

  CellType evalCellType_;

  BasisFamily testBasis_;

  BasisFamily unkBasis_;
      
  /** The curve which might be used for adaptive integration method */
  const ParametrizedCurve globalCurve_;

  /** The curve which might be used for adaptive integration method */
  const Mesh mesh_;

};
}


#endif
