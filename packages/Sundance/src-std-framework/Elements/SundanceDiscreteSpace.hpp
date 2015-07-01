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

#ifndef SUNDANCE_DISCRETESPACE_H
#define SUNDANCE_DISCRETESPACE_H

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceSpectralBasis.hpp"
#include "SundanceDOFMapBase.hpp"
#include "PlayaVectorType.hpp"
#include "PlayaVectorSpaceDecl.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "SundanceDiscreteSpaceTransfBuilder.hpp"


namespace Sundance
{
class FunctionSupportResolver;
}

namespace Sundance
{
  using namespace Teuchos;
  using namespace Playa;

  /** 
   * DiscreteSpace represents a discrete finite-element space (i.e., 
   * a mesh and a basis).
   */
  class DiscreteSpace
  {
  public:
    /** */
    DiscreteSpace(){;}
    /** */
    DiscreteSpace(const Mesh& mesh, const BasisFamily& basis,
      const VectorType<double>& vecType,
      int setupVerb = 0);
    /** */
    DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
                  const VectorType<double>& vecType,
      int setupVerb = 0);

    /** */
    DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
                  const Array<CellFilter>& regions,
                  const VectorType<double>& vecType,
      int setupVerb = 0);


    /** */
    DiscreteSpace(const Mesh& mesh, const BasisFamily& basis,
                  const CellFilter& regions,
                  const VectorType<double>& vecType,
      int setupVerb = 0);


    /** */
    DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
                  const CellFilter& regions,
                  const VectorType<double>& vecType,
      int setupVerb = 0);

    /** */
    DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
                  const RCP<DOFMapBase>& map,
                  const VectorType<double>& vecType,
      int setupVerb = 0);

    /** */
    DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
      const RCP<DOFMapBase>& map,
      const RCP<Array<int> >& bcIndices,
      const VectorType<double>& vecType,
      int setupVerb = 0);

    /** */
    DiscreteSpace(const Mesh& mesh, const BasisFamily& basis,
                  const SpectralBasis& spBasis,
                  const VectorType<double>& vecType,
      int setupVerb = 0);
    /** */
    DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
                  const SpectralBasis& spBasis,
                  const VectorType<double>& vecType,
      int setupVerb = 0);

    /** */
    DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
      const RCP<FunctionSupportResolver>& fsr,
      const VectorType<double>& vecType,
      int setupVerb = 0);

    /** */
    const RCP<DOFMapBase>& map() const {return map_;}

    /** return the number of functions */
    int nFunc() const {return basis_.size();}

    /** */
    const BasisArray& basis() const {return basis_;}

    /** */
    Array<std::pair<int,int> > dimStructure() const {return vectorDimStructure(basis());}

    /** */
    Vector<double> createVector() const {return vecSpace_.createMember();}

    /** */
    VectorSpace<double> vecSpace() const {return vecSpace_;}

    /** */
    const Mesh& mesh() const {return mesh_;}

    /** */
    const VectorType<double>& vecType() const {return vecType_;}

    /** */
    void importGhosts(const Vector<double>& x,
                      RCP<GhostView<double> >& ghostView) const ;

    /** */
    void getAllowedFuncs(const CellFilter& cf, Set<int>& funcs) const ;

    /** */
    const CellFilter& cellFilters(int i) const {return subdomains_[i];}

    /** Return the transformation builder */
    const RCP<DiscreteSpaceTransfBuilder>& getTransformation() const
    		{ return transformationBuilder_; }

  private:

    /** */
    void init(const Array<CellFilter>& regions,
              const BasisArray& basis);

    /** */
    void init(const Array<CellFilter>& regions,
      const BasisArray& basis,
      const RCP<Array<int> >& isBCIndex, 
      bool partitionBCs);

    /** */
    Array<CellFilter> maximalRegions(int n) const ;

    /** */
    void initVectorSpace(
      const RCP<Array<int> >& isBCIndex, 
      bool partitionBCs);
    
    /** */
    void initImporter();

    /** */
    int setupVerb_;

    /** */
    RCP<DOFMapBase> map_;

    /** */
    Mesh mesh_;

    /** */
    Array<CellFilter> subdomains_;

    /** */
    BasisArray basis_;

    /** */
    VectorSpace<double> vecSpace_;

    /** */
    VectorType<double> vecType_;

    /** */
    RCP<GhostImporter<double> > ghostImporter_;

    /** Transformation builder in case when it is needed*/
    RCP<DiscreteSpaceTransfBuilder> transformationBuilder_;

  };

}



#endif
