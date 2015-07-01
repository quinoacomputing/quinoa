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

#ifndef SUNDANCE_DISCRETEFUNCTIONDATA_H
#define SUNDANCE_DISCRETEFUNCTIONDATA_H

#include "SundanceDefs.hpp"
#include "SundanceDiscreteFuncDataStub.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceDiscreteSpace.hpp"
#include "PlayaVectorDecl.hpp"

namespace Sundance
{
using namespace Teuchos;

/** 
 * DiscreteFunctionData 
 */
class DiscreteFunctionData : public DiscreteFuncDataStub
{
public:
  /** */
  DiscreteFunctionData(const DiscreteSpace& space);

  /** */
  DiscreteFunctionData(const DiscreteSpace& space, 
    const Playa::Vector<double>& vec);

  /** */
  DiscreteFunctionData(const DiscreteSpace& space, const double& constantValue);

  /** virtual destructor */
  virtual ~DiscreteFunctionData() {;}

  /** */
  void updateGhosts() const ;

  /** */
  void setVector(const Vector<double>& vec);

  /** */
  const Vector<double>& getVector() const {return vector_;}

  /** */
  const DiscreteSpace& discreteSpace() const {return space_;}

  /** */
  const Mesh& mesh() const {return space_.mesh();}

  /** */
  const RCP<DOFMapBase>& map() const {return space_.map();}

  /** */
  RCP<const MapStructure> getLocalValues(int cellDim, 
    const Array<int>& cellLID,
    Array<Array<double> >& localValues) const ;


  /** */
  RCP<GhostView<double> > ghostView() const 
    {updateGhosts(); return ghostView_;}

  /** */
  const BasisArray& basis() const {return space_.basis();}

  /** */
  static const DiscreteFunctionData* getData(const DiscreteFuncElement* ufe);

  /** */
  static DiscreteFunctionData* getData(DiscreteFuncElement* ufe);


private:

  DiscreteSpace space_;

  Vector<double> vector_;

  mutable RCP<GhostView<double> > ghostView_;

  mutable bool ghostsAreValid_;

};
}



#endif
