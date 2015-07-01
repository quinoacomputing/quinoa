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

#ifndef SUNDANCE_CELLPREDICATE_H
#define SUNDANCE_CELLPREDICATE_H



#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceCellPredicateBase.hpp"
#include "SundancePositionalCellPredicate.hpp"
#include "PlayaHandle.hpp"
#include "SundanceParametrizedCurve.hpp"
#include "SundanceCellCurvePredicate.hpp"

namespace Sundance
{
using namespace Teuchos;
  
  
/** 
 * User-level handle for predicates (deriving from CellPredicateBase)
 * used to decide whether
 * a given cell passes a CellFilter.
 */
class CellPredicate : public Playa::Handle<CellPredicateBase>
{
public:
    
  /* handle boilerplate */
  HANDLE_CTORS(CellPredicate, CellPredicateBase);

  /** construct from a positional cell predicate functor */
  CellPredicate(const RCP<CellPredicateFunctorBase>& func);

  /** construct from a positional cell predicate functor */
  CellPredicate(Playa::Handleable<CellPredicateFunctorBase>* func);

  /** construct from a positional cell predicate functor */
  CellPredicate(ParametrizedCurve& curve, CurveCellFilterMode filterMode);

  /** write to XML */
  XMLObject toXML() const {return ptr()->toXML();}

  /** */
  std::string description() const {return ptr()->description();}



  /** set the mesh on which cells are to be tested */
  void setMesh(const Mesh& mesh, int cellDim) const 
    {ptr()->setMesh(mesh, cellDim);}

  /** compare to another predicate, used for placement in STL containers */
  bool operator<(const CellPredicate& other) const ;


};

}

namespace std
{
inline ostream& operator<<(std::ostream& os, const Sundance::CellPredicate& pred)
{
  os << pred.toXML() << std::endl;
  return os;
}
}


#endif
