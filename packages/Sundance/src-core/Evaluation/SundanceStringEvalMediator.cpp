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


#include "SundanceStringEvalMediator.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceCellDiameterExpr.hpp"
#include "SundanceCellVectorExpr.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"
#include "PlayaExceptions.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;

using namespace Teuchos;
using namespace std;



StringEvalMediator::StringEvalMediator()
  : AbstractEvalMediator() 
{}

void StringEvalMediator::evalCoordExpr(const CoordExpr* expr,
                                       RCP<EvalVector>& vec) const
{
  SUNDANCE_MSG1(verb(), 
    "evaluating coord expr " << expr->toString());
  
  vec->resize(1);
  vec->start()[0]=0.0;
  vec->setString(expr->name());
}

void StringEvalMediator::evalCellDiameterExpr(const CellDiameterExpr* expr,
                                              RCP<EvalVector>& vec) const
{
  SUNDANCE_MSG1(verb(), 
    "evaluating cell diameter expr " << expr->toXML().toString());
  
  vec->resize(1);
  vec->start()[0]=0.0;
  vec->setString(expr->name());
}

void StringEvalMediator::evalCellVectorExpr(const CellVectorExpr* expr,
                                              RCP<EvalVector>& vec) const
{
  SUNDANCE_MSG1(verb(), "evaluating cell vector expr " << expr->toXML().toString());
  
  vec->resize(1);
  vec->start()[0]=0.0;
  vec->setString(expr->name());
}

void StringEvalMediator
::evalDiscreteFuncElement(const DiscreteFuncElement* expr,
                          const Array<MultiIndex>& mi,
                          Array<RCP<EvalVector> >& vec) const 
{
  static Array<string> coordNames;

  if (coordNames.size() != 3)
    {
      coordNames.resize(3);
      coordNames[0] = "x";
      coordNames[1] = "y";
      coordNames[2] = "z";
    }

  std::string funcName = expr->name();
  
  for (int i=0; i<mi.size(); i++)
    {
      vec[i]->resize(1);
      vec[i]->start()[0]=0.0;
      if (mi[i].order()==0)
        {
          vec[i]->setString(funcName);
        }
      else
        {
          int dir = mi[i].firstOrderDirection();
          std::string deriv = "D[" + funcName + ", " + coordNames[dir] + "]";
          vec[i]->setString(deriv);
        }
    }
}
