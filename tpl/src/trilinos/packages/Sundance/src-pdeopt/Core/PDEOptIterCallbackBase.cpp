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


#include "PDEOptIterCallbackBase.hpp"
#include "PDEOptPDEConstrainedObjBase.hpp"
#include "Sundance.hpp"

namespace Sundance
{

DefaultIterCallback::DefaultIterCallback(
  const std::string& filename, 
  const std::string& type,
  int frequency)
  : type_(type), filename_(filename), frequency_(frequency)
{}

void DefaultIterCallback::call(const PDEConstrainedObjBase* obj, 
  int iter) const
{
  if (iter % frequency_ != 0) return;
 
  string name = filename_ + "-iter-" + Teuchos::toString(iter);
  
  FieldWriter writer;

  if (type_=="VTK")
  {
    writer = new VTKWriter(name);
  }
  else if (type_=="Exodus")
  {
    writer = new ExodusWriter(name);
  }
  else if (type_=="Matlab")
  {
    writer = new MatlabWriter(name);
  }
  else 
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, RuntimeError, 
      "writer type [" << type_ << "] not defined");
  }

  Array<Expr> state = obj->stateVars();
  Array<Expr> adjoint = obj->adjointVars();
  Expr design = obj->designVar();

  writer.addMesh(obj->mesh());

  for (int b=0; b<state.size(); b++)
  {
    for (int i=0; i<state[b].size(); i++)
    {
      string tag = "[" + Teuchos::toString(b)
        + "][" + Teuchos::toString(i) + "]";
      writer.addField("state" + tag, new ExprFieldWrapper(state[b][i]));
      writer.addField("adjoint" + tag, new ExprFieldWrapper(adjoint[b][i]));
    }
  }

  for (int i=0; i<design.size(); i++)
  {
    string tag = "[" + Teuchos::toString(i) + "]";
    writer.addField("design" + tag, new ExprFieldWrapper(design[i]));
  }
  
  writer.write();
}


}
