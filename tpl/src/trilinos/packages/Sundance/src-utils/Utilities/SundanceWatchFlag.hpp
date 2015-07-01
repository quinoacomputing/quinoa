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

#ifndef SUNDANCE_WATCHFLAG_H
#define SUNDANCE_WATCHFLAG_H

#include "SundanceDefs.hpp"
#include "SundanceMap.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_ParameterList.hpp"
#include <string>



namespace Sundance
{
using Teuchos::ParameterList;
  /** 
   * Class WatchFlag is used to tag individual expressions for
   * increased verbosity in certain tasks. The tasks that can be marked
   * are listed below. Verbosity values are integers, with 0 being silent, 
   * and typically 5 being excruciatingly detailed.
   * 
<ul>
<li> evaluation
<li> evaluator setup
<li> discrete function evaluation
<li> symbolic preprocessing
<li> equation set setup
<li> assembler setup
<li> assembly loop
<li> dof map setup
<li> dof map access
<li> eval mediator
<li> integration
<li> integral transformation
<li> fill
<li> matrix config
<li> vector config
<li> solve control
</ul>
* The verbosity level for a task is set using the setParam() function.
* The watch can be turned on or off collectively with the deactivate function.
   */
  class WatchFlag
    {
    public:
      /** */
      WatchFlag(const std::string& name="", 
        const ParameterList& params = *defaultParams());

      /** */
      const std::string& name() const {return name_;}

      /** */
      void activate() ;

      /** */
      void deactivate() ;

      /** */
      bool isActive() const ;

      /** */
      bool operator<(const WatchFlag& other) const
        {return name() < other.name();}

      /** */
      XMLObject toXML() const ;

      /** */
      int param(const std::string& name) const ;

      /** */
      void setParam(const std::string& name, int val);

      /** */
      static RCP<ParameterList> defaultParams();


    private:
      std::string name_;

      RCP<ParameterList> params_;

      static Map<std::string, bool>& isActiveMap()
        {
          static Map<std::string, bool> rtn;
          return rtn;
        }

    };
}

#endif
