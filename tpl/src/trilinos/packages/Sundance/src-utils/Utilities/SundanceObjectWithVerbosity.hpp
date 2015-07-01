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

#ifndef SUNDANCE_OBJECTWITHVERBOSITY_H
#define SUNDANCE_OBJECTWITHVERBOSITY_H

#include "SundanceDefs.hpp"
#include "PlayaObjectWithVerbosity.hpp"
#include "SundanceParamUtils.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Sundance
{

using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::FancyOStream;
using Teuchos::ParameterEntry;


/** 
 * Defines traits for getting default verbosity parameters from a class
 */
template <class X>
class VerbosityTraits
{
public:
  static RCP<ParameterList> defaultVerbParams()
    {return X::defaultVerbParams();}
};



/**
 * ObjectWithClassVerbosity and the related verbosity() method
 * provide a method for getting/setting
 * verbosity flags for entire classes.
 *
 * You can set verbosity for a single instance of a class, or for
 * the whole class. To set for an instance, use the verbosity()
 * member function, for example,
 * \code
 * Mesh mesh1 = reader1.getMesh();
 * Mesh mesh2 = reader2.getMesh();
 * Mesh mesh3 = reader3.getMesh();
 * mesh1.verbosity() = 3;
 * \endcode
 * which sets the verbosity of <tt>mesh1</tt> to 3 and leaves
 * those of <tt>mesh2</tt> and <tt>mesh3</tt> unchanged.
 *
 * Alternatively, you can set a default verbosity for an entire
 * class, for example,
 * \code
 * Mesh mesh1 = reader1.getMesh();
 * Mesh mesh2 = reader2.getMesh();
 * Mesh mesh3 = reader3.getMesh();
 * mesh1.verbosity() = 3;
 * verbosity<Mesh>() = 2;
 * \endcode
 * which sets the default verbosity to 2. Since <tt>mesh1</tt>
 * has its own verbosity setting of 3, 
 * it will use it rather than the
 * default, but <tt>mesh2</tt> and <tt>mesh3</tt> will use 2.
 * 
 */
template <class X>
class ObjectWithClassVerbosity : public Playa::ObjectWithVerbosity
{
public:
  /** \deprecated Construct, starting silent */
  ObjectWithClassVerbosity(int verb=classVerbosity())
    : Playa::ObjectWithVerbosity(verb) {;}

  /** \deprecated Writeable access to the default verbosity for the class */
  static int& classVerbosity() 
    {
      static int rtn = 0;
      return rtn;
    }

};



/** 
 * \relates ObjectWithClassVerbosity
 * Global method for setting verbosity of a class
 */
template <class X> int& verbosity() 
{
  return X::classVerbosity();
}

template <class X>
class ParameterControlledObjectWithVerbosity 
  : public ObjectWithClassVerbosity<X>
{
public:
  /** \deprecated Construct, starting silent */
  ParameterControlledObjectWithVerbosity() : ObjectWithClassVerbosity<X>() {;}

  /** Construct with a parameter list controlling the verbosity settings */
  ParameterControlledObjectWithVerbosity(const std::string& objName, const ParameterList& p)
    : ObjectWithClassVerbosity<X>(),
      verbControlParams_() 
    {
      RCP<ParameterList> defaults = VerbosityTraits<X>::defaultVerbParams();
      TEUCHOS_TEST_FOR_EXCEPTION(defaults->name() != objName, std::runtime_error,
        "mismatched ParameterList names for verbosity control: expected "
        << defaults->name() << ", got " << objName);
      TEUCHOS_TEST_FOR_EXCEPTION(defaults->name() != p.name(), std::runtime_error,
        "mismatched ParameterList names for verbosity control: expected "
        << defaults->name() << ", got " << p.name());
      verbControlParams_ = rcp(new ParameterList(mergeParams(*defaults, p)));
    }

  /** */
  int verbLevel(const std::string& context) const 
    {
      const ParameterEntry* pe = verbControlParams_->getEntryPtr(context);
      TEUCHOS_TEST_FOR_EXCEPTION(pe==0, std::runtime_error,
        "parameter with name \"" << context << "\" not found in verbosity "
        "control parameter list " << verbControlParams_);
      TEUCHOS_TEST_FOR_EXCEPTION(pe->getAny().type() != typeid(int),
        std::runtime_error,
        "context parameter name \"" 
        << context << "\" does not have an integer value in verbosity "
        "control parameter list " << verbControlParams_);

      return Teuchos::any_cast<int>(pe->getAny());
    }

  /** */
  const ParameterList& verbSublist(const std::string& name) const 
    {
      TEUCHOS_TEST_FOR_EXCEPTION(!verbControlParams_->isSublist(name),
        std::runtime_error,
        "context parameter name \"" 
        << name << "\" is not a sublist in verbosity "
        "control parameter list " << *verbControlParams_);

      return verbControlParams_->sublist(name);
    }

  /** */
  ParameterList mergeParams(const ParameterList& pDef, const ParameterList& pIn) const
    {
      return mergeParamLists(pDef, pIn);
    }
       
  /** */
  const ParameterList params() const {return *verbControlParams_;}

  /** */
  RCP<ParameterList> modifiableParams() const {return verbControlParams_;}
private:
  RCP<ParameterList> verbControlParams_;
};



}





#endif
