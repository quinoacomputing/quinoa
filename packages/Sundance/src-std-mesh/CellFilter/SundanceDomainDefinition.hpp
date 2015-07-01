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


/*
 * SundanceDomainDefinition.hpp
 *
 *  Created on: Apr 30, 2010
 *      Author: benk
 */

#ifndef SUNDANCEDOMAINDEFINITION_HPP_
#define SUNDANCEDOMAINDEFINITION_HPP_

#include "SundanceDefs.hpp"
#include "PlayaHandleable.hpp"
#include "PlayaHandle.hpp"
#include "SundancePoint.hpp"

#include "SundanceParametrizedCurve.hpp"

namespace Sundance {

using namespace Teuchos;

/** define the predicate also with estimation */
#define MESH_DOMAIN_(name, code) \
  class name : public MeshDomainBase, \
               public Playa::Handleable<MeshDomainBase> \
  { \
  public:\
    name() : MeshDomainBase(){;}            \
	virtual bool isInsideComputationalDomain(const Point& x) const code \
    GET_RCP(MeshDomainBase);\
  }

#define MESH_DOMAIN(name, code) MESH_DOMAIN_(name, code);


/** Base class for mesh refinement , but also to define computational domain
 * (which must not be the same as the mesh domain)*/
class MeshDomainBase  {

public:

	MeshDomainBase() {;}

	virtual ~MeshDomainBase() {;}

	/** Function to answer if the domain is inside the computational domain <br>
	 * The strategy should be that if one point of a cell is inside the domain, then the whole
	 * cell should be considered as in the computational domain.
	 * @param x [in] coordinate of the point
	 * @return true if the point is inside the domain, false otherwise */
	virtual bool isInsideComputationalDomain(const Point& x) const { return true;}
};

// ---------------

/**  Class defines mesh domain based on parametrized curve */
class CurveDomain : public MeshDomainBase ,
                    public Playa::Handleable<MeshDomainBase>{
public:

	/** Ctor with the 2 necessary input arguments */
	CurveDomain(const ParametrizedCurve& curve ,
			    CurveCellFilterMode mode) :
		MeshDomainBase() , curve_(curve) , mode_(mode){;}
    /** empty Dtor */
	virtual ~CurveDomain() {;}

	/**  in or outside the domain */
	virtual bool isInsideComputationalDomain(const Point& x) const {
		if (mode_ == Outside_Curve){
            return (curve_.curveEquation(x) >= -1e-16 );
		}
		else{
			return (curve_.curveEquation(x) <= 1e-16 );
		}
	}

	GET_RCP(MeshDomainBase);

private:

	const ParametrizedCurve& curve_;

	const CurveCellFilterMode mode_;
};

// ---------------

class MeshDomainDef : public Playa::Handle<MeshDomainBase> {
public:

	/* Handle constructors */
	HANDLE_CTORS(MeshDomainDef, MeshDomainBase);

	/** see MeshDomainBase for Docu */
	bool isInsideComputationalDomain(const Point& x) const {
		return ptr()->isInsideComputationalDomain(x);
	}

private:

};


}

#endif /* SUNDANCEDOMAINDEFINITION_HPP_ */
