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


#ifndef SUNDANCE_COORDINATESYSTEM_HPP
#define SUNDANCE_COORDINATESYSTEM_HPP

#include "SundanceExpr.hpp"
#include "SundanceCoordExpr.hpp"
#include "PlayaExceptions.hpp"
#include <cmath>


namespace Sundance
{
using Sundance::Expr;
using namespace Teuchos;

class CoordinateSystemBase : 
  public Playa::Handleable<CoordinateSystemBase>
{
public:
  /** */
  virtual ~CoordinateSystemBase(){}

  /** */
  virtual Expr jacobian() const = 0 ;

  /** */
  virtual bool supportsMeshDimension(int dim) const = 0 ;

  /** */
  virtual std::ostream& print(std::ostream& os) const = 0 ;

  /** */
  double pi() const {return 4.0*::atan(1.0);}
};


class CoordinateSystem : public Playa::Handle<CoordinateSystemBase>
{
public:
  /* boilerplate handle ctors */
  HANDLE_CTORS(CoordinateSystem, CoordinateSystemBase);

  /** */
  Expr jacobian() const {return ptr()->jacobian();}

  /** */
  bool supportsMeshDimension(int dim) const 
    {return ptr()->supportsMeshDimension(dim);}

};


class CartesianCoordinateSystem : public CoordinateSystemBase
{
public:
  /** */
  CartesianCoordinateSystem(){}

  /** */
  Expr jacobian() const {return Expr(1.0);}

  /** */
  bool supportsMeshDimension(int dim) const {return dim > 0;}

  /** */
  std::ostream& print(std::ostream& os) const 
    {
      os << "Cartesian";
      return os;
    };

  
  /** */
  virtual RCP<CoordinateSystemBase> getRcp() {return rcp(this);}
};


class MeridionalCylindricalCoordinateSystem : public CoordinateSystemBase
{
public:
  /** */
  MeridionalCylindricalCoordinateSystem()
    : r_(new CoordExpr(0)) {}

  /** */
  Expr jacobian() const {return pi()*r_;}

  /** */
  bool supportsMeshDimension(int dim) const {return dim > 0 && dim <= 2;}

  /** */
  std::ostream& print(std::ostream& os) const 
    {
      os << "Meridional Cylindrical";
      return os;
    };
  
  /** */
  virtual RCP<CoordinateSystemBase> getRcp() {return rcp(this);}

private:
  Expr r_;
};

class RadialSphericalCoordinateSystem : public CoordinateSystemBase
{
public:
  /** */
  RadialSphericalCoordinateSystem()
    : r_(new CoordExpr(0)) {}

  /** */
  Expr jacobian() const {return pi()*r_*r_;}

  /** */
  bool supportsMeshDimension(int dim) const {return dim == 1;}
  
  /** */
  std::ostream& print(std::ostream& os) const 
    {
      os << "Radial Spherical";
      return os;
    };



  
  /** */
  virtual RCP<CoordinateSystemBase> getRcp() {return rcp(this);}
private:
  Expr r_;
};


class CoordinateSystemBuilder
{
public:
  /** */
  static CoordinateSystem makeCoordinateSystem(const std::string& name);
};


inline std::ostream& operator<<(std::ostream& os, const CoordinateSystem& cs)
{
  return cs.ptr()->print(os);
}

}




#endif
