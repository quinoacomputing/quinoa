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


#ifndef PDEOPT_POINTDATA_H
#define PDEOPT_POINTDATA_H

#include "Sundance.hpp"
#include <fstream>


namespace Sundance
{
/**
 * PointData is a utility class that helps create the Sundance objects required
 * to include point measurements in an inversion problem. 
 */
class PointData 
{
public:
  /** Create a point data object given a set of sensor locations and
   * corresponding sensor readings */
  PointData(const Array<Point>& locations, const Array<double>& values,
    const double& pointComparisonTolerance);

  /** Create a point data object given a set of sensor locations and
   * a field to probe at those points */
  PointData(const Array<Point>& locations, const Expr& field,
    const double& pointComparisonTolerance);

  /** Read point data from an XML file */
  PointData(const XMLObject& xml, const Mesh& mesh);

    

  /** Return an expr that, when called at one of the measurement locations,
   * returns the value of the measurement there. */
  Expr sensorValues() const {return sensorVals_;}

  /** Return a cell filter that identifies the sensor points */
  CellFilter sensorLocations() const {return sensorLocations_;}


  /** Adjust a set of points to the nearest vertices on a mesh */
  static Array<Point> snapToMesh(const Mesh& mesh, const Array<Point>& locations) ;

  /** Find the vertex nearest to a specified point */
  static Point nearestMeshPoint(const Mesh& mesh, const Point& x);
    
    
private:
  void init(const Array<Point>& locations, const Array<double>& values,
    const double& pointComparisonTolerance);

  Expr sensorVals_;
  CellFilter sensorLocations_;
};



/** 
 * Do lexigraphic comparison on points, including a bit of sloppiness 
 */
class SloppyPointComparitor : public std::binary_function<Point, Point, bool>
{
public: 
  SloppyPointComparitor(const double& tol) : tol_(tol) {;}

  bool operator()(const Point& p1, const Point& p2) const ;
private:
  double tol_;
};


/** 
 * This is a functor that allows the sensor readings to be obtained through
 * a Sundance user-defined expression.
 */
class PointDataExprFunctor : public PointwiseUserDefFunctor0
{
public:
  /** */
  PointDataExprFunctor(const Array<Point>& locations, 
    const Array<double>& values,
    const double& pointComparisonTolerance);

  /** */
  virtual void eval0(const double* in, double* out) const ;

  /** */
  int numArgs() const {return dim_;}

private:
  std::map<Point, double, SloppyPointComparitor> pointToValueMap_;
  int dim_;
};
  
/** 
 * This is a functor that allows a positional cell predicate to test
 * against the locations of the sensors.
 */
class PointDataCellPredicateFunctor : public CellPredicateFunctorBase,
                                      public Playa::Handleable<CellPredicateFunctorBase>
{
public:
  /** */
  PointDataCellPredicateFunctor(const Array<Point>& locations, 
    const double& pointComparisonTolerance);
    
  /** */
  virtual bool operator()(const Point& x) const ;

  GET_RCP(CellPredicateFunctorBase);

private:
  std::set<Point, SloppyPointComparitor> pointSet_;
};
  

}

#endif
