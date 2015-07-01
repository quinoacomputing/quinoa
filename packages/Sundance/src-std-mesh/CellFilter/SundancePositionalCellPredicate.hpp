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


#ifndef SUNDANCE_POSITIONALCELLPREDICATE_H
#define SUNDANCE_POSITIONALCELLPREDICATE_H


#include "SundanceDefs.hpp"
#include "SundanceCellPredicateBase.hpp"

namespace Sundance
{
using namespace Teuchos;

#define NEW_CELL_PREDICATE(name)  \
  class name : public CellPredicateFunctorBase,  \
               public Playa::Handleable<CellPredicateFunctorBase>  \
  {  \
  public:  \
    name() : CellPredicateFunctorBase(#name) {} \
    virtual ~name() {}  \
    virtual bool operator()(const Point& x) const; \
    GET_RCP(CellPredicateFunctorBase);  \
  }; \
  \
  bool name::operator()(const Point& x) const


#define CELL_PREDICATE_(name, code) \
  class name : public CellPredicateFunctorBase, \
               public Playa::Handleable<CellPredicateFunctorBase> \
  { \
  public:\
    name() : CellPredicateFunctorBase(#name){;}            \
    virtual ~name(){;}\
    virtual bool operator()(const Point& x) const code \
    GET_RCP(CellPredicateFunctorBase);\
  }


#define CELL_PREDICATE(name, code) CELL_PREDICATE_(name, code);




/** */
class CellPredicateFunctorBase
{
public:
  /** */
  CellPredicateFunctorBase(const std::string& name="Functor(" + Teuchos::toString(topID()) + ")")
    : name_(name) {;}

  /** */
  virtual ~CellPredicateFunctorBase(){;}

  /** */
  virtual bool operator()(const Point& x) const = 0 ;

  /** */
  virtual std::string description() const {return name_;}
private:
  static int& topID() {static int rtn=0; rtn++; return rtn;}
  std::string name_;
};

  
/** 
 * PositionalCellPredicate tests whether the cell's nodes satisfy
 * a condition on their positions.
 */
class PositionalCellPredicate : public CellPredicateBase 
{
public:
      
  /** Construct with a function of positions */
  PositionalCellPredicate(const RCP<CellPredicateFunctorBase>& func) 
    : CellPredicateBase(), func_(func) 
    {;}

  /** virtual dtor */
  virtual ~PositionalCellPredicate(){;}
      
  /** Test the predicate on a batch of cells */
  virtual void testBatch(const Array<int>& cellLID,
    Array<int>& results) const ;

  /** Write to XML */
  virtual XMLObject toXML() const ;

  /** comparison */
  virtual bool lessThan(const CellPredicateBase* other) const ;

  /** */
  virtual std::string description() const {return func_->description();}

  /* */
  GET_RCP(CellPredicateBase);

private:
  RCP<CellPredicateFunctorBase> func_;
};





/** */
class PointCellPredicateFunctor 
  : public CellPredicateFunctorBase
{
public:
  /** */
  PointCellPredicateFunctor(const Point& x, const double& tol=1.0e-12)
    : x_(x), tol_(tol){}

  /** */
  bool operator()(const Point& x) const ;

private:
  Point x_;
  double tol_;
};



/** */
class CoordinateValueCellPredicateFunctor 
  : public CellPredicateFunctorBase
{
public:
  /** */
  CoordinateValueCellPredicateFunctor(
    int direction, const double& value, const double& tol=1.0e-12)
    : direction_(direction), value_(value), tol_(tol) {}

  /** */
  bool operator()(const Point& x) const ;

private:
  int direction_;
  double value_;
  double tol_;
};

/** */
class PointCellPredicate : public PositionalCellPredicate
{
public:
  /** */
  PointCellPredicate(const Point& x, const double& tol=1.0e-12)
    : PositionalCellPredicate(rcp(new PointCellPredicateFunctor(x,tol)))
    {}

  /* */
  GET_RCP(CellPredicateBase);
};

/** */
class CoordinateValueCellPredicate : public PositionalCellPredicate
{
public:
  /** */
  CoordinateValueCellPredicate(int direction,
    const double& value, const double& tol=1.0e-12)
    : PositionalCellPredicate(
      rcp(new CoordinateValueCellPredicateFunctor(direction,value,tol)))
    {}

  /* */
  GET_RCP(CellPredicateBase);
};


}


#endif
