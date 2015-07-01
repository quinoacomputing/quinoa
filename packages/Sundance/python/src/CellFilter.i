// -*- c++ -*-


%{
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceCellFilter.hpp"
#include "SundanceCellSet.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceBoundaryCellFilter.hpp"
#include "SundanceDimensionalCellFilter.hpp"
#include "PySundanceCellPredicate.hpp"

  %}


// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"

namespace Sundance
{

class CellPredicate
{
public:
  CellPredicate();
  ~CellPredicate();


};

%extend CellPredicate
{
  std::string __str__() 
  {
    std::string rtn; 
    std::stringstream os;
    self->print(os);
    rtn = os.str();
    return rtn;
  }
}

class CellSet
{
public:
  CellSet();
  ~CellSet();

  %extend 
    {
      std::string __str__() 
      {
        std::string rtn; 
        std::stringstream os;
        self->print(os);
        rtn = os.str();
        return rtn;
      }
    }
};


class CellFilter
{
public:
  CellFilter();
  ~CellFilter();

    
  CellSet getCells(const Sundance::Mesh& mesh) const ;
  int dimension(const Sundance::Mesh& mesh) const ;

  CellFilter labeledSubset(int label) const ;

  CellFilter intersection(const CellFilter& other) const ;

  CellFilter subset(const CellPredicate& cp) const ;
};

%extend CellFilter
{
  using namespace std;
  std::string __str__() 
  {
    std::string rtn; 
    std::stringstream os;
    self->print(os);
    rtn = os.str();
    return rtn;
  }

  CellFilter subset(PyObject* functor) const
  {
    Teuchos::RCP<Sundance::CellPredicateFunctorBase> f 
      = Teuchos::rcp(new Sundance::PySundanceCellPredicate(functor));
    CellPredicate p = new Sundance::PositionalCellPredicate(f);
    return self->subset(p);
  }

  CellFilter __add__(const CellFilter& other) const
  {
    return self->operator+(other);
  }

  CellFilter __sub__(const CellFilter& other) const
  {
    return self->operator-(other);
  }
}


class CellFilterArray
{
public:
  CellFilterArray();
  ~CellFilterArray();

  int size() const ;

  void append(const CellFilter& x);
};

}

%rename(MaximalCellFilter) makeMaximalCellFilter;
%rename(BoundaryCellFilter) makeBoundaryCellFilter;
%rename(DimensionalCellFilter) makeDimensionalCellFilter;
%rename(PositionalCellPredicate) makePyFunctorCellPredicate;
%rename(CoordinateValueCellPredicate) makeCoordinateValueCellPredicate;


%inline %{
  /* Create a maximal cell filter */
  Sundance::CellFilter makeMaximalCellFilter()
  {
    return new Sundance::MaximalCellFilter();
  }
  %}


%inline %{
  /* Create a boundary cell filter */
  Sundance::CellFilter makeBoundaryCellFilter()
  {
    return new Sundance::BoundaryCellFilter();
  }
  %}

%inline %{
  /* Create a dimensional cell ftiler */
  Sundance::CellFilter makeDimensionalCellFilter(int i)
  {
    return new Sundance::DimensionalCellFilter(i);
  }
  %}

%inline %{
  /*  */
  Sundance::CellPredicate makePyFunctorCellPredicate(PyObject* functor)
  {
    Teuchos::RCP<Sundance::CellPredicateFunctorBase> f 
      = Teuchos::rcp(new Sundance::PySundanceCellPredicate(functor));
    return new Sundance::PositionalCellPredicate(f);
  }
  %}

%inline %{
  /*  */
  Sundance::CellPredicate 
    makeCoordinateValueCellPredicate(int dir, const double& val)
  {
    return new Sundance::CoordinateValueCellPredicate(dir, val);
  }
  %}



%inline %{
  /*  */
  Sundance::CellFilterArray CellFilterList()
  {
    return CellFilterArray();
  }
  /*  */
  Sundance::CellFilterArray CellFilterList(const Sundance::CellFilter& a)
  {
    return Sundance::List(a);
  }

  /*  */
  Sundance::CellFilterArray CellFilterList(const Sundance::CellFilter& a,
    const Sundance::CellFilter& b)
  {
    return Sundance::List(a, b);
  }

  /*  */
  Sundance::CellFilterArray CellFilterList(const Sundance::CellFilter& a,
    const Sundance::CellFilter& b,
    const Sundance::CellFilter& c)
  {
    return Sundance::List(a, b, c);
  }

  /*  */
  Sundance::CellFilterArray CellFilterList(const Sundance::CellFilter& a,
    const Sundance::CellFilter& b,
    const Sundance::CellFilter& c,
    const Sundance::CellFilter& d)
  {
    return Sundance::List(a, b, c, d);
  }

  /*  */
  Sundance::CellFilterArray CellFilterList(const Sundance::CellFilter& a,
    const Sundance::CellFilter& b,
    const Sundance::CellFilter& c,
    const Sundance::CellFilter& d,
    const Sundance::CellFilter& e)
  {
    return Sundance::List(a, b, c, d, e);
  }

  %}


