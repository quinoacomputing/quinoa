// -*- c++ -*-


%{


  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceBasisFamily.hpp"
#include "SundanceLagrange.hpp"
#ifdef HAVE_FIAT
#include "SundanceFIATLagrange.hpp"
#include "PySundanceFIATScalarAdapter.hpp"
#include "PySundanceBasisCheck.hpp"
#endif
  %}


// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"


namespace Sundance
{
  class BasisFamily
  {
  public:
    BasisFamily();
    ~BasisFamily();
  };

  %extend BasisFamily
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
  }


  class BasisArray
  {
  public:
    BasisArray();
    BasisArray(int n);

    void append(const BasisFamily& b);

  };

  %extend BasisArray
  {
    using namespace std;
    std::string __str__() 
    {
      std::string rtn; 
      std::stringstream os;
      os << *self;
      rtn = os.str();
      return rtn;
    }
  }
    

}



%rename(Lagrange) makeLagrange;
#ifdef HAVE_FIAT
%rename(FIATLagrange) makeFIATLagrange;
%rename(FIATScalarAdapter) makeFIATScalarAdapter;
#endif

%inline %{
  /* Create a Lagrange basis function */
  Sundance::BasisFamily makeLagrange(int order)
  {
    return new Sundance::Lagrange(order);
  }
#ifdef HAVE_FIAT
  Sundance::BasisFamily makeFIATLagrange(int order)
  {

    return new Sundance::FIATLagrange(order);
  }

  Sundance::BasisFamily makeFIATScalarAdapter(PyObject *py_basis ,
						    int order)
  {
    return new Sundance::FIATScalarAdapter(py_basis,order);
  }

#endif
  /* */
  Sundance::BasisArray 
    BasisList()
  {
    return BasisArray();
  }

  /* */
  Sundance::BasisArray 
    BasisList(const Sundance::BasisFamily& a)
  {
    return Array<BasisFamily>(tuple(a));
  }

  /* */
  Sundance::BasisArray 
    BasisList(const Sundance::BasisFamily& a,
              const Sundance::BasisFamily& b)
  {
    return Array<BasisFamily>(tuple(a,b));
  }

  /* */
  Sundance::BasisArray 
    BasisList(const Sundance::BasisFamily& a,
              const Sundance::BasisFamily& b,
              const Sundance::BasisFamily& c)
  {
    return  Array<BasisFamily>(tuple(a,b,c));
  }

  /* */
  Sundance::BasisArray 
    BasisList(const Sundance::BasisFamily& a,
              const Sundance::BasisFamily& b,
              const Sundance::BasisFamily& c,
              const Sundance::BasisFamily& d)
  {
    return  Array<BasisFamily>(tuple(a,b,c,d));
  }

  /* */
  Sundance::BasisArray 
    BasisList(const Sundance::BasisFamily& a,
              const Sundance::BasisFamily& b,
              const Sundance::BasisFamily& c,
              const Sundance::BasisFamily& d,
              const Sundance::BasisFamily& e)
  {
    return  Array<BasisFamily>(tuple(a,b,c,d,e));
  }
                                      
  %}




