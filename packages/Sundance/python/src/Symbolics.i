// -*- c++ -*-

%{
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceExpr.hpp"
#include "SundanceUnknownFunction.hpp"
#include "SundanceSpectralBasis.hpp"
#include "SundanceTestFunction.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceIntegral.hpp"
#include "SundanceParameter.hpp"
#include "SundanceFunctionalEvaluator.hpp"
#include "SundanceStdMathOps.hpp"
#include "PlayaVectorImpl.hpp"
  %}


// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"

%feature("autodoc");

namespace Sundance
{
  class Expr
  {
  public:
    Expr();
    ~Expr();
    Expr(const double& c);

    void setParameterValue(const double& val);

    /** Number of elements in top level of list */
    int size() const ;

    /** Total number of elements in list. */
    int totalSize() const ;
    
    /** Append a new element to this list */
    void append(const Expr& expr);
    
    /** Flatten this list */
    Expr flatten() const ;

    /** Return real part of a complex expression */
    Expr real() const ;
    
    /** Return imaginary part of a complex expression */
    Expr imag() const ;
    
    /** Return complex conjugate */
    Expr conj() const ;

  };

  %extend Expr
  {
    using namespace std;
    std::string __str__() 
    {
      std::string rtn = self->toString(); 
      return rtn;
    }

    std::string fullForm() const
    {
      return self->toXML().toString();
    }

    /* Return the Spectral Basis */
    SpectralBasis getSpectralBasis() const
    {
      return Sundance::getSpectralBasis(*self);
    }

    /* Return the coefficient of the nth spectral basis term */
    Expr getSpectralCoeff(int n) const
    {
      return Sundance::getSpectralCoeff(n, *self);
    }

    double integral(const Sundance::CellFilter& domain,
                    const Sundance::Mesh& mesh,
                    const Sundance::QuadratureFamily& quad)
    {
      Expr I = Integral(domain, *self, quad);
      return Sundance::evaluateIntegral(mesh, I);
    }

    Expr __pow__(const double& ex)
    {
      return pow(*self, ex);
    }

    Expr __add__(const Sundance::Expr& other) 
    {
      Sundance::Expr rtn = self->operator+(other);
      return rtn;
    }

    Expr __sub__(const Sundance::Expr& other) 
    {
      Sundance::Expr rtn = self->operator-(other);
      return rtn;
    }

    Expr __mul__(const Sundance::Expr& other) 
    {
      Sundance::Expr rtn = self->operator*(other);
      return rtn;
    }

    Expr __div__(const Sundance::Expr& other) 
    {
      Sundance::Expr rtn = self->operator/(other);
      return rtn;
    }

    

    /* operations with scalars to the right */

    Expr __add__(const double& other) 
    {
      Sundance::Expr rtn = self->operator+(other);
      return rtn;
    }

    Expr __sub__(const double& other) 
    {
      Sundance::Expr rtn = self->operator-(other);
      return rtn;
    }

    Expr __mul__(const double& other) 
    {
      Sundance::Expr rtn = self->operator*(other);
      return rtn;
    }

    Expr __div__(const double& other) 
    {
      Sundance::Expr rtn = self->operator/(other);
      return rtn;
    }

    Expr __add__(const complex<double>& other) 
    {
      Sundance::Expr rtn = self->operator+(other);
      return rtn;
    }

    Expr __sub__(const complex<double>& other) 
    {
      Sundance::Expr rtn = self->operator-(other);
      return rtn;
    }

    Expr __mul__(const complex<double>& other) 
    {
      Sundance::Expr rtn = self->operator*(other);
      return rtn;
    }

    Expr __div__(const complex<double>& other) 
    {
      Sundance::Expr rtn = self->operator/(other);
      return rtn;
    }

    /* operations with scalars to the left */

    Expr __radd__(const double& other) 
    {
      Sundance::Expr rtn = other + *self;
      return rtn;
    }

    Expr __rsub__(const double& other) 
    {
      Sundance::Expr rtn = other - *self;
      return rtn;
    }

    Expr __rmul__(const double& other) 
    {
      Sundance::Expr rtn = other * (*self);
      return rtn;
    }

    Expr __rdiv__(const double& other) 
    {
      Sundance::Expr rtn = other / (*self);
      return rtn;
    }

    Expr __radd__(const complex<double>& other) 
    {
      Sundance::Expr rtn = other + *self;
      return rtn;
    }

    Expr __rsub__(const complex<double>& other) 
    {
      Sundance::Expr rtn = other - *self;
      return rtn;
    }

    Expr __rmul__(const complex<double>& other) 
    {
      Sundance::Expr rtn = other * (*self);
      return rtn;
    }

    Expr __rdiv__(const complex<double>& other) 
    {
      Sundance::Expr rtn = other / (*self);
      return rtn;
    }

    Expr __div__(const double& other) 
    {
      Sundance::Expr rtn = (*self)/other;
      return rtn;
    }

    /* unary operations */
    Expr __neg__() 
    {
      Sundance::Expr rtn = self->operator-();
      return rtn;
    }


    /* list indexing and information */
    
    Expr __getitem__(int i) const
    {
      return self->operator[](i);
    }

    /* get the vector underlying a discrete function */
    Playa::Vector<double> getVector() const 
    {
      /* cast to a discrete function. The validity of the cast
       * is checked within discFunc(). */
      const DiscreteFunction* df = DiscreteFunction::discFunc(*self);
      return df->getVector();
    }

    /* get the vector underlying a discrete function */
    void setVector(const Playa::Vector<double>& vec) 
    {
      /* cast to a discrete function. The validity of the cast
       * is checked within discFunc(). */
      DiscreteFunction* df = DiscreteFunction::discFunc(*self);
      df->setVector(vec);
    }

    /* get the discrete space associated with an expression */
    Sundance::DiscreteSpace discSpace() const
    {
      /* cast to a discrete function. The validity of the cast
       * is checked within discFunc(). */
      const DiscreteFunction* df = DiscreteFunction::discFunc(*self);
      return df->discreteSpace();
    }
    
    

    


    
  }

  Expr Complex(const Expr& real, const Expr& imag);

  Expr conj(const Expr& x);

  Expr Re(const Expr& x);
  Expr Im(const Expr& x);

  Expr List(const Expr& a);
  Expr List(const Expr& a, const Expr& b);
  Expr List(const Expr& a, const Expr& b, const Expr& c);
  Expr List(const Expr& a, const Expr& b, const Expr& c, const Expr& d);
  Expr List(const Expr& a, const Expr& b, const Expr& c, const Expr& d,
            const Expr& e);
  Expr List(const Expr& a, const Expr& b, const Expr& c,
            const Expr& d, const Expr& e, const Expr& f);

  Expr pow(const Expr& x, const double& a);

  Expr fabs(const Expr& x);
  Expr sign(const Expr& x);

  Sundance::Expr Sundance::exp(const Sundance::Expr& x);
  Expr log(const Expr& x);
  Expr sqrt(const Expr& x);

  Expr sin(const Expr& x);
  Expr cos(const Expr& x);
  Expr tan(const Expr& x);
  Expr asin(const Expr& x);
  Expr acos(const Expr& x);
  Expr atan(const Expr& x);

  Expr sinh(const Expr& x);
  Expr cosh(const Expr& x);
  Expr tanh(const Expr& x);
  Expr asinh(const Expr& x);
  Expr acosh(const Expr& x);
  Expr atanh(const Expr& x);

/** \relates Expr */
Expr gradient(int dim);
  
/** \relates Expr */
Expr div(const Expr& f);
  
/** \relates Expr */
Expr cross(const Expr& a, const Expr& b);
  
/** \relates Expr */
Expr curl(const Expr& f);

/** \relates Expr \relates CellVectorExpr */
Expr CellNormalExpr(int dimension, const std::string& name);

/** \relates Expr \relates CellVectorExpr */
Expr CellTangentExpr(int dimension, const std::string& name);

}



%rename(UnknownFunction) makeUnknownFunction;
%rename(TestFunction) makeTestFunction;
%rename(CoordExpr) makeCoordExpr;
%rename(Derivative) makeDerivative;
%rename(Parameter) makeParameter;
%rename(CellDiameterExpr) makeCellDiameterExpr;


%inline %{
  /* Create an unknown function */
  Sundance::Expr makeUnknownFunction(const Sundance::BasisFamily& b,
                                         const std::string& name)
  {
    return new Sundance::UnknownFunction(b, name);
  }
  %}
%inline %{
  /* Create an unknown function */
  Sundance::Expr makeUnknownFunction(const Sundance::BasisFamily& b,
                                         const Sundance::SpectralBasis& sb,
                                         const std::string& name)
  {
    return new Sundance::UnknownFunction(b, sb, name);
  }
  %}
%inline %{
  /* Create an unknown function */
  Sundance::Expr makeUnknownFunction(const Sundance::BasisFamily& b,
                                         const Sundance::SpectralBasis& sb)
  {
    return new Sundance::UnknownFunction(b, sb);
  }
  %}

%inline %{
  /* Create an unknown function */
  Sundance::Expr makeUnknownFunction(const Sundance::BasisFamily& b)
  {
    return new Sundance::UnknownFunction(b);
  }
  %}


%inline %{
  /* Create a test function */
  Sundance::Expr makeTestFunction(const Sundance::BasisFamily& b,
                                      const std::string& name)
  {
    return new Sundance::TestFunction(b, name);
  }
  %}


%inline %{
  /* Create a test function */
  Sundance::Expr makeTestFunction(const Sundance::BasisFamily& b)
  {
    return new Sundance::TestFunction(b);
  }
  %}

%inline %{
  /* Create an unknown function */
  Sundance::Expr makeTestFunction(const Sundance::BasisFamily& b,
                                      const Sundance::SpectralBasis& sb,
                                      const std::string& name)
  {
    return new Sundance::TestFunction(b, sb, name);
  }
  %}
%inline %{
  /* Create an unknown function */
  Sundance::Expr makeTestFunction(const Sundance::BasisFamily& b,
                                      const Sundance::SpectralBasis& sb)
  {
    return new Sundance::TestFunction(b, sb);
  }
  %}

%inline %{
  /* Create a coordinate expression */
  Sundance::Expr makeCoordExpr(int dir)
  {
    return new Sundance::CoordExpr(dir);
  }
  %}

%inline %{
  /* Create a cell diameter expression */
  Sundance::Expr makeCellDiameterExpr()
  {
    return new Sundance::CellDiameterExpr();
  }
  %}


%inline %{
  /* Create a differential operator */
  Sundance::Expr makeDerivative(int dir)
  {
    return new Sundance::Derivative(dir);
  }
  %}


%inline %{
  /* Create a differential operator */
  Sundance::Expr makeParameter(const double& val)
  {
    return new Sundance::Parameter(val);
  }
  %}


%inline %{
  /* Create a differential operator */
  Sundance::Expr makeParameter(const double& val, const std::string& name)
  {
    return new Sundance::Parameter(val, name);
  }
  %}
