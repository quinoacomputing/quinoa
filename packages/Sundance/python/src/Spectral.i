// -*- c++ -*-


%{

#define HAVE_PY_FIAT
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceSpectralBasis.hpp"
#include "SundanceHermiteSpectralBasis.hpp"
  %}


// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"


namespace Sundance
{
  class SpectralBasis
  {
  public:
    SpectralBasis();
    ~SpectralBasis();

    int getDim() const ;
    
    int getOrder() const ;

    int nterms() const ;

    int getElement(int i) const ;

    double expectation(int i, int j, int k) const ;

    std::string toString() const ;
  };

  %extend SpectralBasis
  {
    using namespace std;
    std::string __str__() 
    {
      return self->toString();
    }
  }


}



%rename(HermiteSpectralBasis) makeHermiteSpectralBasis;
%rename(SpectralExpr) makeSpectralExpr;

%inline %{
  /* Create a Hermite basis */
  Sundance::SpectralBasis makeHermiteSpectralBasis(int dim, int order)
  {
    return new Sundance::HermiteSpectralBasis(dim, order);
  }

  /* Create a Hermite basis */
  Sundance::SpectralBasis makeHermiteSpectralBasis(int dim, int order, int nterms)
  {
    return new Sundance::HermiteSpectralBasis(dim, order, nterms);
  }

  /* Create a spectral expression */
  Sundance::Expr makeSpectralExpr(const Sundance::SpectralBasis& sb, 
                                      const Sundance::Expr& coeffs)
  {
    return new Sundance::SpectralExpr(sb, coeffs);
  }


  %}




