// -*- c++ -*-


%{

#define HAVE_PY_FIAT
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceStochBlockJacobiSolver.hpp"

  %}


// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"

namespace Sundance
{
  class StochBlockJacobiSolver
  {
  public:
    StochBlockJacobiSolver(
    const Playa::LinearSolver<double>& diagonalSolver,
    const SpectralBasis& pcBasis, 
    double convTol,
    int maxIters,
    int verbosity);

    void solve(const Teuchos::Array<Playa::LinearOperator<double> >& KBlock,
      const Teuchos::Array<int>& hasNonzeroMatrixBlock,
      const Teuchos::Array<Playa::Vector<double> >& fBlock,
      Teuchos::Array<Playa::Vector<double> >& xBlock) const ;

    
    void solve(const Teuchos::Array<Playa::LinearOperator<double> >& KBlock,
      const Teuchos::Array<Playa::Vector<double> >& fBlock,
      Teuchos::Array<Playa::Vector<double> >& xBlock) const ;
  };

}





%template(IntArray) Teuchos::Array<int>;
%template(LinearOperatorArray) Teuchos::Array<Playa::LinearOperator<double> >;



