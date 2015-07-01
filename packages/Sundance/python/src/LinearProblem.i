// -*- c++ -*-


%{
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceLinearProblem.hpp"
#include "PySundanceLinearSolver.hpp"

  %}

namespace Sundance
{
class Block
  {
  public:
    /** */
    Block(const Sundance::Expr& expr, const Playa::VectorType<double>& vecType);

    /** */
    const Sundance::Expr& expr() const ;

    /** */
    const Playa::VectorType<double>& vecType() const ;
  };

  class BlockArray
  {
  public:
    BlockArray(int n);
  };

  %extend Block
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
    
  %extend BlockArray
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


%inline %{
/* */
  Sundance::BlockArray 
    BlockList(const Sundance::Block& a)
  {
    return Array<Block>(tuple(a));
  }

  /* */
  Sundance::BlockArray 
    BlockList(const Sundance::Block& a,
              const Sundance::Block& b)
  {
    return Array<Block>(tuple(a,b));
  }

  /* */
  Sundance::BlockArray 
    BlockList(const Sundance::Block& a,
              const Sundance::Block& b,
              const Sundance::Block& c)
  {
    return Array<Block>(tuple(a,b,c));
  }

  /* */
  Sundance::BlockArray 
    BlockList(const Sundance::Block& a,
              const Sundance::Block& b,
              const Sundance::Block& c,
              const Sundance::Block& d)
  {
    return Array<Block>( tuple(a,b,c,d) );
  }

  /* */
  Sundance::BlockArray 
    BlockList(const Sundance::Block& a,
              const Sundance::Block& b,
              const Sundance::Block& c,
              const Sundance::Block& d,
              const Sundance::Block& e)
  {
    return Array<Block>(tuple(a,b,c,d,e));
  }
                                      
  %}

// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"

namespace Sundance
{
  
  typedef RCP<DOFMapBase> DOFMap;

  class LinearProblem
  {
  public:
    LinearProblem(const Sundance::Mesh& mesh, 
                  const Sundance::Expr& eqn,
                  const Sundance::Expr& bc,
                  const Sundance::Expr& v, 
                  const Sundance::Expr& u,
                  const Playa::VectorType<double>& vecType);

    LinearProblem(const Sundance::Mesh& mesh, 
                  const Sundance::Expr& eqn,
                  const Sundance::Expr& bc,
                  const Sundance::BlockArray& v, 
                  const Sundance::BlockArray& u);

    Playa::Vector<double> getSingleRHS() const ;

    Playa::LinearOperator<double> getOperator() const ;

    Sundance::Expr solve(const Playa::LinearSolver<double>& solver) const ;

  };

  %extend LinearProblem {
    Sundance::Expr solve(PyObject* pySolver)
    {
      Playa::PySundanceLinearSolver* tmp  
        = new Playa::PySundanceLinearSolver(pySolver);
      
      Teuchos::RCP<Playa::LinearSolverBase<double> > 
        r = rcp(tmp);

      Playa::LinearSolver<double> cxxSolver = r;
      return self->solve(cxxSolver);
    }
  }

  %extend LinearProblem {
    void printRowMap() const 
    {
      for (int b=0; b<self->numBlockRows(); b++)
        {
          self->rowMap(b)->print(cout);
        }
    }
  }

  %extend LinearProblem {
    void printColMap() const 
    {
      for (int b=0; b<self->numBlockCols(); b++)
        {
          self->colMap(b)->print(cout);
        }
    }
  }
 
}
