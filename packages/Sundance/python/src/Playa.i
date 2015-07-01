// -*- c++ -*-

%{
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceDefs.hpp"
#include "PlayaVectorType.hpp"
#include "PlayaVectorSpaceDecl.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaEpetraVectorType.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaSimpleAddedOpDecl.hpp"
#include "PlayaSimpleComposedOpDecl.hpp"
#include "PlayaSimpleBlockOpDecl.hpp"
#include "PlayaSimpleScaledOpDecl.hpp"
#include "PlayaSimpleZeroOpDecl.hpp"
#include "PlayaSimpleIdentityOpDecl.hpp"
#include "PlayaLinearSolverDecl.hpp"
#include "PlayaPreconditionerFactory.hpp"
#include "PlayaPreconditioner.hpp"
#include "PlayaGenericLeftPreconditioner.hpp"
#include "PlayaGenericRightPreconditioner.hpp"
#include "PlayaBlockVectorSpaceDecl.hpp"
#include "PlayaDefaultBlockVectorSpaceDecl.hpp"
#include "PlayaBICGSTABSolverDecl.hpp"
#include "PlayaNOXSolver.hpp"
#include "PlayaLinearSolverBuilder.hpp"
#include "PlayaLinearCombinationDecl.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "PyTeuchos_Utils.hpp"
//#include "PySundanceNOXSolverHandle.hpp"
#include "PySundanceLinearSolver.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorSpaceImpl.hpp"
#include "PlayaVectorImpl.hpp"
#include "PlayaBICGSTABSolverImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaLinearSolverImpl.hpp"
#include "PlayaDefaultBlockVectorSpaceImpl.hpp"
#include "PlayaSimpleBlockOpImpl.hpp"
#include "PlayaSimpleZeroOpImpl.hpp"
#include "PlayaSimpleIdentityOpImpl.hpp"
#include "PlayaSimpleAddedOpImpl.hpp"
#include "PlayaSimpleComposedOpImpl.hpp"
#include "PlayaSimpleScaledOpImpl.hpp"
#endif

#include "PlayaHack.hpp"

  %}




// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"


%template(doubleVector) std::vector<double>;

%rename(Vector) Vec;
%rename(VectorType) VecType;
%rename(VectorSpace) VecSpace;
%rename(LinearOperator) LinOp;
%rename(LinearSolver) LinSol;
%rename(NonlinearOperator) NonlinOp;
//%rename(NOXSolver) makeNOXSolver;
%rename(Preconditioner) Precond;
%rename(PreconditionerFactory) PrecondFactory;



/* --------- vector space ------------ */
namespace Playa
{
  template <class Scalar> class Vector;
  template <class Scalar>
  class VectorSpace
  {
  public:
    Vector<Scalar> createMember();

    int dim() const ;
  };

  %template(VecSpace) VectorSpace<double>;
}

/* --------- vector ------------ */
namespace Playa
{
  template <class Scalar> class Vector
  {
  public:
    Vector();
    ~Vector();

    VectorSpace<Scalar> space() const ;

    Vector<Scalar> copy() const ;

    Vector<Scalar> acceptCopyOf(const Vector<Scalar>& x);

    Vector<Scalar> dotStar(const Vector<Scalar>& other) const ;

    Vector<Scalar> dotSlash(const Vector<Scalar>& other) const ;

    Vector<Scalar> reciprocal() const ;

    Vector<Scalar> abs() const ;

    void setToConstant(const Scalar& alpha) ;

    Scalar norm1() const ;

    Scalar norm2() const ;

    Scalar normInf() const ;

    void zero();

    Scalar max() const;

    Scalar min()const;

    Vector<Scalar> getBlock(int i) const  ;

    void setBlock(int i, const Vector<Scalar>& x) ;


    %extend 
    {
      int numBlocks() const
      {
        return self->space().numBlocks();
      }

      std::string __str__() 
      {
        std::string rtn; 
        std::stringstream os;
        self->print(os);
        rtn = os.str();
        return rtn;
      }

      Vector<Scalar> __add__(const Vector<Scalar>& other) 
      {
        return (*self) + other;
      }

      Vector<Scalar> __sub__(const Vector<Scalar>& other) 
      {
        return (*self) - other;
      }

      Vector<Scalar> __mul__(const Scalar& other) 
      {
        return (*self) * other;
      }

      Vector<Scalar> __div__(const Scalar& other) 
      {
        return (*self) *(1.0/ other);
      }

      Vector<Scalar> __rmul__(const Scalar& other) 
      {
        return (*self) * other;
      }

      Scalar __mul__(const Vector<Scalar>& other) 
      {
        return (*self) * other;
      }

      
      Scalar __getitem__(int localIndex) const 
      {
        return self->operator[](localIndex);
      }
      
      void __setitem__(int localIndex, const Scalar& value)
      {
        self->operator[](localIndex) = value;
      }
    }
  };

  %template(Vec) Vector<double>;

}

/* --------- vector space ------------ */
namespace Playa
{
  template <class Scalar>
  class VectorSpace
  {
  public:
    VectorSpace();
    ~VectorSpace();

    Vector<Scalar> createMember();


    /** return the number of subblocks. */
    int numBlocks() const ;

    /** get the i-th subblock */
    VectorSpace<Scalar> getBlock(int i) const ;


    /** set the i-th subblock */
    void setBlock(int i, const VectorSpace<Scalar>& space);

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

  %template(VecSpace) VectorSpace<double>;
}


%rename(BlockVectorSpace) makeBlockVectorSpace;

%inline %{
  /* Create block vector space */
  Playa::VectorSpace<double> 
    makeBlockVectorSpace(Playa::VectorSpace<double> vs)
  {
    return Playa::blockSpace(vs);
  }

  /* Create block vector space */
  Playa::VectorSpace<double> 
    makeBlockVectorSpace(Playa::VectorSpace<double> vs1,
                         Playa::VectorSpace<double> vs2)
  {
    return Playa::blockSpace(vs1, vs2);
  }

  /* Create block vector space */
  Playa::VectorSpace<double> 
    makeBlockVectorSpace(Playa::VectorSpace<double> vs1,
                         Playa::VectorSpace<double> vs2,
                         Playa::VectorSpace<double> vs3)
  {
    return Playa::blockSpace(vs1, vs2, vs3);
  }

  %}


/* --------- vector type ------------ */
namespace Playa
{
  template <class Scalar>
  class VectorType
  {
  public:
    ~VectorType();
    VectorType();

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

    %extend
    {
      VectorSpace<Scalar> createEvenlyPartitionedSpace(int nLocal) const 
      {
        return self->createEvenlyPartitionedSpace(MPIComm::world(), nLocal);
      }
    }
  };

  %template(VecType) VectorType<double>;

}



/* --------- vector type ------------ */
namespace Playa
{
  enum SolverStatusCode {SolveCrashed, SolveFailedToConverge, SolveConverged};
  
  template <class Scalar>
  class SolverState
  {
  public:
    SolverState();
    SolverState(SolverStatusCode finalState, const std::string& msg, 
                int finalIters, const Scalar& finalResid);
    ~SolverState();
    
    std::string stateDescription() const ;

    /** */
    const Scalar& finalResid() const ;

    /** */
    int finalIters() const ;

    /** */
    const SolverStatusCode& finalState() const ;

    /** */
    const std::string& finalMsg() const ;

    %extend
    {
      std::string __str__() 
      {
        std::string rtn; 
        std::stringstream os;
        os << *self;
        rtn = os.str();
        return rtn;
      }
    }
  };

}



%rename(EpetraVectorType) makeEpetraVectorType;

%inline %{
  /* Create an epetra vector type */
  Playa::VectorType<double> makeEpetraVectorType()
  {
    return Playa::VectorType<double>(new Playa::EpetraVectorType());
  }
  %}


/* --------- linear operator ------------ */
namespace Playa
{
  template <class Scalar>
  class LinearOperator
  {
  public:
    LinearOperator();
    ~LinearOperator();

    /** Return the domain */
    const VectorSpace<Scalar> domain() const ;

    /** Return the range */
    const VectorSpace<Scalar> range() const ;


    /** return number of block rows */
    int numBlockRows() const;
      

    /** return number of block cols */
    int numBlockCols() const;
      

    /** get the (i,j)-th block */
    LinearOperator<Scalar> getBlock(const int &i, const int &j) const ;

    /** set the (i,j)-th block 
     *  If the domain and/or the range are not set, then we
     *  are building the operator
     */
    void setBlock(int i, int j, 
                  const LinearOperator<Scalar>& sub);

    /**
     * Return a TransposeOperator.
     */
    LinearOperator<Scalar> transpose() const ; 


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

      Vector<Scalar> __mul__(const Vector<Scalar>& in) const 
      {
        Vector<Scalar> out;
        self->apply(in, out);
        return out;
      }

      LinearOperator<Scalar> __mul__(const LinearOperator<Scalar>& other) const 
      {
        return (*self) * other;
      }

      LinearOperator<Scalar> __add__(const LinearOperator<Scalar>& other) const 
      {
        return (*self) + other;
      }

      LinearOperator<Scalar> __sub__(const LinearOperator<Scalar>& other) const 
      {
        return (*self) + -1.0*other;
      }
      
      LinearOperator<Scalar> __mul__(const Scalar& other) const 
      {
        return other * (*self);
      }
      
      LinearOperator<Scalar> __rmul__(const Scalar& other) const 
      {
        return other * (*self);
      }

      
    }
  };

  %template(LinOp) LinearOperator<double>;

}

%inline %{
  /* Create block operator */
  Playa::LinearOperator<double> 
    BlockOperator(const Playa::VectorSpace<double>& domain,
      const Playa::VectorSpace<double>& range)
  {
    return makeBlockOperator(domain, range);
  }

  %}



%rename(IdentityOperator) makeIdentityOperator;

%inline %{
  /* Create block operator */
  Playa::LinearOperator<double> 
    makeIdentityOperator(const Playa::VectorSpace<double>& space)
  {
    return identityOperator(space);
  }

  %}

%rename(ZeroOperator) makeZeroOperator;

%inline %{
  /* Create block operator */
  Playa::LinearOperator<double> 
    makeZeroOperator(const Playa::VectorSpace<double>& domain,
      const Playa::VectorSpace<double>& range)
  {
    return zeroOperator(domain, range);
  }

  %}


/* --------- linear solver ------------ */
namespace Playa
{
  template <class Scalar>
  class LinearSolver
  {
  public:
    LinearSolver();
    ~LinearSolver();

    SolverState<Scalar> solve(const LinearOperator<Scalar>& op,
                              const Vector<Scalar>& rhs,
                              Vector<Scalar>& soln) const ;

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

  %template(LinSol) LinearSolver<double>;


 //  %extend LinOp 
//   {
//     LinOp inverse(const LinearSolver<Scalar>& solver) const
//     {
//       return LinOp(new InverseOperator<double>(*self, solver));
//     }
//   }

}



%rename(BICGSTABSolver) makeBICGSTABSolver;



%inline %{
  Playa::LinearSolver<double> makeBICGSTABSolver(const Teuchos::ParameterList& params)
  {
    return LinearSolver<double>(new Playa::BICGSTABSolver<double>(params));
  }
  %}


%inline %{
  Playa::LinearSolver<double> makeBICGSTABSolver(const Teuchos::ParameterList& params,
                                                    const Playa::PreconditionerFactory<double>& precond)
  {
    return LinearSolver<double>(new Playa::BICGSTABSolver<double>(params, precond));
  }
  %}


 

%inline %{
  /* Read a linear solver from an XML file */
  Playa::LinearSolver<double> readSolver(const std::string& filename)
  {
    Teuchos::ParameterXMLFileReader reader(filename);
    Teuchos::ParameterList solverParams = reader.getParameters();
    Playa::LinearSolver<double> solver 
      = Playa::LinearSolverBuilder::createSolver(solverParams);
    return solver;
  }
  %}

%inline %{
  /* Read a linear solver from a parameter list */
  Playa::LinearSolver<double> buildSolver(const Teuchos::ParameterList& params)
  {
    Playa::LinearSolver<double> solver ;
    try
      {
        solver = Playa::LinearSolverBuilder::createSolver(params);
      }
    catch(std::exception& e)
      {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        std::cerr << "detected exception "
                  << e.what() << " in buildSolver()" << std::endl;
      }
    return solver;
  }
  %}


/* --------- preconditioner ------------ */
namespace Playa
{
  template <class Scalar>
  class Preconditioner
  {
  public:
    Preconditioner();
    ~Preconditioner();

    
    
    /** Left preconditioner */
    LinearOperator<Scalar> left() const ;
    
    /** Right preconditioner */
    LinearOperator<Scalar> right() const ;
    
    /** return true if this preconditioner has both left and
     * right components. */
    bool isTwoSided() const ;
    
    /** return true if this preconditioner has a nontrivial left component */
    bool hasLeft() const ;
    
    /** return true if this preconditioner has
     * a nontrivial right component */
    bool hasRight() const ;
    
    /** return true if this preconditioner has neither left nor
     * right operators defined */
    bool isIdentity() const ;

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

  %template(Precond) Preconditioner<double>;
}


/* --------- preconditioner factory ------------ */
namespace Playa
{
  template <class Scalar>
  class PreconditionerFactory
  {
  public:
    PreconditionerFactory();
    ~PreconditionerFactory();

    Preconditioner<Scalar> createPreconditioner(const LinearOperator<Scalar>& A) const ;

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

  %template(PrecondFactory) PreconditionerFactory<double>;
}


%rename(ILUKPreconditionerFactory) makeILUKPreconditionerFactory;
%inline %{
  /* Create ILUK preconditioner factory  */
  Playa::PreconditionerFactory<double> 
    makeILUKPreconditionerFactory(const Teuchos::ParameterList& params)
  {
    return Playa::PreconditionerFactory<double>(new Playa::ILUKPreconditionerFactory<double>(params));
  }
  %}

%rename(GenericLeftPreconditioner) makeGenericLeftPreconditioner;

%inline %{
  /* Create generic left preconditioner  */
  Playa::Preconditioner<double> 
    makeGenericLeftPreconditioner(const Playa::LinearOperator<double>& left)
  {
    return Playa::Preconditioner<double>(new Playa::GenericLeftPreconditioner<double>(left));
  }
  %}

%rename(GenericRightPreconditioner) makeGenericRightPreconditioner;

%inline %{
  /* Create generic left preconditioner  */
  Playa::Preconditioner<double> 
    makeGenericRightPreconditioner(const Playa::LinearOperator<double>& op)
  {
    return Playa::Preconditioner<double>(new Playa::GenericRightPreconditioner<double>(op));
  }
  %}





/* --------- nonlinear operator ------------ */
namespace Playa
{
  template <class Scalar>
  class NonlinearOperator
  {
  public:
    NonlinearOperator();
    ~NonlinearOperator();

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

  %template(NonlinOp) NonlinearOperator<double>;

}

namespace NOX
{
  namespace StatusTest
  {
    enum StatusType {Unevaluated, Unconverged, Converged, Failed};
  }
}


namespace Playa
{
class NOXSolver 
{
public:
  /** */
  NOXSolver();
  /** */
  NOXSolver(const Teuchos::ParameterList& params);

  %extend {
    NOXSolver(PyObject* dict)
    {
      Teuchos::ParameterList params = dict2ParameterList(dict);
      return new NOXSolver(params);
    }
    }
  
  /** */
  SolverState<double> solve(const NonlinearOperator<double>& F, 
    Vector<double>& soln) const ;

  /** */
  const LinearSolver<double>& linSolver() const ;

};

%extend NOXSolver {
  NOXSolver(PyObject* dict)
  {
    Teuchos::ParameterList params = dict2ParameterList(dict);
    return NOXSolver(params);
  }
}
}



%inline %{
  /* Create a nonlinear solver from a Python dictionary */
  Playa::NOXSolver makeNOXSolver(PyObject* dict)
  {
    Teuchos::ParameterList params = dict2ParameterList(dict);
    return NOXSolver(params);
  }
  /* Create a nonlinear solver from a Python dictionary */
  Playa::NOXSolver makeNOXSolver(const Teuchos::ParameterList& params)
  {
    return NOXSolver(params);
  }
  %}




%inline %{
  namespace Playa
  {
    SolverState<double> 
    PySundanceLinearSolver_solve(const PySundanceLinearSolver* solver,
                                 const LinearOperator<double>& op,
                                 const Vector<double>& rhs,
                                 Vector<double>& soln)
    {
      swig_type_info* opType = SWIG_TypeQuery("Playa::LinearOperator<double>*");
      TEUCHOS_TEST_FOR_EXCEPTION(opType==0, runtime_error,
                         "swig could not find a match for type name "
                         "[Playa::LinearOperator<double>]");


      swig_type_info* vecType = SWIG_TypeQuery("Playa::Vector<double>*");
      TEUCHOS_TEST_FOR_EXCEPTION(vecType==0, runtime_error,
                         "swig could not find a match for type name "
                         "[Playa::Vector<double>]");


      swig_type_info* stateType = SWIG_TypeQuery("Playa::SolverState<double>*");
      TEUCHOS_TEST_FOR_EXCEPTION(stateType==0, runtime_error,
                         "swig could not find a match for type name "
                         "[Playa::SolverState<double>]");


      PyObject* opObj = SWIG_NewPointerObj( (void*) &op, opType, 0);
      PyObject* rhsObj = SWIG_NewPointerObj( (void*) &rhs, vecType, 0);
      PyObject* x0Obj = SWIG_NewPointerObj( (void*) &soln, vecType, 0);

      PyObject* result = solver->pySolve(opObj, rhsObj, x0Obj);

      if (0 == result) {
        PyErr_Print();
        return SolverState<double>(SolveCrashed, "null result from PySundanceLinearSolver",
                                   1, 0.0);
      }

      PyObject* solnObj = 0;
      PyObject* stateObj = 0 ;

      Vector<double>* x = 0 ;
      SolverState<double>* state = 0 ;

      int isTuple = PyTuple_Check(result);

      if (isTuple)
        {
          int size = PyTuple_Size(result);
          switch(size)
            {
            case 2:
              stateObj = PyTuple_GetItem(result, 1);
              TEUCHOS_TEST_FOR_EXCEPTION(stateObj==0, runtime_error,
                                 "null solver state in PySundanceLinearSolver_solve()");
              SWIG_Python_ConvertPtr(stateObj, (void**) &state, stateType,  
                                     SWIG_POINTER_EXCEPTION | 0);
            case 1:
              solnObj = PyTuple_GetItem(result, 0);
              TEUCHOS_TEST_FOR_EXCEPTION(solnObj==0, runtime_error,
                                 "null solution object in PySundanceLinearSolver_solve()");
              SWIG_Python_ConvertPtr(solnObj, (void**) &x, vecType,  
                                     SWIG_POINTER_EXCEPTION | 0);
              break;
            default:
              TEUCHOS_TEST_FOR_EXCEPTION(size < 1 || size > 2, runtime_error,
                                 "invalid return value size " << size 
                                 << " in PySundanceLinearSolver_solve()");
            }
        }
      else
        {
          SWIG_Python_ConvertPtr(result, (void**) &x, vecType,  
                                 SWIG_POINTER_EXCEPTION | 0);
        }

      TEUCHOS_TEST_FOR_EXCEPTION(x==0, runtime_error, "null return vector in "
                         " PySundanceLinearSolver_solve()");
      soln = *x;

      SolverState<double> rtn(SolveConverged, "unknown solve state", 1, 0);
      if (state!=0) rtn = *state;
      

      Py_DECREF(result); // All done with returned result object

      return rtn;
    }
  }
  %}



