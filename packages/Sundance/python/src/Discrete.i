// -*- c++ -*-

%{
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceBlock.hpp"
#include "SundanceL2Projector.hpp"
#include "SundanceDiscreteSpace.hpp"
#include "SundanceDiscreteFunction.hpp"
  %}


// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"



%{
  Sundance::BasisArray pyListToBasisArray(PyObject* lst)
  {
    PyObject_Print(lst, stderr, Py_PRINT_RAW);
    TEUCHOS_TEST_FOR_EXCEPTION(!PyList_Check(lst), std::runtime_error, 
                       "Expecting a python list as argument to conversion to basis array");
    int n = PyList_Size(lst);
    Sundance::BasisArray rtn(n);
    
    for (int i=0; i<n; i++)
    {
      PyObject *obj_i = PyList_GetItem(lst,i);
      Sundance::BasisFamily *basis_i = 0;
      SWIG_Python_ConvertPtr(obj_i, (void**) &basis_i, 
                             SWIGTYPE_p_Sundance__BasisFamily, 
                             SWIG_POINTER_EXCEPTION | 0);
      rtn[i] = *basis_i;
    }

    return rtn;
  }

  Sundance::Mesh pyObjToMesh(PyObject* obj)
  {
    Sundance::Mesh rtn;
    Sundance::Mesh* meshPtr = 0;
    SWIG_Python_ConvertPtr(obj, (void**) &meshPtr, 
                             SWIGTYPE_p_Sundance__Mesh, 
                             SWIG_POINTER_EXCEPTION | 0);
    rtn = *meshPtr;
    return rtn;
  }

  Playa::VectorType<double> pyObjToVectorType(PyObject* obj)
  {
    Playa::VectorType<double>  rtn;
    Playa::VectorType<double> * vPtr = 0;
    SWIG_Python_ConvertPtr(obj, (void**) &vPtr, 
                             SWIGTYPE_p_Playa__VectorTypeTdouble_t, 
                             SWIG_POINTER_EXCEPTION | 0);
    rtn = *vPtr;
    return rtn;
  }

%}

/*
%typemap(in) (const Sundance::BasisArray& basis)(Sundance::BasisArray basis)
{
  std::cerr << "in basis array typemap" << std::endl;
  basis = pyListToBasisArray($input);

  $1 = &basis;
}



%typemap(in) (const Sundance::Mesh& mesh, 
              const Sundance::BasisArray& basis,
              const Playa::VectorType<double>& vecType)
  (Sundance::Mesh mesh,
   Sundance::BasisArray basis,
   Playa::VectorType<double> vecType)
{
  std::cerr << "in (mesh, basis, vecType) typemap" << std::endl;
  TEUCHOS_TEST_FOR_EXCEPTION(!PyTuple_Check($input), std::runtime_error,
                     "expecting a tuple");
  TEUCHOS_TEST_FOR_EXCEPTION(PyTuple_Size($input) != 3, std::runtime_error,
                     "expecting a tuple of length 3");
  mesh = pyObjToMesh(PyTuple_GetItem($input, 0));
  basis = pyListToBasisArray(PyTuple_GetItem($input, 1));
  vecType = pyObjToVectorType(PyTuple_GetItem($input, 2));
  
  $1 = &mesh;
  $2 = &basis;
  $3 = &vecType;
}
*/
 


namespace Sundance
{

  


  class DiscreteSpace
  {
  public:
    /* */
    DiscreteSpace(const Sundance::Mesh& mesh, 
                  const BasisFamily& basis,
                  const Playa::VectorType<double>& vecType);
    /* */
    DiscreteSpace(const Sundance::Mesh& mesh, 
                  const Sundance::BasisArray& basis,
                  const Playa::VectorType<double>& vecType);
    /* */
    DiscreteSpace(const Sundance::Mesh& mesh, 
                  const Sundance::BasisArray& basis,
                  const Sundance::CellFilterArray& domains,
                  const Playa::VectorType<double>& vecType);
    /* */
    DiscreteSpace(const Sundance::Mesh& mesh, 
                  const BasisFamily& basis,
                  const SpectralBasis& sb,
                  const Playa::VectorType<double>& vecType);
    /** */
    DiscreteSpace(const Sundance::Mesh& mesh, 
                  const BasisFamily& basis,
                  const CellFilter& regions,
                  const Playa::VectorType<double>& vecType);


    /** */
    DiscreteSpace(const Sundance::Mesh& mesh, 
                  const BasisArray& basis,
                  const CellFilter& regions,
                  const Playa::VectorType<double>& vecType);
    /* */
    DiscreteSpace(const Sundance::Mesh& mesh, 
                  const Sundance::BasisArray& basis,
                  const SpectralBasis& sb,
                  const Playa::VectorType<double>& vecType);
    /* */
    ~DiscreteSpace();

    /* */
    const Sundance::Mesh& mesh() const ;

    /* */
    Playa::VectorSpace<double> vecSpace() const ;

    /* */
    Playa::VectorType<double> vecType() const ;

    
  };

  class L2Projector
  {
  public:
    /* */
    L2Projector(const DiscreteSpace& space, 
                const Sundance::Expr& expr);
    /* */
    L2Projector(const DiscreteSpace& space, 
                const Sundance::Expr& expr,
                const Playa::LinearSolver<double>& solver);

    /* */
    ~L2Projector();

    /* */
    Sundance::Expr project() const ;

    /* */
    const LinearProblem& prob() const ;
    
  };
}



%inline %{
  void printVecBasis(const Sundance::BasisArray& basis)
  {
    std::cerr << "vector basis = " << basis << std::endl;
  }
  %}


/*
%inline %{
  void printVecBasis(int i, const Sundance::BasisArray& basis)
  {
    std::cerr << i << " vector basis = " << basis << std::endl;
  }
  %}

%inline %{
  void printVecBasis(const Sundance::BasisArray& basis, int i)
  {
    std::cerr << "vector basis = " << basis << " " << i << std::endl;
  }
  %}
*/





%rename(DiscreteFunction) makeDiscreteFunction;

%inline %{
  /* Create a discrete function */
  Sundance::Expr makeDiscreteFunction(const Sundance::DiscreteSpace& space,
                                          const Playa::Vector<double>& vec)
  {
    return new Sundance::DiscreteFunction(space, vec);
  }
  %}

%inline %{
  /* Create a discrete function */
  Sundance::Expr makeDiscreteFunction(const Sundance::DiscreteSpace& space,
                                          const double& val)
  {
    return new Sundance::DiscreteFunction(space, val);
  }
  %}


%inline %{
  /* Create a discrete function */
  Sundance::Expr makeDiscreteFunction(const Sundance::DiscreteSpace& space,
                                          const Playa::Vector<double>& vec,
                                          const std::string& name)
  {
    return new Sundance::DiscreteFunction(space, vec, name);
  }
  %}

%inline %{
  /* Create a discrete function */
  Sundance::Expr makeDiscreteFunction(const Sundance::DiscreteSpace& space,
                                          const double& val,
                                          const std::string& name)
  {
    return new Sundance::DiscreteFunction(space, val, name);
  }
  %}



