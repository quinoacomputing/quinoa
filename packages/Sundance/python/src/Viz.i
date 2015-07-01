// -*- c++ -*-

%{
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceFieldWriter.hpp"
#include "SundanceVTKWriter.hpp"
#include "SundanceTriangleWriter.hpp"
#include "SundanceMatlabWriter.hpp"
#include "SundanceExprFieldWrapper.hpp"

  %}




// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"


namespace Sundance
{
  class FieldWriter
  {
  public:
    FieldWriter();
    ~FieldWriter();

    void addMesh(const Mesh& mesh);

    %extend 
    {
      void addField(const std::string& name, const Sundance::Expr& f)
      {
        self->addField(name, new Sundance::ExprFieldWrapper(f));
      }
    }

    void write();

    void setUndefinedValue(const double& x);

  };

}

%rename(VTKWriter) makeVTKWriter;
%rename(MatlabWriter) makeMatlabWriter;
%rename(TriangleWriter) makeTriangleWriter;


%inline %{
  /* Create a VTK writer */
  Sundance::FieldWriter makeVTKWriter(const std::string& filename)
  {
    return new Sundance
      ::VTKWriter(filename);
  }
  %}

%inline %{
  /* Create a Triangle writer */
  Sundance::FieldWriter makeTriangleWriter(const std::string& filename)
  {
    return new Sundance
      ::TriangleWriter(filename);
  }
  %}

%inline %{
  /* Create a Matlab writer */
  Sundance::FieldWriter makeMatlabWriter(const std::string& filename)
  {
    return new Sundance
      ::MatlabWriter(filename);
  }
  %}



