// -*- c++ -*-

%{
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceMesh.hpp"
#include "SundanceMeshSource.hpp"
#include "SundancePartitionedLineMesher.hpp"
#include "SundancePartitionedRectangleMesher.hpp"
#include "SundanceTriangleMeshReader.hpp"
#include "SundanceExodusNetCDFMeshReader.hpp"
#include "SundanceBasicSimplicialMeshType.hpp"
#include "SundanceGeomUtils.hpp"

  %}




// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"


namespace Sundance
{
  class Point
  {
  public:
    Point(const double& x);
    Point(const double& x, const double& y);
    Point(const double& x, const double& y, const double& z);

    int dim();

    
  };

  %extend Point
  {
    using namespace std;
    std::string __str__() 
    {
      return self->toString();
    }
  }
}

namespace Sundance
{
  class Mesh
  {
  public:
    Mesh();
    ~Mesh();
    int spatialDim() const ;
    int numCells(int dim) const ;
    void dump(const std::string& str) const ;
    bool checkConsistency(const std::string& str) const ;
  };

  %extend Mesh
  {
    using namespace std;
    std::string __str__() 
    {
      std::string rtn; 
      rtn = "Mesh[dim=" + Teuchos::toString(self->spatialDim()) + "]";
      return rtn;
    }
  }

  class MeshSource
  {
  public:
    MeshSource();
    ~MeshSource();
    Mesh getMesh() const ;
  };

}

%rename(PartitionedLineMesher) makePartitionedLineMesher;
%rename(PartitionedRectangleMesher) makePartitionedRectangleMesher;
%rename(ExodusNetCDFMeshReader) makeExodusNetCDFMeshReader;
%rename(ExodusMeshReader) makeExodusMeshReader;
%rename(TriangleMeshReader) makeTriangleMeshReader;

%inline %{
  /* Create a line mesher */
  Sundance::MeshSource makePartitionedLineMesher(double ax, 
                                  double bx, 
                                  int nx)
  {
    return new Sundance
      ::PartitionedLineMesher(ax, bx, nx, 
                              new Sundance::BasicSimplicialMeshType());
  }
  %}

%inline %{
  /* Create a rectangle mesher */
  Sundance::MeshSource makePartitionedRectangleMesher(double ax, double bx, 
                                                             int nx, int npx,
                                                             double ay, double by,
                                                             int ny, int npy)
  {
    return new Sundance
      ::PartitionedRectangleMesher(ax, bx, nx, npx, ay, by, ny, npy,
                                   new Sundance::BasicSimplicialMeshType());
  }
  %}


%inline %{
  /* Create an exodus reader */
  Sundance::MeshSource makeExodusNetCDFMeshReader(const std::string& filename)
  {
    return new Sundance
      ::ExodusNetCDFMeshReader(filename,
                               new Sundance::BasicSimplicialMeshType());
  }
  %}


%inline %{
  /* Create an exodus reader */
  Sundance::MeshSource makeExodusMeshReader(const std::string& filename)
  {
    return new Sundance
      ::ExodusMeshReader(filename,
        new Sundance::BasicSimplicialMeshType());
  }
  %}




%inline %{
  /* Create a triangle reader */
  Sundance::MeshSource makeTriangleMeshReader(const std::string& filename)
  {
    return new Sundance
      ::TriangleMeshReader(filename,
                               new Sundance::BasicSimplicialMeshType());
  }
  %}

