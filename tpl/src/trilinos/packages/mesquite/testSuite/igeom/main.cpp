/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2009 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    (2010) kraftche@cae.wisc.edu    

  ***************************************************************** */





#include "meshfiles.h"

#include "Mesquite.hpp"
#include "MsqIGeom.hpp"
#include "MeshImpl.hpp"
#include "MsqError.hpp"
#include "ShapeImprovementWrapper.hpp"
#include "MsqVertex.hpp"
#include "QualityAssessor.hpp"
#include "SphericalDomain.hpp"

#include <memory>
#include <iostream>

#include "MsqIBase.hpp"

using namespace Mesquite;

// characteristics of geometry
const double SPHERE_RADIUS = 3.0;
const Vector3D SPHERE_CENTER( 2.0, 2.0, 0.0 );

#define CHKIGEOM \
  if (chk_igeom_error(ierr, __FILE__, __LINE__)) return 0
  
bool chk_igeom_error( int ierr, const char* file, int line )
{
  if (iBase_SUCCESS == ierr) return false;
  std::cerr << "iGeom call failed at " << file << ":" << line << std::endl;
  std::cerr << process_itaps_error(ierr) << std::endl;
  return true;
}

const char* const default_file_name = MESH_FILES_DIR "2D/vtk/quads/untangled/quads_on_sphere_529.vtk";

void usage()
{
  std::cout << "main [-N] [filename] [-o <output_file>]" << std::endl;
  std::cout << "  -N : Use native representation instead of iGeom implementation\n" << std::endl;
  std::cout << "  If no file name is specified, will use \"" 
       << default_file_name << '"' << std::endl;
  exit (1);
}

void run_smoother( Mesh& mesh, MeshDomain* domain, MsqError& err );
bool check_results( Mesh& mesh, MeshDomain* domain, MsqError& err );
MeshDomain* get_itaps_domain();
MeshDomain* get_mesquite_domain();

int main(int argc, char* argv[])
{
    // command line arguments
  const char* file_name = 0;
  const char* output_file = 0;
  bool use_native = false, opts_done = false;
  for (int arg = 1; arg < argc; ++arg)
  {
    if (!opts_done && argv[arg][0] == '-')
    {
      if (!strcmp( argv[arg], "-N"))
        use_native = true;
      else if (!strcmp( argv[arg], "-o")) 
        output_file = argv[arg++];
      else if (!strcmp( argv[arg], "--"))
        opts_done = true;
      else
        usage();
    }
    else if (!file_name)
      file_name = argv[arg];
    else
      usage();
  }
  if (!file_name)
  {
    file_name = default_file_name;
    std::cout << "No file specified: using default: " << default_file_name << std::endl;
  }  
  
  
  MsqError err;
  MeshImpl mesh;
  mesh.read_vtk( file_name, err );
  if (MSQ_CHKERR(err)) {
    std::cerr << err << std::endl;
    std::cerr << "Failed to read input file: " << file_name << std::endl;
    return 1;
  }

  std::auto_ptr<MeshDomain> dom(use_native ? 
                                get_mesquite_domain() :
                                get_itaps_domain());

  if (!dom.get()) {
    MSQ_SETERR(err)("Domain creation failed", MsqError::INTERNAL_ERROR);
    std::cerr << err << std::endl;
    return 1;
  }
  
  run_smoother( mesh, dom.get(), err );
  if (MSQ_CHKERR(err)) {
    std::cerr << err << std::endl;
    std::cerr << "Smoother failed." << std::endl;
    return 1;
  }
  
  if (output_file) {
    mesh.write_vtk( output_file, err );
    if (MSQ_CHKERR(err)) {
      std::cerr << err << std::endl;
      std::cerr << "Failed to write file: " << output_file << std::endl;
      return 1;
    }
  }
  
  bool valid = check_results( mesh, dom.get(), err );
  if (MSQ_CHKERR(err)) {
    std::cerr << err << std::endl;
    std::cerr << "Error while validating results." << std::endl;
    return 1;
  }
  
  return !valid;
}

void run_smoother( Mesh& mesh, MeshDomain* dom, MsqError& err )
{
  ShapeImprovementWrapper smoother;
  smoother.run_instructions( &mesh, dom, err );
  MSQ_CHKERR(err);
  
  if (smoother.quality_assessor().invalid_elements()) {
    MSQ_SETERR(err)("Smoothing produced invalid elements", MsqError::INVALID_MESH);
    return;
  }
}


bool check_results( Mesh& mesh, MeshDomain* dom, MsqError& err )
{
  double EPS = 1e-4;

  std::vector<Mesh::VertexHandle> handles;
  mesh.get_all_vertices( handles, err ); MSQ_ERRZERO(err);
  std::vector<MsqVertex> coords(handles.size());
  mesh.vertices_get_coordinates( arrptr(handles), arrptr(coords), handles.size(), err );
  MSQ_ERRZERO(err);
  
  bool valid = true;
  size_t c = 0;
  for (size_t i = 0; i < coords.size(); ++i) {
    double d = (coords[i] - SPHERE_CENTER).length();
    if (fabs(d - SPHERE_RADIUS) > EPS) 
      ++c;
  }
  if (c) {
    std::cerr << c << " vertices not on domain" << std::endl;
    valid = false;
  }
  
  std::vector<unsigned short> dof(handles.size()), exp_dof(handles.size(),2);
  dom->domain_DoF( arrptr(handles), arrptr(dof), handles.size(), err ); MSQ_ERRZERO(err);
  if (dof != exp_dof) {
    std::cerr << "Invalid domain dimension for one or more vertices" << std::endl;
    valid = false;
  }
  
  c = 0;
  std::vector<Vector3D> normals( coords.begin(), coords.end() );
  dom->vertex_normal_at( arrptr(handles), arrptr(normals), handles.size(), err );
  MSQ_ERRZERO(err);
  for (size_t i = 0; i < handles.size(); ++i) {
    Vector3D exp_norm = ~(coords[i] - SPHERE_CENTER);
    Vector3D act_norm = ~normals[i];
    Vector3D diff = exp_norm - act_norm;
    if (diff.length() > EPS) 
      ++c;
  }
  if (c) {
    std::cerr << c << " invalid vertex normals" << std::endl;
    valid = false;
  }
  
  return valid;
}


MeshDomain* get_itaps_domain()
{
  const double EPS = 1e-6;
  int ierr;
  iGeom_Instance igeom;
  iGeom_newGeom( "", &igeom, &ierr, 0 ); CHKIGEOM;
  
  iBase_EntityHandle sphere_vol;
  iGeom_createSphere( igeom, SPHERE_RADIUS, &sphere_vol, &ierr ); CHKIGEOM;
  iGeom_moveEnt( igeom, sphere_vol, 
                 SPHERE_CENTER[0], SPHERE_CENTER[1], SPHERE_CENTER[2],
                 &ierr ); CHKIGEOM;
  
  iBase_EntityHandle sphere_surf;
  iBase_EntityHandle* ptr = &sphere_surf;
  int size = 0, alloc = 1;
  iGeom_getEntAdj( igeom, sphere_vol, iBase_FACE, &ptr, &alloc, &size, &ierr );
  CHKIGEOM;
  if (size != 1) {
    std::cerr << "Failed to get spherical surface from iGeom instance" << std::endl;
    return 0;
  }
  
  Vector3D bmin, bmax;
  iGeom_getEntBoundBox( igeom, sphere_surf, &bmin[0], &bmin[1], &bmin[2],
                        &bmax[0], &bmax[1], &bmax[2], &ierr ); CHKIGEOM;
  Vector3D center = 0.5 * (bmin + bmax);
  Vector3D rad = 0.5 * (bmax - bmin);
  if ((center - SPHERE_CENTER).length() > EPS ||
      fabs(rad[0] - SPHERE_RADIUS) > EPS ||
      fabs(rad[1] - SPHERE_RADIUS) > EPS ||
      fabs(rad[2] - SPHERE_RADIUS) > EPS) {
    std::cerr << "iGeom implementation returned invalid box for sphere" << std::endl
              << "  Expected Sphere Center: " << SPHERE_CENTER << std::endl
              << "  Expected Sphere Radius: " << SPHERE_RADIUS << std::endl
              << "  Actual Box Minimum    : " << bmin << std::endl
              << "  Actual Box Maximum    : " << bmax << std::endl;
    return 0;
  }
  
  return new MsqIGeom( igeom, sphere_surf );
}

MeshDomain* get_mesquite_domain()
{
  return new SphericalDomain( SPHERE_CENTER, SPHERE_RADIUS );
}

