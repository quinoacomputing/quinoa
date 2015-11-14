/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2010 Sandia National Laboratories.  Developed at the
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
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA

    (2010) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file untangle.cpp
 *  \brief tast untangle wrapper
 *  \author Jason Kraftcheck 
 */


#include "meshfiles.h"

#include "Mesquite_MeshImpl.hpp"
#include "Mesquite_MsqError.hpp"
#include "Mesquite_MsqVertex.hpp"
#include "Mesquite_PlanarDomain.hpp"
#include "Mesquite_UntangleWrapper.hpp"
#include "Mesquite_QualityAssessor.hpp"

#include <iostream>
using std::cout;
using std::endl;
#include <cstdlib>

using namespace Mesquite;

#define VTK_2D_DIR MESH_FILES_DIR "2D/vtk/"

// Test untangle wrapper
// Assumes all meshes lie in a plane for which the normal is [0,0,1].
int uwt( bool skip,
         UntangleWrapper::UntangleMetric metric,
         const char* input_file,
         int expected_number_of_remaining_inverted_elems,
         bool flip_domain = false );

bool brief_output = false;
bool write_output = false;
double mu_sigma = -1;
double beta = -1;

void usage( const char* argv0 ) {
  std::cerr << "Usage: " << argv0 << " [-<flags>] [-c <sigma>] [-b <beta>]" << std::endl
            << "       " << argv0 << " -h" << std::endl;
}
void help( const char* argv0 ) {
  std::cout << "Usage: " << argv0 << " [-<flags>] [-c <sigma>] [-b <beta>]" << std::endl
            << "Flags: -q : brief output" << std::endl
            << "       -w : write result meshes" << std::endl
            << "       -B : skip tests using untangle beta target metric" << std::endl
            << "       -Z : skip tests using size untangle target metric" << std::endl
            << "       -P : skip tests using shapesize untangle target metric" << std::endl
            << "       -H : skip tests using 'tangled_horse1.vtk' as input" << std::endl
            << "       -Q : skip tests using 'hole_in_square_tanglap.vtk' as input" << std::endl
            << "       -I : skip tests using 'inverted-hole-2.vtk' as input" << std::endl
            << "       -S : skip tests using 'shest_grid32.vtk' as input" << std::endl
            << "       -c : specify sigma value for untangle metrics" << std::endl
            << "       -b : specify beta value for untangle beta metric" << std::endl
            << std::endl;
}

int main( int argc, char* argv[] )
{
  bool skip_beta = false;
  bool skip_size = false;
  bool skip_shape = false;
  bool skip_horse = false;
  bool skip_hole = false;
  bool skip_invrt = false;
  bool skip_shest = false;
  std::list<double*> expected;

  for (int i = 1; i < argc; ++i) {
    if (!expected.empty()) {
      char* endptr;
      *expected.front() = strtod( argv[i], &endptr );
      if (*endptr || *expected.front() <= 0) {
        std::cerr << "Expected positive number, found \"" << argv[i] << '"' << std::endl;
        usage(argv[0]);
        return 1;
      }
      expected.pop_front();
    }
    else if (argv[i][0] == '-' && argv[i][1]) {
      for (int j = 1; argv[i][j]; ++j) {
        switch (argv[i][j]) {
          case 'q': brief_output = true; break;
          case 'w': write_output = true; break;
          case 'B': skip_beta = true; break;
          case 'Z': skip_size = true; break;
          case 'P': skip_shape = true; break;
          case 'H': skip_horse = true; break;
          case 'Q': skip_hole = true; break;
          case 'I': skip_invrt = true; break;
          case 'S': skip_shest = true; break;
          case 'c': expected.push_back(&mu_sigma); break;
          case 'b': expected.push_back(&beta); break;
          case 'h': help(argv[0]); return 0;
          default:
            std::cerr << "Invalid flag: -" << argv[i][j] << std::endl;
            usage(argv[0]);
            return 1;
        }
      }
    }
    else {
      std::cerr << "Unexpected argument: \"" << argv[i] << '"' << std::endl;
      usage(argv[0]);
      return 1;
    }
  }

  int result = 0;
  
  result += uwt( skip_beta||skip_horse, UntangleWrapper::BETA, "quads/tangled/tangled_horse1.vtk", 0 );
  result += uwt( skip_beta||skip_hole , UntangleWrapper::BETA, "quads/tangled/hole_in_square_tanglap.vtk", 0, true );
  result += uwt( skip_beta||skip_invrt, UntangleWrapper::BETA, "quads/tangled/inverted-hole-2.vtk", 0 );
  result += uwt( skip_beta||skip_shest, UntangleWrapper::BETA, "quads/untangled/shest_grid32.vtk", 0 );

  result += uwt( skip_size||skip_horse, UntangleWrapper::SIZE, "quads/tangled/tangled_horse1.vtk", 0 );
  result += uwt( skip_size||skip_hole , UntangleWrapper::SIZE, "quads/tangled/hole_in_square_tanglap.vtk", 6, true );
  result += uwt( skip_size||skip_invrt, UntangleWrapper::SIZE, "quads/tangled/inverted-hole-2.vtk", 0  );
  result += uwt( skip_size||skip_shest, UntangleWrapper::SIZE, "quads/untangled/shest_grid32.vtk", 0 );

  result += uwt( skip_shape||skip_horse, UntangleWrapper::SHAPESIZE, "quads/tangled/tangled_horse1.vtk", 0 );
  result += uwt( skip_shape||skip_hole , UntangleWrapper::SHAPESIZE, "quads/tangled/hole_in_square_tanglap.vtk", 0, true );
  result += uwt( skip_shape||skip_invrt, UntangleWrapper::SHAPESIZE, "quads/tangled/inverted-hole-2.vtk", 8  );
  result += uwt( skip_shape||skip_shest, UntangleWrapper::SHAPESIZE, "quads/untangled/shest_grid32.vtk", 0 );
  
  return result;
}

const char* tostr( UntangleWrapper::UntangleMetric m )
{
  static const char BETA[] = "BETA";
  static const char SIZE[] = "SIZE";
  static const char SHAPESIZE[] = "SHAPESIZE";
  switch (m) {
    case UntangleWrapper::BETA:      return BETA;
    case UntangleWrapper::SIZE:      return SIZE;
    case UntangleWrapper::SHAPESIZE: return SHAPESIZE;
  }
  return 0;
}

int uwt( bool skip,
         UntangleWrapper::UntangleMetric metric,
         const char* input_file_base,
         int expected,
         bool flip_domain )
{
  if (skip)
    return 0;

  if (!brief_output)
    std::cout << std::endl << "**********************************************" << std::endl;
  std::cout << "Running \"" << input_file_base << "\" for " << tostr(metric) << std::endl;
  if (!brief_output)
    std::cout << "**********************************************" << std::endl << std::endl;

    // get mesh
  MsqError err;
  MeshImpl mesh;
  std::string input_file( VTK_2D_DIR );
  input_file += input_file_base;
  mesh.read_vtk( input_file.c_str(), err );
  if (err) {
    std::cerr << err << std::endl;
    std::cerr << "ERROR: " << input_file << " : failed to read file" << std::endl;
    return 1;
  }
    // get domain
  std::vector<Mesh::VertexHandle> verts;
  mesh.get_all_vertices( verts, err );
  if (err || verts.empty()) abort();
  MsqVertex coords;
  mesh.vertices_get_coordinates( arrptr(verts), &coords, 1, err );
  if (err) abort();
  Vector3D norm(0,0,flip_domain ? -1 : 1);
  PlanarDomain domain( norm, coords );
    // run wrapper
  UntangleWrapper wrapper( metric );
  wrapper.set_vertex_movement_limit_factor( 0.005 );
  double constant = (metric == UntangleWrapper::BETA) ? beta : mu_sigma;
  if (constant > 0)
    wrapper.set_metric_constant( constant );
  if (brief_output)
    wrapper.quality_assessor().disable_printing_results();
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, &domain);
  wrapper.run_instructions( &mesh_and_domain, err );
  if (err) {
    std::cerr << err << std::endl;
    std::cerr << "ERROR: optimization failed" << std::endl;
    return 1;
  }
    // write output file
  if (write_output) {
    std::string result_file(tostr(metric));
    result_file += "-";
    result_file += input_file_base;
    mesh.write_vtk( result_file.c_str(), err );
    if (err) {
      std::cerr << err << std::endl;
      std::cerr << "ERROR: " << result_file << " : failed to write file" << std::endl;
      err.clear();
    }
    else {
      std::cerr << "Wrote file: " << result_file << std::endl;
    }
  }
  
    // test number of inverted elements
  int count, junk;
  wrapper.quality_assessor().get_inverted_element_count( count, junk, err );
  if (err) abort();
  if (count < expected) {
    std::cout << "WARNING: expected " << expected 
              << " inverted elements but finished with only " 
              << count << std::endl
              << "Test needs to be updated?" << std::endl << std::endl;
    return 0;
  }
  else if (count == expected) {
    std::cout << "Completed with " << count << " inverted elements remaining" 
              << std::endl << std::endl;
    return 0;
  }
  else {
    std::cerr << "ERROR: expected " << expected 
              << " inverted elements but finished with " 
              << count << std::endl << std::endl;
    return 1;
  }
}
