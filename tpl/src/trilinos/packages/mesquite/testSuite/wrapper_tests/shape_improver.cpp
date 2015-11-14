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
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA

    (2009) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file main.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifdef TEST_OLD_WRAPPER
#  include "Mesquite_ShapeImprovementWrapper.hpp"
#else
#  include "Mesquite_ShapeImprover.hpp"
#endif
#include "Mesquite_QualityAssessor.hpp"
#include "Mesquite_MeshImpl.hpp"
#include "Mesquite_MsqError.hpp"
#include "Mesquite_PlanarDomain.hpp"
#include "Mesquite_IdealWeightInverseMeanRatio.hpp"
#include "meshfiles.h"

#include <iostream>
#include <stdlib.h>
#include <string.h>

using namespace Mesquite;


const char DEFAULT_INPUT[] = MESH_FILES_DIR "/2D/vtk/quads/tangled/inverted-hole-1.vtk";

void usage( const char* argv0, bool help = false ) {
  std::ostream& str = help ? std::cout : std::cerr;
  str << "Usage: " << argv0 << "[input_file] [output_file]" << std::endl;
  str << "Usage: " << argv0 << "-h" << std::endl;
  if (help) {
    str << "Default input file: " << DEFAULT_INPUT << std::endl
        << "Default output file: (none)" << std::endl
        << "Input meshes are expected to lie in the XY-plane" << std::endl
        << "Input meshes are expected to be smoothable to ideal elements" << std::endl;
  }
  exit(!help);
}
    
    
int main( int argc, char* argv[] )
{
  const char* input_file = 0;
  const char* output_file = 0;
  for (int i = 1; i < argc; ++i) {
    if (!strcmp("-h", argv[i]))
      usage(argv[0],true);
    else if (!input_file)
      input_file = argv[i];
    else if (!output_file)
      output_file = argv[i];
    else
      usage(argv[0]);
  }
  if (!input_file)
    input_file = DEFAULT_INPUT;
  
  MsqError err;
  MeshImpl mesh;
  mesh.read_vtk( input_file, err );
  if (err) {
    std::cerr << err << std::endl
              << input_file << ": failed to read file" << std::endl;
    return 3;
  }
  
  PlanarDomain plane(PlanarDomain::XY);
#ifdef TEST_OLD_WRAPPER
  ShapeImprovementWrapper smoother;
#else
  ShapeImprover smoother;
#endif
  IdealWeightInverseMeanRatio extra_metric;
  smoother.quality_assessor().add_quality_assessment(&extra_metric);
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, &plane);
  smoother.run_instructions( &mesh_and_domain, err );
  if (err) {
    std::cerr << err << std::endl
              << input_file << ": smoother failed" << std::endl;
    return 2;
  }
  
  if (output_file) {
    mesh.write_vtk( output_file, err );
    if (err) {
      std::cerr << err << std::endl
                << output_file << ": failed to write file" << std::endl;
      return 3;
    }
  }

  if (smoother.quality_assessor().invalid_elements()) {
    std::cerr << "Resulting mesh contains invalid elements: untangler did not succeed" << std::endl;
    return 4;
  }
  
  const QualityAssessor::Assessor* quality = 
    smoother.quality_assessor().get_results(&extra_metric);
  if (!quality) {
    std::cerr << "Failed to get quality stats for IMR metric" << std::endl;
    return 2;
  }
  
  if (fabs(1 - quality->get_average()) > 1e-3) {
    std::cerr << "Average quality is not optimal." << std::endl;
    return 4;
  }
  
  if (quality->get_stddev() > 1e-3) {
    std::cerr << "Not all elements have optimal quality." << std::endl;
    return 4;
  }
  
  return 0;
}

