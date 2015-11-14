#include "Mesquite_ShapeImprover.hpp"
#include "Mesquite_MeshImpl.hpp"
#include "Mesquite_CLArgs.hpp"
#include "Mesquite_MsqError.hpp"

#include "Mesquite_domain.hpp"

#include <iostream>

using namespace Mesquite;

int main( int argc, char* argv[] )
{
  const double zero = 0.0, one = 1.0;
  double default_cpu_time = 0.0;
  CLArgs::DoubleRangeArg movement_beta( &zero, &one ), cpu_time( default_cpu_time, &zero, 0 );
  
  CLArgs args( "msqshape",
               "Run Shape Improvement smoother for input mesh.",
               "Read VTK file, smooth, and re-write file." );
  args.double_flag( 'b', "Vtx Movement Beta", "fraction of mean edge length to use as termination criterion on vertex movement", &movement_beta );
  args.double_flag( 't', "Cpu Seconds", "time-out", &cpu_time );
                    
  add_domain_args( args );
  args.add_required_arg( "input_file" );
  args.add_required_arg( "output_file" );

  std::vector<std::string> files;
  if (!args.parse_options( argc, argv, files, std::cerr )) {
    args.print_usage( std::cerr );
    exit(1);
  }
  std::string input_file = files[0];
  std::string output_file = files[1];
  
  MsqError err;
  MeshImpl mesh;
  mesh.read_vtk( input_file.c_str(), err );
  if (err) {
    std::cerr << "ERROR READING FILE: " << input_file << std::endl
                    << err << std::endl;
    return 2;
  }
  MeshDomain* domain = process_domain_args( &mesh );

  ShapeImprover smoother;
  if (movement_beta.seen())
    smoother.set_vertex_movement_limit_factor( movement_beta.value() );
  if (cpu_time.seen())
    smoother.set_cpu_time_limit( cpu_time.value() );
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, domain);
  smoother.run_instructions( &mesh_and_domain, err );
  if (err) {
    std::cerr << err << std::endl;
    return 3;
  }
  
  mesh.write_vtk( output_file.c_str(), err );
  if (err) {
    std::cerr << "ERROR WRITING FILE: " << output_file << std::endl
                    << err << std::endl;
    return 2;
  }
  
  return 0;
}

