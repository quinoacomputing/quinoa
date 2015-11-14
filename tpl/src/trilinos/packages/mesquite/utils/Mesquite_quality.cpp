#include "Mesquite_CLArgs.hpp"
#include "Mesquite_QualityAssessor.hpp"
#include "Mesquite_IdealWeightInverseMeanRatio.hpp"
#include "Mesquite_SizeMetric.hpp"
#include "Mesquite_TQualityMetric.hpp"
#include "Mesquite_IdealShapeTarget.hpp"
#include "Mesquite_TShapeNB1.hpp"
#include "Mesquite_InstructionQueue.hpp"
#include "Mesquite_MsqError.hpp"
#include "Mesquite_MeshImpl.hpp"
//#include "Mesquite_QuadLagrangeShape.hpp"
//#include "Mesquite_TetLagrangeShape.hpp"
//#include "Mesquite_TriLagrangeShape.hpp"
#include "Mesquite_ElementMaxQM.hpp"
#include "Mesquite_domain.hpp"

using namespace Mesquite;

#include <vector>
#include <algorithm>

int main( int argc, char* argv[] )
{
  int two;
  CLArgs::ToggleArg freeonly;
  CLArgs::IntRangeArg histogram( &two );
  CLArgs args( "msqquality", "Assess mesh quality",
               "Caculate various quality metrics for a mesh," 
               "and optinally export a VTK file for which quality "
               "values are stored as attribute data." );
  args.int_flag( 'H', "ints", "Print histograms with specified number of intervals", &histogram );
  args.toggle_flag( 'f', "Assess quality only for elements with at least one free vertex", &freeonly );
  args.add_required_arg( "input_file" );
  args.add_optional_arg( "output_file" );
  add_domain_args( args );
  std::vector<std::string> files;
  if (!args.parse_options( argc, argv, files, std::cerr )) {
    args.print_usage( std::cerr );
    return 1;
  }
  
  MsqError err;
  MeshImpl mesh;
  mesh.read_vtk( files.front().c_str(), err );
  if (err) {
    std::cerr << err << std::endl 
                    << "Failed to read file: " << files.front() << std::endl;
    return 2;
  }
  
  MeshDomain* domain = process_domain_args( &mesh );

  QualityAssessor qa( true, freeonly.value(), "INVERTED" );
  IdealWeightInverseMeanRatio imr;
  SizeMetric size;
  IdealShapeTarget tc;
  TShapeNB1 tm;
  TQualityMetric tmp( &tc, &tm );
  ElementMaxQM max_tmp( &tmp );
  
  int intervals = histogram.seen() ? histogram.value() : 0;
  qa.add_quality_assessment( &imr, intervals, 0.0, "InverseMeanRatio" );
  qa.add_quality_assessment( &size, intervals, 0.0, "Size" );
  qa.add_quality_assessment( &max_tmp, intervals, 0.0, "TMP_Shape" );
  qa.tag_fixed_elements( "FIXED_ELEMS" );

//  QuadLagrangeShape quad;
//  TriLagrangeShape tri;
//  TetLagrangeShape tet;
  InstructionQueue q;
//  q.set_mapping_function( &quad );
//  q.set_mapping_function( &tri );
//  q.set_mapping_function( &tet );
  
  q.add_quality_assessor( &qa, err );
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, domain);
  q.run_instructions( &mesh_and_domain, err );
  delete domain;
  if (err) {
    std::cerr << err << std::endl;
    return 3;
  }
  
  if (files.size() > 1) {
    mesh.write_vtk( files[1].c_str(), err );
    if (err) {
      std::cerr << err << std::endl 
                      << "Failed to write file: " << files[1] << std::endl;
      return 2;
    }
  }
  
  return 0;
}
