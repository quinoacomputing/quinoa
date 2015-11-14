/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retian certain rights to this software.

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
 
    (2006) kraftche@cae.wisc.edu
   
  ***************************************************************** */

#include "meshfiles.h"

#include "Mesquite.hpp"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::ostream;
#include <memory>
using std::auto_ptr;
#include <ctype.h>

#include "Mesquite_MeshImpl.hpp"
#include "Mesquite_MsqError.hpp"
#include "Mesquite_InstructionQueue.hpp"
#include "Mesquite_TerminationCriterion.hpp"
#include "Mesquite_QualityAssessor.hpp"
#include "Mesquite_PMeanPTemplate.hpp"
#include "Mesquite_PatchPowerMeanP.hpp"
#include "Mesquite_ConjugateGradient.hpp"
#include "Mesquite_PlanarDomain.hpp"
#include "Mesquite_IdealShapeTarget.hpp"
#include "Mesquite_ConditionNumberQualityMetric.hpp"
#include "Mesquite_ReferenceMesh.hpp"
#include "Mesquite_RefMeshTargetCalculator.hpp"

#include "Mesquite_TQualityMetric.hpp"
#include "Mesquite_ElementPMeanP.hpp"
#include "Mesquite_VertexPMeanP.hpp"
  
#include "Mesquite_TSquared.hpp"
#include "Mesquite_TShapeNB1.hpp"
#include "Mesquite_TShapeB1.hpp"
#include "Mesquite_TShapeOrientNB1.hpp"
#include "Mesquite_TShapeOrientNB2.hpp"
#include "Mesquite_TShapeOrientB1.hpp"
#include "Mesquite_TShapeSize2DNB1.hpp"
#include "Mesquite_TShapeSizeB1.hpp"
#include "Mesquite_TShapeSize2DB2.hpp"
#include "Mesquite_TShapeSizeB3.hpp"
#include "Mesquite_TShapeSizeOrientNB1.hpp"
#include "Mesquite_TShapeSizeOrientB1.hpp"
#include "Mesquite_TShapeSizeOrientB2.hpp"

using namespace Mesquite;

static const struct { TMetric* u; const char* n; }
metrics[] = { { new TSquared,           "TSquared"               },
              { new TShapeNB1,          "Shape"                  },
              { new TShapeB1,           "ShapeBarrier"           },
              { new TShapeOrientNB1,    "ShapeOrient1"           },
              { new TShapeOrientNB2,    "ShapeOrient2"           },
              { new TShapeOrientB1,     "ShapeOrientBarrier"     },
              { new TShapeSize2DNB1,    "ShapeSize"              },
              { new TShapeSizeB1,       "ShapeSizeBarrier1"      },
              { new TShapeSize2DB2,     "ShapeSizeBarrier2"      },
              { new TShapeSizeB3,       "ShapeSizeBarrier3"      },
              { new TShapeSizeOrientB1, "ShapeSizeOrient1"       },
              { new TShapeSizeOrientB1, "ShapeSizeOrientBarrier1"},
              { new TShapeSizeOrientB2, "ShapeSizeOrientBarrier2"},
              { 0, 0 } };

enum AveragingScheme { NONE = 0, ELEMENT, VERTEX, PATCH };
const char* const averaging_names[] = { "none", "element", "vertex", "patch" };

// default values
const double DEFAULT_OF_POWER = 1.0;
const unsigned DEFAULT_METRIC_IDX = 0;
const AveragingScheme DEFAULT_AVG_SCHEME = NONE;
const char DEFAULT_INPUT_FILE[] = MESH_FILES_DIR "2D/vtk/quads/untangled/quads_4by2_bad.vtk";
const char DEFAULT_OUTPUT_FILE[] = "./out.vtk";

static PlanarDomain make_domain( Mesh* mesh, MsqError& );

static int do_smoother( const char* input_file, 
                        const char* output_file, 
                        const char* ref_mesh_file,
                        double of_power, 
                        unsigned metric_idx,
                        AveragingScheme avg_scheme )
{
  MsqPrintError err(cerr);
  
  TMetric *const target_metric = metrics[metric_idx].u;
  cout << "Input file:  " << input_file << endl;
  cout << "Metric:      ";
  if (avg_scheme != NONE)
    cout << averaging_names[avg_scheme] << " average of ";
  cout << metrics[metric_idx].n << endl;
  cout << "Of Power:    " << of_power << endl;
  
  
  auto_ptr<TargetCalculator> tc;
  auto_ptr<MeshImpl> ref_mesh_impl;
  auto_ptr<ReferenceMesh> ref_mesh;
  if (ref_mesh_file) {
    ref_mesh_impl.reset(new MeshImpl);
    ref_mesh_impl->read_vtk( ref_mesh_file, err );
    if (MSQ_CHKERR(err)) return 2;
    ref_mesh.reset( new ReferenceMesh( ref_mesh_impl.get() ));
    tc.reset( new RefMeshTargetCalculator( ref_mesh.get() ) );
  }
  else {
    tc.reset( new IdealShapeTarget( ) );
  }
    
  TQualityMetric jacobian_metric( tc.get(), target_metric );
  ElementPMeanP elem_avg( of_power, &jacobian_metric );
  VertexPMeanP vtx_avg( of_power, &jacobian_metric );
  QualityMetric* mmetrics[] = { &jacobian_metric, &elem_avg, &vtx_avg, &jacobian_metric };
  QualityMetric* metric = mmetrics[avg_scheme];

  TerminationCriterion outer, inner;
  outer.add_iteration_limit( 1 );
  inner.add_absolute_vertex_movement( 1e-4 );
  inner.add_iteration_limit( 100 );
  
  PMeanPTemplate obj1( of_power, metric );
  PatchPowerMeanP obj2( of_power, metric );
  ObjectiveFunction& objective = *((avg_scheme == PATCH) ? (ObjectiveFunction*)&obj2 : (ObjectiveFunction*)&obj1);
  
  ConjugateGradient solver( &objective, err );
  if (MSQ_CHKERR(err)) return 1;
  solver.set_inner_termination_criterion( &inner );
  solver.set_outer_termination_criterion( &outer );
  solver.use_global_patch();
  
  ConditionNumberQualityMetric qm_metric;
  QualityAssessor before_assessor;
  QualityAssessor after_assessor;
  before_assessor.add_quality_assessment( metric, 10);
  before_assessor.add_quality_assessment( &qm_metric );
  after_assessor.add_quality_assessment( metric, 10 );

  InstructionQueue q;
  q.add_quality_assessor( &before_assessor, err );
  q.set_master_quality_improver( &solver, err );
  q.add_quality_assessor( &after_assessor, err );
  
  MeshImpl mesh;
  mesh.read_vtk( input_file, err );
  if (MSQ_CHKERR(err)) return 2;
  PlanarDomain geom = make_domain( &mesh, err );
  if (MSQ_CHKERR(err)) return 1;
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, &geom);
  q.run_instructions( &mesh_and_domain, err );
  if (MSQ_CHKERR(err)) return 3;
  mesh.write_vtk( output_file, err );
  if (MSQ_CHKERR(err)) return 2;
  cout << "Wrote: " << output_file << endl;

  before_assessor.scale_histograms(&after_assessor);
  
  return 0;
}

static PlanarDomain make_domain( Mesh* mesh, MsqError& err )
{
    // calculate bounding box of mesh vertices
  Vector3D minimum(  HUGE_VAL,  HUGE_VAL,  HUGE_VAL );
  Vector3D maximum( -HUGE_VAL, -HUGE_VAL, -HUGE_VAL );
  std::vector<Mesh::VertexHandle> vertices;
  mesh->get_all_vertices( vertices, err ); 
  if (MSQ_CHKERR(err)) { return PlanarDomain( minimum, maximum ); }
  if (vertices.empty()) {
    std::cerr << "Mesh contains no vertices" << std::endl;
    exit(1);
  }
  std::vector<MsqVertex> coords(vertices.size());
  mesh->vertices_get_coordinates( arrptr(vertices), arrptr(coords), vertices.size(), err );
  if (MSQ_CHKERR(err)) { return PlanarDomain( minimum, maximum ); }
  std::vector<MsqVertex>::const_iterator i;
  for (i = coords.begin(); i != coords.end(); ++i) {
    const MsqVertex& v = *i;
    for (unsigned j = 0; j < 3; ++j) {
      if (v[j] < minimum[j])
        minimum[j] = v[j];
      if (v[j] > maximum[j])
        maximum[j] = v[j];
    }
  }
  
    // Look for a "zero" plane
  int k;
  maximum -= minimum;
  for (k = 2; k >= 0 && maximum[k] > 1e-6; --k);
  if (k < 0) {
    MSQ_SETERR(err)("Cannot determine plane of mesh.", MsqError::INVALID_STATE );
    return PlanarDomain( minimum, maximum ); 
  }
    
  Vector3D point( 0.0, 0.0, 0.0 ), normal( 0.0, 0.0, 0.0 );
  normal[k] = 1.0;
  point[k] = minimum[k];
  return PlanarDomain( normal, point );
}


static void usage( const char* argv0, bool help = false )
{
  ostream& str = help ? cout : cerr;
  str << argv0 << " [-p <power>] [-m metric_name] [-a averaging] [-e]"
               << " -r [ref_mesh] [input_file] [output_file]" << endl
      << argv0 << " <-l|-h>" << std::endl;
  if (help) {
    str << "     -p  : specify exponent value for p-mean^p OF template (default: " << DEFAULT_OF_POWER << ")" << endl
        << "     -m  : specify 2D target metric to use (default: " << metrics[DEFAULT_METRIC_IDX].n << ")" << endl
        << "     -a  : specify metric averaging scheme (default: " << averaging_names[DEFAULT_AVG_SCHEME] << ")" << endl
        << "              (none,vertex,element,patch)" << endl
        << "     -e  : sample at mid-edge points (default is corners only)" << endl
        << "     -r  : use reference mesh instead of ideal elements for targets" << endl
        << "     -l  : list available metrics" << endl
        << "     -h  : this help output" << endl
        << " default input file:  " << DEFAULT_INPUT_FILE << endl
        << " default output file: " << DEFAULT_OUTPUT_FILE << endl
        << endl;
  }
  exit (!help);
}

static void check_next_arg( int argc, char* argv[], int& i )
{
  if (i == argc) {
    cerr << "Expected argument following \"" << argv[i] << '"' << endl;
    usage( argv[0] );
  }
  ++i;
}

static double parse_double( int argc, char* argv[], int& i )
{
  check_next_arg( argc, argv, i );
  char* endptr;
  double result = strtod( argv[i], &endptr );
  if (endptr == argv[i] || *endptr) {
    cerr << "Expected real number following \"" << argv[i-1] << '"' << endl;
    usage( argv[0] );
  }
  return result;
}

static int comp_string_start( const char* p, const char* s ) {
  int i;
  for (i = 0; p[i]; ++i) 
    if (tolower(p[i]) != tolower(s[i]))
      return 0;
  return s[i] ? -1 : 1;
}

static AveragingScheme parse_averaging( int argc, char* argv[], int& i )
{
  check_next_arg( argc, argv, i );
  for (unsigned j = 0; j < 4; ++j) 
    if (comp_string_start( argv[i], averaging_names[j]))
      return (AveragingScheme)j;
  cerr << "Expected one of { ";
  for (unsigned j = 0; j < 4; ++j)
    cerr << '"' << averaging_names[j] << "\" ";
  cerr << "} following \"" << argv[i-1] << '"' << endl;
  usage( argv[0] );
  return NONE;
}

static unsigned parse_metric( int argc, char* argv[], int& i )
{
  check_next_arg( argc, argv, i );
  unsigned part = 0, all = 0, count = 0, have_all = 0;
  for (unsigned j = 0; metrics[j].u; ++j) {
    if (unsigned k = comp_string_start( argv[i], metrics[j].n)) {
      if (k > 0) {
        all = j;
        have_all = 1;
      }
      else {
        part = j;
        ++count;
      }
    }
  }
  
  if (have_all)
    return all;
  
  if (count) {
    if (count == 1)
      return part;
    cerr << "Ambiguous metric name: \"" << argv[i] << '"' << endl;
    usage( argv[0] );
  }
  
  cerr << "Invalid metric name following \"" << argv[i-1] << "\" option" << endl;
  usage( argv[0] );
  return (unsigned)-1;
}

static void list_metrics( )
{
  cout << "Available metrics:" << endl;
  for (unsigned j = 0; metrics[j].u; ++j) 
    cout << "\t" << metrics[j].n << endl;
  exit(0);
}

int main( int argc, char* argv[] )
{
  MsqPrintError err(cout);
  
    // CL settings
  double of_power = DEFAULT_OF_POWER;
  unsigned metric_idx = DEFAULT_METRIC_IDX;
  AveragingScheme avg_scheme = DEFAULT_AVG_SCHEME;
  const char* input_file = 0;
  const char* output_file = 0;
  const char* ref_mesh_file = 0;
  
  bool proc_opts = true;
  for (int i = 1; i < argc; ++i) {
    if (!proc_opts || argv[i][0] != '-') {
      if (output_file) {
        cerr << "Unexpected file name: \"" << argv[i] << '"' << endl;
        usage(argv[0]);
      }
      else if (input_file) 
        output_file = argv[i];
      else
        input_file = argv[i];
      continue;
    }
    
    if (!argv[i][1] || argv[i][2]) {
      cerr << "Invalid option: \"" << argv[i] << '"' << endl;
      usage(argv[0]);
    }
    
    switch( argv[i][1] ) {
      case 'p': of_power     = parse_double   ( argc, argv, i );    break;
      case 'm': metric_idx   = parse_metric   ( argc, argv, i );    break;
      case 'a': avg_scheme   = parse_averaging( argc, argv, i );    break;
      case 'r': check_next_arg(argc,argv,i); ref_mesh_file=argv[i]; break;
      case '-': proc_opts    = false;                               break;
      case 'h':                usage( argv[0], true );              break;
      case 'l':                list_metrics();                      break;
      default:
        cerr << "Invalid option: \"" << argv[i] << '"' << endl;
        usage(argv[0]);
    }
  }
  
  if (!input_file)
    input_file = DEFAULT_INPUT_FILE;
  if (!output_file)
    output_file = DEFAULT_OUTPUT_FILE;
  
  return do_smoother( input_file, 
                      output_file, 
                      ref_mesh_file,
                      of_power, 
                      metric_idx,
                      avg_scheme );
}
