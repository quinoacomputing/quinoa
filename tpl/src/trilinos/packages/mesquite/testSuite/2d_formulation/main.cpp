/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2007 Sandia National Laboratories.  Developed at the
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

    (2007) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file main.cpp
 *  \brief Try examples from "Formulation of a Target-Matrix Paradigm
 *         for Mesh Optimization", Patrick Knupp.
 *  \author Jason Kraftcheck 
 */
 
#define USE_GLOBAL_PATCH

#include "Mesquite.hpp"

#include "Mesquite_PMeanPTemplate.hpp"
#include "Mesquite_AffineMapMetric.hpp"
#include "Mesquite_ConjugateGradient.hpp"
#include "Mesquite_TerminationCriterion.hpp"
#include "Mesquite_ElementPMeanP.hpp"
#include "Mesquite_MsqError.hpp"
#include "Mesquite_TSquared.hpp"
#include "Mesquite_MeshImpl.hpp"
#include "Mesquite_PlanarDomain.hpp"
#include "Mesquite_InstructionQueue.hpp"
#include "Mesquite_TargetCalculator.hpp"
#include "Mesquite_MetricWeight.hpp"
#include "Mesquite_InverseMetricWeight.hpp"
#include "Mesquite_TargetWriter.hpp"
#include "Mesquite_WeightReader.hpp"

#include <iostream>
#include <stdlib.h>

using namespace Mesquite;
using namespace std;

const double epsilon = 2e-2;
const bool write_results = true;

#define CHKERR(A) if (A) { cerr << (A) << endl; exit(1); }

enum Grouping { SAMPLE, ELEMENT, QUADRANT, HALF };
enum Weight { UNIT, METRIC, INV_METRIC };


class IdentityTarget : public TargetCalculator {
public:
 bool get_3D_target( PatchData& , 
                     size_t ,
                     Sample ,
                     MsqMatrix<3,3>& W_out,
                     MsqError&  )
  { W_out = MsqMatrix<3,3>(1.0); return true; }

 bool get_surface_target( PatchData& , 
                          size_t ,
                          Sample ,
                          MsqMatrix<3,2>& W_out,
                          MsqError&  )
  { W_out = MsqMatrix<3,2>(1.0); return true; }

 bool get_2D_target( PatchData& , 
                     size_t ,
                     Sample ,
                     MsqMatrix<2,2>& W_out,
                     MsqError&  )
  { W_out = MsqMatrix<2,2>(1.0); return true; }

 bool have_surface_orient() const
  { return false; }
 
};

void run_test( Grouping grouping, int of_power, Weight w, const string filename )
{
  MsqError err;
  
  IdentityTarget target;
  TSquared target_metric;
  AffineMapMetric qual_metric( &target, &target_metric );
  ElementPMeanP elem_metric( of_power, &qual_metric );
  QualityMetric* qm_ptr = (grouping == ELEMENT) ? (QualityMetric*)&elem_metric : (QualityMetric*)&qual_metric;

  PMeanPTemplate OF( of_power, qm_ptr );
  ConjugateGradient solver( &OF );
  TerminationCriterion tc;
  TerminationCriterion itc;
  tc.add_absolute_vertex_movement( 1e-4 );
  itc.add_iteration_limit( 2 );
#ifdef USE_GLOBAL_PATCH
  solver.use_global_patch();
  solver.set_inner_termination_criterion( &tc );
#else
  solver.use_element_on_vertex_patch();
  solver.set_inner_termination_criterion( &itc );
  solver.set_outer_termination_criterion( &tc );
#endif
  
  MeshImpl mesh, expected_mesh;
  mesh.read_vtk( SRCDIR "/initial.vtk", err ); CHKERR(err)
//  expected_mesh.read_vtk( (filename + ".vtk").c_str(), err ); CHKERR(err)
  
  PlanarDomain plane( PlanarDomain::XY );

  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, &plane);

  MetricWeight mw( &qual_metric );
  InverseMetricWeight imw( &qual_metric );
  WeightReader reader;
  if (w == METRIC) {
    TargetWriter writer( 0, &mw );
    InstructionQueue tq;
    tq.add_target_calculator( &writer, err );
    tq.run_instructions( &mesh_and_domain, err ); CHKERR(err);
    qual_metric.set_weight_calculator( &reader );
  }
  else if (w == INV_METRIC) {
    TargetWriter writer( 0, &imw );
    InstructionQueue tq;
    tq.add_target_calculator( &writer, err );
    tq.run_instructions( &mesh_and_domain, err ); CHKERR(err);
    qual_metric.set_weight_calculator( &reader );
  }

  InstructionQueue q;
  q.set_master_quality_improver( &solver, err );
  q.run_instructions( &mesh_and_domain, err ); CHKERR(err)
/*  
  vector<Mesh::VertexHandle> vemain.cpprts;
  vector<MsqVertex> mesh_coords, expected_coords;
  mesh.get_all_vertices( verts, err ); CHKERR(err)
  mesh_coords.resize(verts.size());
  mesh.vertices_get_coordinates( arrptr(verts), arrptr(mesh_coords), verts.size(), err ); CHKERR(err)
  expected_mesh.get_all_vertices( verts, err ); CHKERR(err)
  expected_coords.resize(verts.size());
  expected_mesh.vertices_get_coordinates( arrptr(verts), arrptr(expected_coords), verts.size(), err ); CHKERR(err)
  if (expected_coords.size() != mesh_coords.size()) {
    cerr << "Invlid expected mesh.  Vertex count doesn't match" << endl;
    exit(1);
  }
  
  unsigned error_count = 0;
  for (size_t i = 0; i < mesh_coords.size(); ++i) 
    if ((expected_coords[i] - mesh_coords[i]).length_squared() > epsilon*epsilon)
      ++error_count;

  if (!error_count) 
    cout << filename << " : SUCCESS" << endl;
  else
    cout << filename << " : FAILURE (" << error_count 
         << " vertices differ by more than " << epsilon << ")" << endl;
*/
  if (write_results)
    mesh.write_vtk( (filename + ".results.vtk").c_str(), err ); CHKERR(err)
}

int main()
{
  run_test( SAMPLE, 1, UNIT, "1-1" );
  run_test( SAMPLE, 2, UNIT, "1-2" );
  run_test( SAMPLE, 4, UNIT, "1-4" );
  run_test( SAMPLE, 8, UNIT, "1-8" );
  
  run_test(  SAMPLE, 1, UNIT, "2-NW" );
  run_test( ELEMENT, 1, UNIT, "2-NE" );
  
  run_test( SAMPLE, 1,       UNIT, "3-Left"  );
  run_test( SAMPLE, 1,     METRIC, "3-Mid"   );
  run_test( SAMPLE, 1, INV_METRIC, "3-Right" );
}
