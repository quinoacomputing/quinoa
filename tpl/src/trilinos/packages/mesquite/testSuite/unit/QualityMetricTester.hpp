/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
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
 
    (2006) kraftche@cae.wisc.edu
   
  ***************************************************************** */

#ifndef QUALITY_METRIC_TESTER_HPP
#define QUALITY_METRIC_TESTER_HPP

#include "Mesquite.hpp"
#include "Mesquite_PlanarDomain.hpp"
#include "Mesquite_Settings.hpp"
#include <algorithm>

using namespace Mesquite;

namespace MESQUITE_NS {
class PatchData;
class QualityMetric;
class ElemSampleQM;
class EdgeQM;
}

class QualityMetricTester
{
public:
  void get_ideal_tris ( PatchData& pd, bool unit_area );
  void get_ideal_quads( PatchData& pd );
  void get_ideal_hexes( PatchData& pd );
  void get_ideal_element( EntityTopology type, 
                          bool unit_area, 
                          PatchData& pd,
                          bool first_vertex_fixed = false );
  void get_ideal_element( EntityTopology type, 
                          bool unit_area, 
                          PatchData& pd,
                          int free_vertex_index );
  void get_nonideal_element  ( EntityTopology type, PatchData& pd, bool first_vertex_fixed = false );
  void get_nonideal_element  ( EntityTopology type, PatchData& pd, int free_vertex_index );
  void get_degenerate_element( EntityTopology type, PatchData& pd );
  void get_zero_element      ( EntityTopology type, PatchData& pd );
  void get_inverted_element  ( EntityTopology type, PatchData& pd );
                              
  enum ElemTypeGroup { SIMPLICIES, // triangle and tetrahedron
                       NON_MIXED_FE, // tri, quad, tet, hex
                       TWO_D,        // tri, quad, polygon
                       TWO_D_FE,     // tri, quad
                       THREE_D,      // tet, hex, pyr, wedge, septahedron, polyhedron
                       THREE_D_FE,   // tet, hex, pyr, wedge, septahedron
                       THREE_D_NON_MIXED_FE,   // tet, hex, 
                       THREE_D_FE_EXCEPT_SEPTAHEDRON, // tet, hex, pyr, wedge
                       ALL_FE_EXCEPT_SEPTAHEDRON, // tri, quad, tet, hex, pyr, wedge
                       ALL_FE, // everything except polygon and polyhedron
                       ALL }; // everything (including polyhedron)

  QualityMetricTester( ElemTypeGroup group, 
                       const Settings* set = 0 );

  QualityMetricTester( const EntityTopology* supported_elem_types,
                       size_t supported_elem_types_len,
                       const Settings* set = 0 );

    /** Ideal pyramids should be considerd to have a heigth
     *  equal to the length of a side, rather than the default which
     *  is equilateral triangle faces.
     */
  inline void ideal_pyramid_base_equals_height( bool flag ) 
    { degenHexPyramid = flag; }

    /** Test that metric evaluation succeeds for all supported element
     *  types and fails for all unsupported types */
  void test_supported_element_types( QualityMetric* qm );
  
    /** Test that metric value increases (or decreases if negate_flag is -1)
     *  as the quality of an element worsens.  Compares metric values for
     *  ideal element with unit edge length to the value for the same
     *  element with one corner vertex moved 1/2 of the distance towards
     *  the element centroid.  This test is applicable only for element-based
     *  metrics.
     */
  void test_measures_quality( QualityMetric* qm );
  
    /** Test that metric value increases (or decreases if negate_flag is -1)
     *  as the quality of the element worsen.  Compares metric values for
     *  ideal elements with unit edge length to the value for the same
     *  elements with the shared vertex moved 1/2 of the length of one of the
     *  adjacent edges.   This test is done only for TRIANGLE, QUADRILATERAL,
     *  and HEXAHEDRAL element types.
     */
  void test_measures_vertex_quality( QualityMetric* qm );
  
    /** Test measures deviation from domain */
  void test_domain_deviation_quality( QualityMetric* qm );
  
    /** Test that metric value increases (or decreases if negate_flag is -1)
     *  as the quality of an element worsens.  Compares gradient values for
     *  ideal element with unit edge length to the value for the same
     *  element with one corner vertex moved 1/2 of the distance towards
     *  the element centroid.  This test is applicable only for element-based
     *  metrics.
     */
  void test_gradient_reflects_quality( QualityMetric* qm );
  
    /** Test that metric value increases (or decreases if negate_flag is -1)
     *  as the quality of the element worsen.  Compares gradient values for
     *  ideal elements with unit edge length to the value for the same
     *  elements with the shared vertex moved 1/2 of the length of one of the
     *  adjacent edges.  This test is done only for TRIANGLE, QUADRILATERAL,
     *  and HEXAHEDRAL element types.
     */
  void test_vertex_gradient_reflects_quality( QualityMetric* qm );
  
    /** Test gradient reflects deviation from domain */
  void test_domain_deviation_gradient( QualityMetric* qm );

    /** Test evaluation of a single ideal element with unit area/volume */
  void test_evaluate_unit_element( QualityMetric* qm, EntityTopology type, double value );
    /** Test evaluation of a single ideal element with unit edge length */
  void test_evaluate_unit_edge_element( QualityMetric* qm, EntityTopology type, double value );
    /** Test evaluation of metric over patch with one free vertex surrounded 
     *  by 6 ideal unit-area tris */
  void test_evaluate_unit_tris_about_vertex( QualityMetric* qm, double expected_val );
    /** Test evaluation of metric over patch with one free vertex surrounded 
     *  by 4 unit-area quads */
  void test_evaluate_unit_quads_about_vertex( QualityMetric* qm, double expected_val );
    /** Test evaluation of metric over patch with one free vertex surrounded 
     *  by 8 unit-volume hexes */
  void test_evaluate_unit_hexes_about_vertex( QualityMetric* qm, double expected_val );
    /** Test evaluation of metric over patch with one free vertex surrounded 
     *  by 6 ideal unit-edge-length tris */
  void test_evaluate_unit_edge_tris_about_vertex( QualityMetric* qm, double expected_val );

    /** Test that evaluation of the metric for an inverted element does
     *  not pass back an error condition, and that the returned boolean
     *  feasible value is the specified value */
  void test_evaluate_inverted_element( QualityMetric* qm, bool should_succeed );

    /** Test that evaluation of the metric for degenerate element does
     *  not pass back an error condition, and that the returned boolean
     *  feasible value is the specified value */
  void test_evaluate_degenerate_element( QualityMetric* qm, bool should_succeed );

    /** Test that evaluation of the metric for zero area/volume element does
     *  not pass back an error condition, and that the returned boolean
     *  feasible value is the specified value */
  void test_evaluate_zero_element( QualityMetric* qm, bool should_succeed );

    /** Test get_evaluatinos() method for element-based metric */
  void test_get_element_evaluations(QualityMetric* qm  );
    /** Test get_evaluatinos() method for vertex-based metric */
  void test_get_vertex_evaluations( QualityMetric* qm );
    /** Test get_evaluatinos() method for sample-based metric */
  void test_get_sample_evaluations( QualityMetric* qm );
    /** Test method to get samples in an element in the EdgeQM base class */
  void test_get_edge_evaluations( EdgeQM* qm );
    /** Test method to get samples in an element in the ElemSampleQM base class */
  void test_get_in_element_evaluations( ElemSampleQM* qm );

    /** Test that evaluate_with_indices returns indices only for free vertices */
  void test_get_indices_fixed( QualityMetric* qm );

    /** Test indices from evaluate_with_indices() method */
  void test_get_element_indices( QualityMetric* qm );
    /** Test indices from evaluate_with_indices() method, assuming
     *  quality at a vertex depends only on edge-connected vertices. 
     */
  void test_get_vertex_indices( QualityMetric* qm );
  void test_get_edge_indices( EdgeQM* qm );
    /** Test indices from evaluate_with_indices() method */
  void test_get_sample_indices( QualityMetric* qm );

    /** compare results of evaluate() and evaluate_with_indices() methods */
  void compare_eval_and_eval_with_indices( QualityMetric* qm );
  void compare_eval_and_eval_with_indices( QualityMetric* qm, PatchData& pd );
    /** compare results of evaluate_with_indices() and evaluate_with_gradient() methods */
  void compare_eval_with_indices_and_eval_with_gradient( QualityMetric* qm );
  void compare_eval_with_indices_and_eval_with_gradient( QualityMetric* qm, PatchData& pd );
    /** compare results of evaluate_with_indices() and evaluate_with_Hessian() methods */
  void compare_eval_with_indices_and_eval_with_hessian( QualityMetric* qm );
  void compare_eval_with_indices_and_eval_with_hessian( QualityMetric* qm, PatchData& pd );
    /** compare results of evaluate_with_indices() and evaluate_with_Hessian_diagonal() methods */
  void compare_eval_with_indices_and_eval_with_diagonal( QualityMetric* qm );
  void compare_eval_with_indices_and_eval_with_diagonal( QualityMetric* qm, PatchData& pd );
    /** compare results of evaluate_with_gradient() and evaluate_with_Hessian() methods */
  void compare_eval_with_grad_and_eval_with_hessian( QualityMetric* qm );
  void compare_eval_with_grad_and_eval_with_hessian( QualityMetric* qm, PatchData& pd );
    /** compare results of evaluate_with_gradient() and evaluate_with_Hessian_diagonal() methods */
  void compare_eval_with_grad_and_eval_with_diagonal( QualityMetric* qm );
  void compare_eval_with_grad_and_eval_with_diagonal( QualityMetric* qm, PatchData& pd );
    /** compare results of evaluate_with_Hessian_diagonal() and evaluate_with_Hessian() methods */
  void compare_eval_with_diag_and_eval_with_hessian( QualityMetric* qm );
  void compare_eval_with_diag_and_eval_with_hessian( QualityMetric* qm, PatchData& pd );
    /** compare analytical and numerical gradient results */
  void compare_analytical_and_numerical_gradients( QualityMetric* qm );
  void compare_analytical_and_numerical_gradients( QualityMetric* qm, PatchData& pd );
    /** compare analytical and numerical Hessian results */
  void compare_analytical_and_numerical_hessians( QualityMetric* qm );
  void compare_analytical_and_numerical_hessians( QualityMetric* qm, PatchData& pd );
    /** compare analytical and numerical Hessian diagonal results */
  void compare_analytical_and_numerical_diagonals( QualityMetric* qm );
  void compare_analytical_and_numerical_diagonals( QualityMetric* qm, PatchData& pd );

    /** compare gradient w/ no fixed vertices to gradient
     *  for element with all but one vertex fixed.
     */
  void test_gradient_with_fixed_vertex( QualityMetric* qm, const Settings* settings = 0 );
  void test_gradient_with_fixed_vertex( EntityTopology type, QualityMetric* qm,
                                        const Settings* settings = 0 );
    /** compare Hessian w/ no fixed vertices to Hessian
     *  for element with all but one vertex fixed.
     */
  void test_hessian_with_fixed_vertex( QualityMetric* qm, const Settings* settings = 0 );
  void test_hessian_with_fixed_vertex( EntityTopology type, QualityMetric* qm,
                                       const Settings* settings = 0 );
    /** compare Hessian diagonal w/ no fixed vertices to Hessian
     *  for element with all but one vertex fixed.
     */
  void test_diagonal_with_fixed_vertex( QualityMetric* qm, const Settings* settings = 0 );
  void test_diagonal_with_fixed_vertex( EntityTopology type, QualityMetric* qm,
                                        const Settings* settings = 0 );

    /** Test that gradient values are zero for an ideal element.
     *  If 'unit_area' is true, then ideal elements have unit measure,
     *  otherwise they have unit edge lengths.  This test is applicable
     *  only to element-based metrics.
     */
  void test_ideal_element_zero_gradient( QualityMetric* qm, bool unit_area );
    /** Test that gradient values are zero at the shared vertex in a
     *  patch of containing ideal elements.
     *  If 'unit_area' is true, then ideal elements have unit measure,
     *  otherwise they have unit edge lengths.  Test is done only for
     *  TRIANGLE, QUADRILATERAL, and HEXAHEDRON element types.
     */
  void test_ideal_element_zero_vertex_gradient( QualityMetric* qm, bool unit_area );

    /** Test that Hessian is positive-definite for ideal elements */
  void test_ideal_element_positive_definite_Hessian( QualityMetric* qm, bool unit_area );

    /** Test that diagonal bocks of Hessian are symetrical */
  void test_symmetric_Hessian_diagonal_blocks( QualityMetric* qm );

    /** test that metric value is consistent for element translation */
  void test_location_invariant( QualityMetric* qm, bool untangler = false );
    /** test that metric value is consistent for element scalaing */
  void test_scale_invariant( QualityMetric* qm, bool untangler = false );
    /** test that metric value is consistent for element rotation */
  void test_orient_invariant( QualityMetric* qm, bool untangler = false );

    /** test that gradient values don't change with element translation */
  void test_grad_location_invariant( QualityMetric* qm, bool untangler = false );
    /** test that gradient values rotate with element rotation */
  void test_grad_orient_invariant( QualityMetric* qm, bool untangler = false );
    /** test that Hessian values don't change with element translation */
  void test_hessian_location_invariant( QualityMetric* qm, bool untangler = false );

    /** test that metric inceases (decreases) as size deviates from ideal */
  void test_measures_size( QualityMetric* qm, bool unit_area );
    /** test that metric value increases as element orientation changes from ideal */
  void test_measures_in_plane_orientation( QualityMetric* qm );
    /** test that metric value increases as element orientation changes from ideal */
  void test_measures_out_of_plane_orientation( QualityMetric* qm );

  class PatchXform { 
    public: 
      virtual ~PatchXform() {}
      virtual void xform(PatchData& pd, PlanarDomain* dom ) = 0; 
      virtual void xform_grad(std::vector<Vector3D>& grads) = 0;
  };
  
  void test_transform_invariant( QualityMetric* qm, 
                                 PatchXform& transform,
                                 bool untangler );
  
  void test_grad_transform_invariant( QualityMetric* qm, 
                                      PatchXform& transform,
                                      bool untangler );
  
  void test_hessian_transform_invariant( QualityMetric* qm, 
                                         PatchXform& transform,
                                         bool untangler );
                                 
  void test_measures_transform( QualityMetric* qm,
                                PatchXform& transform,
                                bool unit_area );
private:
  inline bool type_is_supported( EntityTopology type )
    { return std::find( types.begin(), types.end(), type ) != types.end(); }

  void test_type_is_supported( EntityTopology type, QualityMetric* qm );
  void test_type_is_not_supported( EntityTopology type, QualityMetric* qm );

  bool degenHexPyramid; //!< See: ideal_pyramid_base_equals_height()
  std::vector<EntityTopology> types;
  const Settings* mSettings;
  PlanarDomain geomPlane;
};

#endif

