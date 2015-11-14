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

    (2011) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file voshell.cpp
 *  \brief Implement some of the examples from N. Voshell.
 *  \author Jason Kraftcheck 
 *  \author Nick Voshell
 */

#include "Mesquite_ShapeImprover.hpp"
#include "Mesquite_UntangleWrapper.hpp"
#include "Mesquite_MeshImpl.hpp"
#include "Mesquite_MsqError.hpp"
#include "Mesquite_QualityAssessor.hpp"

#include "Mesquite_PlanarDomain.hpp"
#include "Mesquite_CylinderDomain.hpp"
#include "Mesquite_SphericalDomain.hpp"
#include "Mesquite_MeshDomain1D.hpp"
#include "Mesquite_DomainClassifier.hpp"

#include "Mesquite_HexLagrangeShape.hpp"

#include <iostream>
#include <string>
#include <cstdlib>

using namespace MESQUITE_NS;

std::string get_homogenious_example( DomainClassifier& geom, MeshImpl& mesh, MsqError& err );

std::string get_part_example_tri( DomainClassifier& geom, MeshImpl& mesh, MsqError& err );

std::string get_part_example_quad( DomainClassifier& geom, MeshImpl& mesh, MsqError& err );

std::string get_sphere_cube_example( DomainClassifier& geom, MeshImpl& mesh, MsqError& err );

std::string get_cut_cube_example( DomainClassifier& geom, MeshImpl& mesh, MsqError& err );

std::string get_sphere_cylinder_example( DomainClassifier& geom, MeshImpl& mesh, MsqError& err );

std::string get_hex_3d_part_example( DomainClassifier& geom, MeshImpl& mesh, MsqError& err );


typedef std::string (*example_setup_t)( DomainClassifier& geom, MeshImpl& mesh, MsqError& err );
struct Example {
  char flag;
  example_setup_t func;
  const char* desc;
};

int run_example( const Example& e, bool write_output_file );


const Example examples[] = {
  { 'H', &get_homogenious_example, "homogenous mesh (curvey surface) example" },
  { 't', &get_part_example_tri,    "part mesh with triangles" },
  { 'q', &get_part_example_quad,   "part mesh with quads" },
  { 'c', &get_sphere_cube_example, "cube with subtracted hemisphere" },
  { 'l', &get_cut_cube_example,    "cube with slot removed" },
  { 's', &get_sphere_cylinder_example, "cylinder with subtracted sphere" },
  { 'x', &get_hex_3d_part_example, "rectangular prism with two cylindrical slots" }
};
const int num_examples = sizeof(examples)/sizeof(examples[0]);



void usage( const char* argv0, bool help = false )
{
  std::cerr << "Usage: " << argv0 << "[-w] ";
  for (int i = 0; i < num_examples; ++i)
    std::cerr << " [-" << examples[i].flag << "]";
  std::cerr << std::endl;
  
  if (!help) {
    std::cerr << "       " << argv0 << " -h" << std::endl;
    std::exit(1);
  }
  std::cerr << "  -w : Write final meshes to VTK file" << std::endl;
  std::cerr << "NOTE: If no options are specified, all examples will be run" << std::endl;
  for (int i = 0; i < num_examples; ++i)
    std::cerr << "  -" << examples[i].flag << " : Run " << examples[i].desc << std::endl;
  std::cerr << std::endl;
  std::exit(0);
}



int main( int argc, char* argv[] )
{
  bool write_final_meshes = false;
  std::vector<Example> list;
  for (int i = 1; i < argc; ++i) {
    if (argv[i][0] != '-') {
      std::cerr << "Invalid argument: \"" << argv[i] << '"' << std::endl;
      usage(argv[0]);
    }
    for (int j = 1; argv[i][j]; ++j) {
      if (argv[i][j] == 'w') {
        write_final_meshes = true;
      }
      else {
        int idx;
        for (idx = 0; idx < num_examples; ++idx)
          if (examples[idx].flag == argv[i][j])
            break;
        if (idx == num_examples) {
          if (argv[i][j] == 'h')
            usage(argv[0],true);
          else {
            std::cerr << "Invalid flag: '" << argv[i][j] << "'" << std::endl;
            usage(argv[0]);
          }
        }
        list.push_back(examples[idx]);
      }
    }
  }
  
  const Example* exlist;
  int exlistlen;
  if (list.empty()) {
    exlist = examples;
    exlistlen = num_examples;
  }
  else {
    exlist = &list[0];
    exlistlen = list.size();
  }
  
  int result = 0;
  for (int i = 0; i < exlistlen; ++i)
    result += run_example( exlist[i], write_final_meshes );
  
  return result;
}



int run_example( const Example& e, bool write_output_file )
{
  MsqPrintError err(std::cerr);
  MeshImpl mesh;
  DomainClassifier domain;
  HexLagrangeShape hex27;
  
  std::cout << std::endl 
            << "--------------------------------------------------------------------"
            << std::endl
            << e.desc << std::endl
            << "--------------------------------------------------------------------"
            << std::endl;
   
  std::string name = e.func( domain, mesh, err );
  if (MSQ_CHKERR(err)) return 2;
  std::cout << "Loaded mesh from: " << name << std::endl;
  
  UntangleWrapper untangler;
  untangler.set_slaved_ho_node_mode( Settings::SLAVE_NONE );
  untangler.set_mapping_function( &hex27 );
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, &domain, false, true);
  untangler.run_instructions( &mesh_and_domain, err ); 
  if (MSQ_CHKERR(err)) return 1;
  ShapeImprover smoother;
  smoother.set_slaved_ho_node_mode( Settings::SLAVE_NONE );
  smoother.set_mapping_function( &hex27 );
  smoother.set_vertex_movement_limit_factor( 0.05 );
  smoother.run_instructions( &mesh_and_domain, err );
  if (MSQ_CHKERR(err)) return 1;
  
  if (write_output_file) {
    size_t idx = name.find( ".vtk" );
    if (idx != std::string::npos) {
      std::string newname( name.substr(0, idx) );
      newname += ".out";
      newname += name.substr(idx);
      name.swap(newname);
    }
    else {
      name += ".out";
    }
    
    mesh.write_vtk( name.c_str(), err ); MSQ_CHKERR(err);
    std::cout << "Write mesh to file: " << name << std::endl;
  }
  
  return smoother.quality_assessor().invalid_elements();
}

void get_planar_example( const char* filename,
                         DomainClassifier& geom,
                         MeshImpl& mesh,
                         MsqError& err )
{
  static PlanarDomain z(PlanarDomain::XY);
  mesh.read_vtk(filename, err); MSQ_ERRRTN(err);
  mesh.mark_skin_fixed(err); MSQ_ERRRTN(err);
  DomainClassifier::DomainSet set(&z);
  mesh.get_all_vertices( set.vertices, err ); MSQ_ERRRTN(err);
  mesh.get_all_elements( set.elements, err ); MSQ_ERRRTN(err);
  DomainClassifier::classify_by_handle( geom, &mesh, &set, 1, err ); MSQ_ERRRTN(err);
}

std::string get_homogenious_example( DomainClassifier& geom, MeshImpl& mesh, MsqError& err )
{
  const char filename[] = SRCDIR "homogeneousPart.vtk";
  get_planar_example( filename, geom, mesh, err ); MSQ_CHKERR(err);
  return filename;
}


std::string get_part_example_tri( DomainClassifier& geom, MeshImpl& mesh, MsqError& err )
{
  const char filename[] = SRCDIR "triPart.vtk";
  get_planar_example( filename, geom, mesh, err ); MSQ_CHKERR(err);
  return filename;
}


std::string get_part_example_quad( DomainClassifier& geom, MeshImpl& mesh, MsqError& err )
{
  const char filename[] = SRCDIR "quadPart.vtk";
  get_planar_example( filename, geom, mesh, err ); MSQ_CHKERR(err);
  return filename;
}


std::string get_sphere_cube_example( DomainClassifier& geom, MeshImpl& mesh, MsqError& err )
{
  const char filename[] = SRCDIR "sphereCube.vtk";

  const Vector3D vec_i(5,0,0), vec_ni(-5,0,0);
  const Vector3D vec_j(0,5,0), vec_nj(0,-5,0);
  const Vector3D vec_k(0,0,5), vec_nk(0,0,-5);

  const Vector3D vec_tfl(5,-5,5), vec_tfr(5,5,5);
  const Vector3D vec_tbl(-5,-5,5), vec_tbr(-5,5,5);
  const Vector3D vec_bfl(5,-5,-5), vec_bfr(5,5,-5);
  const Vector3D vec_bbl(-5,-5,-5), vec_bbr(-5,5,-5);

  const Vector3D vec_ts(5,0,2), vec_bs(5,0,-2);


  //++ 0D domains ++

  static PointDomain pt_tfl(vec_tfl), pt_tfr(vec_tfr);
  static PointDomain pt_tbl(vec_tbl), pt_tbr(vec_tbr);
  static PointDomain pt_bfl(vec_bfl), pt_bfr(vec_bfr);
  static PointDomain pt_bbl(vec_bbl), pt_bbr(vec_bbr);


  //++ 1D domains ++

  //top square
  static LineDomain ln_tf(vec_tfl,vec_j), lin_tb(vec_tbr,vec_nj);
  static LineDomain ln_tl(vec_tfl,vec_ni), lin_tr(vec_tbr,vec_i);
  //bottom square
  static LineDomain ln_bf(vec_bfl,vec_j), lin_bb(vec_bbr,vec_nj);
  static LineDomain ln_bl(vec_bfl,vec_ni), lin_br(vec_bbr,vec_i);
  //sides
  static LineDomain ln_lf(vec_bfl,vec_k), lin_rf(vec_bfr,vec_k);
  static LineDomain ln_lb(vec_bbl,vec_k), lin_rb(vec_bbr,vec_k);

  //top sphere
  //CircleDomain cr_tf(vec_ts,vec_i,1.0), cr_tn(vec_ts,vec_k,1.0);
  //bottom sphere
  //CircleDomain cr_bf(vec_bs,vec_i,1.0), cr_bn(vec_bs,vec_k,1.0);
  static CircleDomain cr_ct(vec_k,vec_k,1.0);

  //++ 2D domains ++

  //cube
  static PlanarDomain sf_i(vec_i,vec_i), sf_ni(vec_ni,vec_ni);
  static PlanarDomain sf_j(vec_j,vec_j), sf_nj(vec_nj,vec_nj);
  static PlanarDomain sf_k(vec_k,vec_k), sf_nk(vec_nk,vec_nk);
  //cut
  //CylinderDomain cy(1.0,vec_k,vec_bs);
  static SphericalDomain sp_t(vec_k,1.0);//, sp_b(vec_bs,1.0);

  MeshDomain* base_domains[] = {
    &pt_tfl, &pt_tfr,
    &pt_tbl, &pt_tbr,
    &pt_bfl, &pt_bfr,
    &pt_bbl, &pt_bbr,

    &ln_tf, &lin_tb,
    &ln_tl, &lin_tr,
    &ln_bf, &lin_bb,
    &ln_bl, &lin_br,
    &ln_lf, &lin_rf,
    &ln_lb, &lin_rb,

    //  &cr_tf, &cr_tn,
    //  &cr_bf, &cr_bn,
    &cr_ct,

    &sf_i, &sf_ni,
    &sf_j, &sf_nj,
    &sf_k, &sf_nk,
    //&cy,
    &sp_t//, &sp_b
  };
  const int NDOM = sizeof(base_domains)/sizeof(base_domains[0]);

  int dim_array[NDOM] = {
    0,0,0,0,0,0,0,0,
    1,1,1,1,1,1,1,1,
    1,1,1,1,1,//1,1,1,
    2,2,2,2,2,2,2//,2,2
  };

  //++ MESH & DOMAIN ++

  mesh.read_vtk(filename, err); MSQ_ERRZERO(err);
  DomainClassifier::classify_skin_geometrically (geom, &mesh, 0.1, base_domains, dim_array, NDOM, err);
  MSQ_ERRZERO(err);
  mesh.set_skin_flags( false, false, true, err ); MSQ_ERRZERO(err);
  
  return filename;
}


std::string get_cut_cube_example( DomainClassifier& geom, MeshImpl& mesh, MsqError& err )
{
  const char filename[] = SRCDIR "cutCube.vtk";

  const Vector3D vec_i(5,0,0), vec_ni(-5,0,0);
  const Vector3D vec_j(0,5,0), vec_nj(0,-5,0);
  const Vector3D vec_k(0,0,5), vec_nk(0,0,-5);

  const Vector3D vec_tfl(5,-5,5), vec_tfr(5,5,5);
  const Vector3D vec_tbl(-5,-5,5), vec_tbr(-5,5,5);
  const Vector3D vec_bfl(5,-5,-5), vec_bfr(5,5,-5);
  const Vector3D vec_bbl(-5,-5,-5), vec_bbr(-5,5,-5);

  const Vector3D vec_ts(5,0,2), vec_bs(5,0,-2);


  //++ 0D domains ++

  static PointDomain pt_tfl(vec_tfl), pt_tfr(vec_tfr);
  static PointDomain pt_tbl(vec_tbl), pt_tbr(vec_tbr);
  static PointDomain pt_bfl(vec_bfl), pt_bfr(vec_bfr);
  static PointDomain pt_bbl(vec_bbl), pt_bbr(vec_bbr);


  //++ 1D domains ++

  //top square
  static LineDomain ln_tf(vec_tfl,vec_j), lin_tb(vec_tbr,vec_nj);
  static LineDomain ln_tl(vec_tfl,vec_ni), lin_tr(vec_tbr,vec_i);
  //bottom square
  static LineDomain ln_bf(vec_bfl,vec_j), lin_bb(vec_bbr,vec_nj);
  static LineDomain ln_bl(vec_bfl,vec_ni), lin_br(vec_bbr,vec_i);
  //sides
  static LineDomain ln_lf(vec_bfl,vec_k), lin_rf(vec_bfr,vec_k);
  static LineDomain ln_lb(vec_bbl,vec_k), lin_rb(vec_bbr,vec_k);

  //top sphere
  static CircleDomain cr_tf(vec_ts,vec_i,1.0), cr_tn(vec_ts,vec_k,1.0);
  //bottom sphere
  static CircleDomain cr_bf(vec_bs,vec_i,1.0), cr_bn(vec_bs,vec_k,1.0);


  //++ 2D domains ++

  //cube
  static PlanarDomain sf_i(vec_i,vec_i), sf_ni(vec_ni,vec_ni);
  static PlanarDomain sf_j(vec_j,vec_j), sf_nj(vec_nj,vec_nj);
  static PlanarDomain sf_k(vec_k,vec_k), sf_nk(vec_nk,vec_nk);
  //cut
  static CylinderDomain cy(1.0,vec_k,vec_bs);
  static SphericalDomain sp_t(vec_ts,1.0), sp_b(vec_bs,1.0);


  MeshDomain* base_domains[] = {
    &pt_tfl, &pt_tfr,
    &pt_tbl, &pt_tbr,
    &pt_bfl, &pt_bfr,
    &pt_bbl, &pt_bbr,

    &ln_tf, &lin_tb,
    &ln_tl, &lin_tr,
    &ln_bf, &lin_bb,
    &ln_bl, &lin_br,
    &ln_lf, &lin_rf,
    &ln_lb, &lin_rb,
    &cr_tf, &cr_tn,
    &cr_bf, &cr_bn,

    &sf_i, &sf_ni,
    &sf_j, &sf_nj,
    &sf_k, &sf_nk,
    &cy,
    &sp_t, &sp_b
  };
  const int NDOM = sizeof(base_domains)/sizeof(base_domains[0]);

  int dim_array[NDOM] = {
    0,0,0,0,0,0,0,0,
    1,1,1,1,1,1,1,1,
    1,1,1,1,1,1,1,1,
    2,2,2,2,2,2,2,2,2
  };

  //++ MESH & DOMAIN ++

  mesh.read_vtk(filename, err); MSQ_ERRZERO(err);
  DomainClassifier::classify_skin_geometrically (geom, &mesh, 0.1, base_domains, dim_array, NDOM, err);
  MSQ_ERRZERO(err);
  mesh.set_skin_flags( false, false, true, err ); MSQ_ERRZERO(err);
  
  return filename;
}


std::string get_sphere_cylinder_example( DomainClassifier& geom, MeshImpl& mesh, MsqError& err )
{
  const char filename[] = SRCDIR "sphereCylinder_1194_inv.vtk";

  const Vector3D vec_k(0,0,8), vec_nk(0,0,-8);
  const Vector3D vec_c(0,0,5);

  //++ 0D domains ++

  //++ 1D domains ++

  //top circle
  static CircleDomain cr_to(vec_k,vec_k,8.0), cr_ti(vec_k,vec_k,4.0);
  //bottom circle
  static CircleDomain cr_bo(vec_nk,vec_nk,8.0);

  //++ 2D domains ++

  static PlanarDomain sf_t(vec_k,vec_k), sf_b(vec_nk,vec_nk);
  static CylinderDomain cy(8.0,vec_k,vec_k);
  static SphericalDomain sp(vec_c,5.0);

  MeshDomain* base_domains[] = {
    &cr_to, &cr_ti, &cr_bo,
    &sf_t, &sf_b, &cy, &sp
  };
  const int NDOM = sizeof(base_domains)/sizeof(base_domains[0]);

  int dim_array[NDOM] = {
    1, 1, 1, 
    2, 2, 2, 2
  };

  //++ MESH & DOMAIN ++

  mesh.read_vtk(filename, err); MSQ_ERRZERO(err);
  DomainClassifier::classify_skin_geometrically (geom, &mesh, 0.001, base_domains, dim_array, NDOM, err);
  MSQ_ERRZERO(err);
  mesh.set_skin_flags( false, false, true, err ); MSQ_ERRZERO(err);
  
  return filename;
}


std::string get_hex_3d_part_example( DomainClassifier& geom, MeshImpl& mesh, MsqError& err )
{
  const char filename[] = SRCDIR "hex3Dpart.vtk";

  //2D domains

  const Vector3D vec_i(1,0,0), vec_j(0,1,0), vec_k(0,0,1);
  const Vector3D vec_ni(-1,0,0), vec_nj(0,-1,0), vec_nk(0,0,-1);

  const Vector3D vec_left(-11.5,0,0), vec_right(11.5,0,0);
  const Vector3D vec_top(0,5,0), vec_bottom(0,-5,0);
  const Vector3D vec_front(0,0,5), vec_back(0,0,-5);

  //1D domains

  const Vector3D pt_top_front(0,5,5),pt_top_back(0,5,-5);
  const Vector3D pt_bottom_front(0,-5,5),pt_bottom_back(0,-5,-5);

  const Vector3D pt_top_left(-11.5,5,0),pt_top_right(11.5,5,0);
  const Vector3D pt_bottom_left(-11.5,-5,0),pt_bottom_right(11.5,-5,0);

  const Vector3D pt_left_front(-11.5,0,5),pt_left_back(-11.5,0,-5);
  const Vector3D pt_right_front(11.5,0,5),pt_right_back(11.5,0,-5);


  const Vector3D cpt_top_left(-1.5,5,0),cpt_top_right(1.5,5,0);
  const Vector3D cpt_bottom_left(-1.5,-5,0),cpt_bottom_right(1.5,-5,0);

  //0D domains

  const Vector3D pt_tlf(-11.5,5,5), pt_tlb(-11.5,5,-5);
  const Vector3D pt_trf(11.5,5,5), pt_trb(11.5,5,-5);
  const Vector3D pt_blf(-11.5,-5,5), pt_blb(-11.5,-5,-5);
  const Vector3D pt_brf(11.5,-5,5), pt_brb(11.5,-5,-5);

  const Vector3D pt_c_tlf(-1.5,5,5), pt_c_tlb(-1.5,5,-5);
  const Vector3D pt_c_trf(1.5,5,5), pt_c_trb(1.5,5,-5);
  const Vector3D pt_c_blf(-1.5,-5,5), pt_c_blb(-1.5,-5,-5);
  const Vector3D pt_c_brf(1.5,-5,5), pt_c_brb(1.5,-5,-5);

  static PointDomain p00(pt_tlf);
  static PointDomain p01(pt_tlb);
  static PointDomain p02(pt_trf);
  static PointDomain p03(pt_trb);
  static PointDomain p04(pt_blf);
  static PointDomain p05(pt_blb);
  static PointDomain p06(pt_brf);
  static PointDomain p07(pt_brb);

  static PointDomain p08(pt_c_tlf);
  static PointDomain p09(pt_c_tlb);
  static PointDomain p10(pt_c_trf);
  static PointDomain p11(pt_c_trb);
  static PointDomain p12(pt_c_blf);
  static PointDomain p13(pt_c_blb);
  static PointDomain p14(pt_c_brf);
  static PointDomain p15(pt_c_brb);

  static CircleDomain c00(pt_top_front,vec_k,1.5);
  static CircleDomain c01(pt_top_back,vec_k,1.5);
  static CircleDomain c02(pt_bottom_front,vec_k,1.5);
  static CircleDomain c03(pt_bottom_back,vec_k,1.5);

  static LineDomain l00(cpt_top_left,vec_k);
  static LineDomain l01(cpt_top_right,vec_k);
  static LineDomain l02(cpt_bottom_left,vec_k);
  static LineDomain l03(cpt_bottom_right,vec_k);

  static LineDomain l04(pt_top_front,vec_i);
  static LineDomain l05(pt_top_back,vec_i);
  static LineDomain l06(pt_bottom_front,vec_i);
  static LineDomain l07(pt_bottom_back,vec_i);
  static LineDomain l08(pt_top_left,vec_k);
  static LineDomain l09(pt_top_right,vec_k);
  static LineDomain l10(pt_bottom_left,vec_k);
  static LineDomain l11(pt_bottom_right,vec_k);
  static LineDomain l12(pt_left_front,vec_j);
  static LineDomain l13(pt_left_back,vec_j);
  static LineDomain l14(pt_right_front,vec_j);
  static LineDomain l15(pt_right_back,vec_j);

  static CylinderDomain C00(1.5,vec_k,vec_top,false);
  static CylinderDomain C01(1.5,vec_k,vec_bottom,false);
  static PlanarDomain P00(vec_ni,vec_left);
  static PlanarDomain P01(vec_i,vec_right);
  static PlanarDomain P02(vec_j,vec_top);
  static PlanarDomain P03(vec_nj,vec_bottom);
  static PlanarDomain P04(vec_k,vec_front);
  static PlanarDomain P05(vec_nk,vec_back);

  const int NDOM = 44;
  MeshDomain* base_domains[NDOM] = {
    &p00,
    &p01,
    &p02,
    &p03,
    &p04,
    &p05,
    &p06,
    &p07,

    &p08,
    &p09,
    &p10,
    &p11,
    &p12,
    &p13,
    &p14,
    &p15,

    &c00,
    &c01,
    &c02,
    &c03,

    &l00,
    &l01,
    &l02,
    &l03,

    &l04,
    &l05,
    &l06,
    &l07,
    &l08,
    &l09,
    &l10,
    &l11,
    &l12,
    &l13,
    &l14,
    &l15,

    &C00,
    &C01,
    &P00,
    &P01,
    &P02,
    &P03,
    &P04,
    &P05
  };

  int dim_array[NDOM] = { 
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    1, 1, 1, 1, 
    1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    2, 2, 2, 2, 2, 2, 2, 2
  };


  //++ MESH & DOMAIN ++

  mesh.read_vtk(filename, err); MSQ_ERRZERO(err);
  DomainClassifier::classify_skin_geometrically (geom, &mesh, 0.1, base_domains, dim_array, NDOM, err);
  MSQ_ERRZERO(err);
  mesh.set_skin_flags( false, false, true, err ); MSQ_ERRZERO(err);
  
  return filename;
}



