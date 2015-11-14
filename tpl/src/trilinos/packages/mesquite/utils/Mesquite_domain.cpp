#include "Mesquite_MeshDomain1D.hpp"
#include "Mesquite_PlanarDomain.hpp"
#include "Mesquite_CylinderDomain.hpp"
#include "Mesquite_ConicDomain.hpp"
#include "Mesquite_SphericalDomain.hpp"
#include "Mesquite_DomainClassifier.hpp"
#include "Mesquite_MeshImpl.hpp"
#include "Mesquite_CLArgs.hpp"
#include "Mesquite_MsqError.hpp"
#include "Mesquite_TopologyInfo.hpp"

#include "Mesquite_domain.hpp"
#include <iostream>
#include <algorithm>
#include <stdlib.h>

using namespace Mesquite;

class SphereDomainArg : public CLArgs::DoubleListArgI
{
  private:
  std::vector<MeshDomain*>& domList;
  std::vector<int>& dimList;
  public:
  SphereDomainArg( std::vector<MeshDomain*>& domlist,
                   std::vector<int>& dims ) 
    : domList(domlist), dimList( dims ) {}
  virtual bool value( const std::vector<double>& list );
};
bool SphereDomainArg::value( const std::vector<double>& list )
{
  double rad = list[0];
  if (rad <= 0.0)
    return false;
  Vector3D center(0,0,0);
  if (list.size() == 4)
    center.set( list[1], list[2], list[3] );
  domList.push_back( new SphericalDomain( center, rad ) );
  dimList.push_back( 2 );
  return true;
}

class ConicDomainArg : public CLArgs::DoubleListArgI
{
  private:
  std::vector<MeshDomain*>& domList;
  std::vector<int>& dimList;
  public:
  ConicDomainArg( std::vector<MeshDomain*>& domlist,
                     std::vector<int>& dims ) 
    : domList(domlist), dimList( dims ) {}
  virtual bool value( const std::vector<double>& list );
};
bool ConicDomainArg::value( const std::vector<double>& vals )
{
  double base_rad = vals[0];
  double height = vals[1];
  Vector3D axis( 0, 0, 1 );
  if (vals.size() >= 5) {
    axis[0] = vals[2];
    axis[1] = vals[3];
    axis[2] = vals[4];
  }
  Vector3D point( 0, 0, 0 );
  if (vals.size() == 8) {
    axis[0] = vals[5];
    axis[1] = vals[6];
    axis[2] = vals[7];
  }
    
  domList.push_back( new ConicDomain( base_rad, height, axis, point ) );
  dimList.push_back( 2 );
  return true;
}

class CylinderDomainArg : public CLArgs::DoubleListArgI
{
  private:
  std::vector<MeshDomain*>& domList;
  std::vector<int>& dimList;
  public:
  CylinderDomainArg( std::vector<MeshDomain*>& domlist,
                     std::vector<int>& dims ) 
    : domList(domlist), dimList( dims ) {}
  virtual bool value( const std::vector<double>& list );
};
bool CylinderDomainArg::value( const std::vector<double>& vals )
{
  double rad = vals[0];
  Vector3D normal( vals[1], vals[2], vals[3] );
  Vector3D point(0,0,0);
  if (vals.size() == 7)
    point.set( vals[4], vals[5], vals[6] );
  domList.push_back( new CylinderDomain( rad, normal, point ) );
  dimList.push_back( 2 );
  return true;
}

class PlanarDomainArg : public CLArgs::DoubleListArgI
{
  private:
  std::vector<MeshDomain*>& domList;
  std::vector<int>& dimList;
  public:
  PlanarDomainArg( std::vector<MeshDomain*>& domlist,
                   std::vector<int>& dims ) 
    : domList(domlist), dimList( dims ) {}
  virtual bool value( const std::vector<double>& list );
};
bool PlanarDomainArg::value( const std::vector<double>& list )
{
  Vector3D normal( list[0], list[1], list[2] );
  Vector3D point(0,0,0);
  if (list.size() == 6)
    point.set( list[3], list[4], list[5] );
  domList.push_back( new PlanarDomain( normal, point ) );
  dimList.push_back( 2 );
  return true;
}

class LineDomainArg : public CLArgs::DoubleListArgI
{
  private:
  std::vector<MeshDomain*>& domList;
  std::vector<int>& dimList;
  public:
  LineDomainArg( std::vector<MeshDomain*>& domlist,
                   std::vector<int>& dims ) 
    : domList(domlist), dimList( dims ) {}
  virtual bool value( const std::vector<double>& list );
};
bool LineDomainArg::value( const std::vector<double>& vals )
{
  Vector3D dir( vals[0], vals[1], vals[2] );
  Vector3D point(0,0,0);
  if (vals.size() == 6)
    point.set( vals[3], vals[4], vals[5] );
  LineDomain* pdom = new LineDomain( point, dir );
  domList.push_back( pdom );
  dimList.push_back( 1 );
  return true;
}

class CircleDomainArg : public CLArgs::DoubleListArgI
{
  private:
  std::vector<MeshDomain*>& domList;
  std::vector<int>& dimList;
  public:
  CircleDomainArg( std::vector<MeshDomain*>& domlist,
                   std::vector<int>& dims ) 
    : domList(domlist), dimList( dims ) {}
  virtual bool value( const std::vector<double>& list );
};
bool CircleDomainArg::value( const std::vector<double>& vals )
{
  double rad = vals[0];
  Vector3D normal( vals[1], vals[2], vals[3] );
  Vector3D point(0,0,0);
  if (vals.size() == 7)
    point.set( vals[4], vals[5], vals[6] );
  CircleDomain* pdom = new CircleDomain( point, normal, rad );
  domList.push_back( pdom );
  dimList.push_back( 1 );
  return true;
}

class PointDomainArg : public CLArgs::DoubleListArgI
{
  private:
  std::vector<MeshDomain*>& domList;
  std::vector<int>& dimList;
  public:
  PointDomainArg( std::vector<MeshDomain*>& domlist,
                  std::vector<int>& dims ) 
    : domList(domlist), dimList( dims ) {}
  virtual bool value( const std::vector<double>& list );
};
bool PointDomainArg::value( const std::vector<double>& vals )
{
  Vector3D point( vals[0], vals[1], vals[2] );
  PointDomain* pdom = new PointDomain( point );
  domList.push_back( pdom );
  dimList.push_back( 0 );
  return true;
}

std::vector<MeshDomain*> domains;
std::vector<int> domain_dims;
SphereDomainArg     sphere_arg( domains, domain_dims );
ConicDomainArg       conic_arg( domains, domain_dims );
CylinderDomainArg cylinder_arg( domains, domain_dims );
PlanarDomainArg      plane_arg( domains, domain_dims );
CircleDomainArg     circle_arg( domains, domain_dims );
LineDomainArg         line_arg( domains, domain_dims );
PointDomainArg       point_arg( domains, domain_dims );
CLArgs::ToggleArg skin_mesh( false );

const char* SPHERE_VALUES[] = { "rad", "x", "y", "z" };
const char* CYLINDER_VALUES[] = { "rad", "i", "j", "k", "x", "y", "z" };
const char* CONE_VALUES[] = { "rad", "h", "i", "j", "k", "x", "y", "z" };

void add_domain_args( CLArgs& args )
{
  args.toggle_flag( SKIN_FLAG, "Mark boundary vertices as fixed (default if no domain specified)", &skin_mesh );
  args.double_list_flag( SPHERE_FLAG, "Spherical domain as center and radius", &sphere_arg );
  args.limit_list_flag( SPHERE_FLAG, 4, SPHERE_VALUES );
  args.limit_list_flag( SPHERE_FLAG, 1, SPHERE_VALUES );
  args.double_list_flag( PLANE_FLAG, "Planar domain as normal and point", &plane_arg );
  args.limit_list_flag( PLANE_FLAG, 3, CYLINDER_VALUES+1 );
  args.limit_list_flag( PLANE_FLAG, 6, CYLINDER_VALUES+1 );
  args.double_list_flag( CYLINDER_FLAG, "Cylindrical radius, axis, and point", &cylinder_arg );
  args.limit_list_flag( CYLINDER_FLAG, 4, CYLINDER_VALUES );
  args.limit_list_flag( CYLINDER_FLAG, 7, CYLINDER_VALUES );
  args.double_list_flag( CONE_FLAG, "Conic domain as base radius, height, axis, and base center", &conic_arg );
  args.limit_list_flag( CONE_FLAG, 2, CONE_VALUES );
  args.limit_list_flag( CONE_FLAG, 5, CONE_VALUES );
  args.limit_list_flag( CONE_FLAG, 8, CONE_VALUES );
  args.double_list_flag( LINE_FLAG, "Linear domain as direction and point", &line_arg );
  args.limit_list_flag( LINE_FLAG, 3, CYLINDER_VALUES+1 );
  args.limit_list_flag( LINE_FLAG, 6, CYLINDER_VALUES+1 );
  args.double_list_flag( CIRCLE_FLAG, "Circular domain as radius, normal, and center", &circle_arg );
  args.limit_list_flag( CIRCLE_FLAG, 4, CYLINDER_VALUES );
  args.limit_list_flag( CIRCLE_FLAG, 7, CYLINDER_VALUES );
  args.double_list_flag( POINT_FLAG, "Point domain", &point_arg );
  args.limit_list_flag( POINT_FLAG, 3, SPHERE_VALUES+1 );
}

MeshDomain* process_domain_args( MeshImpl* mesh )
{
  MsqPrintError err(std::cerr);
  MeshDomain* rval = 0;
  
  if (!domains.empty()) {
    int max_domain_dim = *std::max_element( domain_dims.begin(), domain_dims.end() );
    std::vector<Mesh::ElementHandle> elems;
    mesh->get_all_elements( elems, err );
    std::vector<EntityTopology> types( elems.size() );
    mesh->elements_get_topologies( arrptr(elems), arrptr(types), elems.size(), err );
    EntityTopology max_type = *std::max_element( types.begin(), types.end() );
    int max_elem_dim = TopologyInfo::dimension( max_type );
  
    if (max_domain_dim == max_elem_dim && domains.size() == 1) {
      rval = domains.front();
    }
    else {
      DomainClassifier* result = new DomainClassifier();
      if (max_domain_dim < max_elem_dim) {
        DomainClassifier::classify_skin_geometrically( *result,
                                                mesh,
                                                1e-4,
                                                arrptr(domains),
                                                arrptr(domain_dims),
                                                domains.size(),
                                                err );
      }
      else {
        DomainClassifier::classify_geometrically( *result,
                                                mesh,
                                                1e-4,
                                                arrptr(domains),
                                                arrptr(domain_dims),
                                                domains.size(),
                                                err );
      }
      rval = result;
    }
  }
  if (skin_mesh.value()) {
    mesh->mark_skin_fixed( err, false );
  }
  if (err) {
    std::cerr << err << std::endl;
    exit( 3 );
  }
  
  return rval;
}
