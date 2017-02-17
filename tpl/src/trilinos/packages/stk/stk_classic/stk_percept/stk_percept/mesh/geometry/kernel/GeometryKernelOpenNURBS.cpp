#include <iostream>
#include <typeinfo>

#include "GeometryKernelOpenNURBS.hpp"
#include <string>

#define DEBUG_GEOM_ON 0

GeometryKernelOpenNURBS::GeometryKernelOpenNURBS()
{
  ON::Begin();
}

GeometryKernelOpenNURBS::~GeometryKernelOpenNURBS()
{
  ON::End();
}

bool GeometryKernelOpenNURBS::read_file
(
  const std::string& file_name,
  std::vector<GeometryHandle>& geometry_entities
)
{
  FILE* archive_fp = ON::OpenFile( file_name.c_str(), "rb");
  if ( !archive_fp )
  {
    return false;
  }

  ON_BinaryFile archive( ON::read3dm, archive_fp );

  // read the contents of the file into "model"
  bool rc = onModel.Read( archive );

  // close the file
  ON::CloseFile( archive_fp );

  if ( !onModel.IsValid() )
    return false;

#if DEBUG_GEOM_ON
  std::cout << "tmp geom onModel.m_object_table.Count() = " << onModel.m_object_table.Count() << std::endl;
  for (int i=0; i<onModel.m_object_table.Count(); i++)
  {
    std::cout << "tmp geom object[" << i << "] = "
              << (onModel.m_object_table[i].m_object ? typeid(*onModel.m_object_table[i].m_object).name() : "null") 
              << " name = " << (onModel.m_object_table[i].m_object ?  get_attribute(i) : " no name ")
              << std::endl;
  }
#endif

  for (int i=0; i<onModel.m_object_table.Count(); i++)
  {
    if ( is_curve(i) )
      geometry_entities.push_back(i);
  }
  for (int i=0; i<onModel.m_object_table.Count(); i++)
  {
    if ( is_surface(i) )
      geometry_entities.push_back(i);
  }

  return rc;
}

std::string GeometryKernelOpenNURBS::get_attribute(GeometryHandle geom)
{
  ON_wString key = "geometry";
  ON_wString geometry_name;

  geometry_name = onModel.m_object_table[geom].m_attributes.m_name;

  std::string geom_name;
  geom_name.assign(geometry_name.Array(),
                   geometry_name.Array()+geometry_name.Length());

  return geom_name;
}

void GeometryKernelOpenNURBS::snap_to
(
  KernelPoint& point,
  GeometryHandle geom,
  double *converged_tolerance,
  double *uvw_computed,
  double *uvw_hint
)
{
  const ON_Surface* surface = dynamic_cast<const ON_Surface*>(onModel.m_object_table[geom].m_object);
  const ON_Curve* curve = dynamic_cast<const ON_Curve*>(onModel.m_object_table[geom].m_object);
  if (surface)
  {
    ON_3dPoint p(point);
    double u, v;

    surface->GetClosestPoint(p, &u, &v, 0.0, NULL, NULL, converged_tolerance, uvw_computed);
    if(uvw_computed)
    {
      uvw_computed[0] = u;
      uvw_computed[1] = v;
    }
    surface->EvPoint(u, v, p);
    point[0] = p.x;
    point[1] = p.y;
    point[2] = p.z;
  }
  else if (curve)
  {
    ON_3dPoint p(point);
    double u;

    curve->GetClosestPoint(p, &u);
    curve->EvPoint(u, p);
    point[0] = p.x;
    point[1] = p.y;
    point[2] = p.z;
  }
}

void GeometryKernelOpenNURBS::normal_at(KernelPoint& point, GeometryHandle geom, std::vector<double>& normal)
{
  const ON_Surface* surface = dynamic_cast<const ON_Surface*>(onModel.m_object_table[geom].m_object);
  const ON_Curve* curve = dynamic_cast<const ON_Curve*>(onModel.m_object_table[geom].m_object);
  if (surface)
  {
    ON_3dPoint p(point);
    double u, v;

    surface->GetClosestPoint(p, &u, &v);
    //surface->EvPoint(u, v, p);
    ON_3dVector norm;
    surface->EvNormal(u, v, norm);
    normal[0] = norm[0];
    normal[1] = norm[1];
    normal[2] = norm[2];
  }
  else if (curve)
  {
    ON_3dPoint p(point);
    double u;

    curve->GetClosestPoint(p, &u);
    //curve->EvPoint(u, p);
    ON_3dVector tangent;
    ON_3dVector kappa;
    //curve->EvTangent(u, p, tangent);
    curve->EvCurvature(u, p, tangent, kappa);

    // FIXME
    normal[0] = kappa[0];
    normal[1] = kappa[1];
    normal[2] = kappa[2];
  }
}

bool GeometryKernelOpenNURBS::is_curve
(
  GeometryHandle geom
) const
{
  const ON_Curve* curve = dynamic_cast<const ON_Curve*>(onModel.m_object_table[geom].m_object);
  return curve ? true : false;
}

bool GeometryKernelOpenNURBS::is_surface
(
  GeometryHandle geom
) const
{
  const ON_Surface* surface = dynamic_cast<const ON_Surface*>(onModel.m_object_table[geom].m_object);
  return surface ? true : false;
}

