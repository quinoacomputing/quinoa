#include "Mesquite_MeshImpl.hpp"
#include "Mesquite_MsqError.hpp"
#include "Mesquite_MsqVertex.hpp"

#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <map>
#include <vector>
#include <algorithm>

const double PERTURB_FRACT = 0.4;

using namespace Mesquite;

void usage( const char* argv0 )
{
  std::cerr << "Usage: " << argv0 << " <depth> <input_file> <output_file>" << std::endl;
  exit(1);
}

// Routine to create initial mesh for test.
// o Marks vertices at a greater topological depth than the specified
//   value as slaved.  
// o Perturbs higher-order vertices on skin towards element center
// o Marks skin vertices as fixed
int main( int argc, char* argv[] )
{
  if (argc != 4)
    usage(argv[0]);
  
  char* endptr = 0;
  const long n = strtol( argv[1], &endptr, 0 );
  if (*endptr || n < 0) 
    usage(argv[0]);
  
    // read input mesh
  MeshImpl mesh;
  MsqPrintError err(std::cerr);
  mesh.read_vtk( argv[2], err );
  if (err) return 1;
  
    // get skin vertices
  mesh.mark_skin_fixed( err, true );
  if (err) return 1;
  std::vector<Mesh::VertexHandle> verts;
  mesh.get_all_vertices( verts, err );
  if (err) return 1;
  std::vector<bool> fixed;
  mesh.vertices_get_fixed_flag( arrptr(verts), fixed, verts.size(), err );
  if (err) return 1;
  std::vector<Mesh::VertexHandle> skin;
  for (size_t i = 0; i < verts.size(); ++i)
    if (fixed[i])
      skin.push_back( verts[i] );
  
    // create map for vertex depth, and initialize to 0 for skin vertices
  std::map<Mesh::VertexHandle,int> depth;
  std::map<Mesh::VertexHandle,int>::iterator d_iter;
  for (size_t i = 0; i < skin.size(); ++i)
    depth[skin[i]] = 0;
  
    // get all elements
  std::vector<Mesh::ElementHandle> curr, next;
  std::vector<Mesh::ElementHandle> conn;
  std::vector<size_t> off;
  mesh.get_all_elements( next, err );

    // build sorted list of higher-order vertices
  std::vector<Mesh::VertexHandle> higher_order;
  for (size_t i = 0; i < next.size(); ++i) {
    Mesh::ElementHandle elem = next[i];
    conn.clear();
    mesh.elements_get_attached_vertices( &elem, 1, conn, off, err );
    if (err) return 1;
    EntityTopology type;
    mesh.elements_get_topologies( &elem, &type, 1, err );
    std::copy( conn.begin() + TopologyInfo::corners(type), conn.end(), 
               std::back_inserter( higher_order ) );
  }
  std::sort( higher_order.begin(), higher_order.end() );
  higher_order.erase( std::unique( higher_order.begin(), higher_order.end() ), 
                      higher_order.end() );

    // build depth map for all vertices
  while (!next.empty()) {
    curr.swap( next );
    next.clear();
    while (!curr.empty()) {
      Mesh::ElementHandle elem = curr.back();
      curr.pop_back();
      
      conn.clear();
      mesh.elements_get_attached_vertices( &elem, 1, conn, off, err );
      if (err) return 1;
      
      int min = std::numeric_limits<int>::max();
      for (size_t i = 0; i < conn.size(); ++i) {
        d_iter = depth.find( conn[i] );
        if (d_iter != depth.end() && d_iter->second < min)
          min = d_iter->second;
      }
      
      if (min == std::numeric_limits<int>::max()) {
        next.push_back( elem );
        continue;
      }
      
      for (size_t i = 0; i < conn.size(); ++i) {
        d_iter = depth.find( conn[i] );
      
        if (d_iter == depth.end() || d_iter->second > min+1)
          depth[conn[i]] = min+1;
      }
    }
  }
  
    // write depth map to tag for debugging purposes
  std::vector<int> depth_vals(verts.size());
  for (size_t i = 0; i < verts.size(); ++i)
    depth_vals[i] = depth[verts[i]];
  TagHandle tag = mesh.tag_create( "depth", Mesh::INT, 1, 0, err );
  if (err) return 1;
  mesh.tag_set_vertex_data( tag, verts.size(), arrptr(verts), arrptr(depth_vals), err );
  if (err) return 1;
  
  
    // set tag specifying slaved vertices
  for (size_t i = 0; i < verts.size(); ++i)
    if (std::binary_search( higher_order.begin(), higher_order.end(), verts[i] ))
      depth_vals[i] = depth[verts[i]] > n;
    else
      depth_vals[i] = 0;
  tag = mesh.tag_create( "slaved", Mesh::INT, 1, 0, err );
  if (err) return 1;
  mesh.tag_set_vertex_data( tag, verts.size(), arrptr(verts), arrptr(depth_vals), err );
  if (err) return 1;
  
    // perturb mid-edge nodes along boundary
  std::vector<MsqVertex> coords;
  for (size_t i = 0; i < skin.size(); ++i) {
    if (!std::binary_search( higher_order.begin(), higher_order.end(), skin[i]))
      continue;
  
    curr.clear();
    mesh.vertices_get_attached_elements( &skin[i], 1, curr, off, err );
    if (err) return 1;
    assert(curr.size() == 1);
    conn.clear();
    mesh.elements_get_attached_vertices( arrptr(curr), 1, conn, off, err );
    if (err) return 1;
    
    // estimate element center
    coords.resize( conn.size() );
    mesh.vertices_get_coordinates( arrptr(conn), arrptr(coords), conn.size(), err );
    if (err) return 1;
    
    Vector3D mean(0.0);
    for (size_t j = 0; j < coords.size(); ++j)
      mean += coords[j];
    mean /= coords.size();
    
    size_t idx = std::find( conn.begin(), conn.end(), skin[i] ) - conn.begin();
    assert(idx < conn.size());
    Vector3D init = coords[idx];
    Vector3D pos = (1 - PERTURB_FRACT) * init + PERTURB_FRACT * mean;
    mesh.vertex_set_coordinates( skin[i], pos, err );
    if (err) return 1;
  }
  
  mesh.write_vtk( argv[3], err );
  if (err) return 1;
  
  return 0;
}

