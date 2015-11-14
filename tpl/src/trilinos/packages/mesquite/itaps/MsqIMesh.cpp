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


/** \file MsqIMesh.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "MsqIMesh.hpp"
#include "MsqError.hpp"
#include "MeshInterface.hpp"
#include "MsqDebug.hpp"
#include "MsqVertex.hpp"
#include <assert.h>
#include "MsqIBase.hpp"
#include <algorithm>

#ifdef IMESH_MAJOR_VERSION
# define IMESH_VERSION_ATLEAST(MAJOR,MINOR) \
           1000*IMESH_MAJOR_VERSION+IMESH_MINOR_VERSION <= \
           1000*MAJOR+MINOR
#else
# define IMESH_VERSION_ATLEAST(MAJOR,MINOR) 0
#endif

namespace MESQUITE_NS {


/*************************************************************************
 *                          Mesh Definition
 ************************************************************************/

  MsqIMesh::MsqIMesh( iMesh_Instance mesh, 
                      iBase_EntitySetHandle meshset, 
		      iBase_EntityType type,
                      MsqError& err,
		      const iBase_TagHandle* fixed_tag,
		      const iBase_TagHandle* slaved_tag)
  : meshInstance(mesh), 
    inputSetType( iBase_ALL_TYPES ),
    inputSet(0),
    byteTag(0), 
    createdByteTag(false),
    geometricDimension(0)
  {
    init_active_mesh( mesh, err, fixed_tag, slaved_tag ); 
    MSQ_ERRRTN(err);  
    set_active_set( meshset, type, err );
    MSQ_ERRRTN(err);  
  }

  MsqIMesh::MsqIMesh( iMesh_Instance mesh, 
		      iBase_EntityType type,
                      MsqError& err,
		      const iBase_TagHandle* fixed_tag,
		      const iBase_TagHandle* slaved_tag)
  : meshInstance(mesh), 
    inputSetType( iBase_ALL_TYPES ),
    inputSet(0),
    byteTag(0), 
    createdByteTag(false),
    geometricDimension(0)
  {
    init_active_mesh( mesh, err, fixed_tag, slaved_tag ); 
    MSQ_ERRRTN(err);  
    
    int ierr;
    iBase_EntitySetHandle root_set;
    iMesh_getRootSet( mesh, &root_set, &ierr );
    if (ierr != iBase_SUCCESS) {
      MSQ_SETERR(err)("Invalid iMesh instance.", MsqError::INVALID_STATE );
      return;
    }
    set_active_set( root_set, type, err );
    MSQ_ERRRTN(err);  
    
    
  }

  MsqIMesh::~MsqIMesh() 
  {
    int ierr;
    if (createdByteTag)
      iMesh_destroyTag( meshInstance, byteTag, true, &ierr );
  }

  iMesh_Instance MsqIMesh::get_imesh_instance() const
  {
    return meshInstance;
  }

  iBase_EntitySetHandle MsqIMesh::get_entity_set() const
  {
    return inputSet;
  }

iBase_TagValueType MsqIMesh::check_valid_flag_tag( iBase_TagHandle tag,
                                                   const char* which_flag,
                                                   MsqError& err )
{
  int ierr, size, type;
  const int MAX_NAME_LEN = 127;
  char name[MAX_NAME_LEN+1];

  std::fill( name, name+sizeof(name), '\0' );
  iMesh_getTagName( meshInstance, tag, name, &ierr, MAX_NAME_LEN );
  name[MAX_NAME_LEN-1] = '\0'; // make sure strings are null-terminated
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)(MsqError::INVALID_ARG,
                    "Invalid tag handle for vertex %s flag",
                    which_flag );
    return iBase_ENTITY_HANDLE;
  }
  iMesh_getTagSizeBytes( meshInstance, tag, &size, &ierr );
  if (iBase_SUCCESS != ierr || size != sizeof(int)) {
    MSQ_SETERR(err)( MsqError::INVALID_STATE,
		     "Tag \"%s\" exists with invalid size", 
		     name );
    return iBase_ENTITY_HANDLE;
  }

  iMesh_getTagType( meshInstance, tag, &type, &ierr );
  if (iBase_SUCCESS != ierr || (type != iBase_INTEGER && type != iBase_BYTES)) {
    MSQ_SETERR(err)( MsqError::INVALID_STATE,
		     "Tag \"%s\" exists with invalid type", 
		     name );
    return iBase_ENTITY_HANDLE;
  }
  return static_cast<iBase_TagValueType>(type);
}

void MsqIMesh::init_active_mesh( iMesh_Instance mesh, 
                                 MsqError& err,
				 const iBase_TagHandle* fixed_tag,
				 const iBase_TagHandle* slaved_tag )
{
  int ierr;

  // Initialize topology map 
  
  const size_t mapsize = sizeof(topologyMap) / sizeof(Mesquite::EntityTopology);
  if (mapsize < iMesh_ALL_TOPOLOGIES)
  {
    MSQ_SETERR(err)("MsqIMesh needs to be updated for new iMesh element topologies.",
		    MsqError::INTERNAL_ERROR);
  }
  
  for (size_t i = 0; i <= iMesh_ALL_TOPOLOGIES; ++i)
    topologyMap[i] = Mesquite::MIXED;
  
  topologyMap[iMesh_TRIANGLE     ] = Mesquite::TRIANGLE;
  topologyMap[iMesh_QUADRILATERAL] = Mesquite::QUADRILATERAL;
  topologyMap[iMesh_TETRAHEDRON  ] = Mesquite::TETRAHEDRON;
  topologyMap[iMesh_HEXAHEDRON   ] = Mesquite::HEXAHEDRON;
  topologyMap[iMesh_PRISM        ] = Mesquite::PRISM;
  topologyMap[iMesh_PYRAMID      ] = Mesquite::PYRAMID;
  
      // Check that fixed tag is valid
  haveFixedTag = false;
  if (fixed_tag) {
    fixedTagType = check_valid_flag_tag( *fixed_tag, "fixed", err );
    MSQ_ERRRTN(err);
    haveFixedTag = true;
    fixedTag = *fixed_tag;
  }
  
      // Check that slaved tag is valid
  haveSlavedTag = false;
  if (slaved_tag) {
    slavedTagType = check_valid_flag_tag( *slaved_tag, "slaved", err );
    MSQ_ERRRTN(err);
    haveSlavedTag = true;
    slavedTag = *slaved_tag;
  }
  
    // Get/create tag for vertex byte
  iMesh_getTagHandle( meshInstance, 
                      VERTEX_BYTE_TAG_NAME,
                      &byteTag, &ierr,
                      strlen(VERTEX_BYTE_TAG_NAME) );
  if (iBase_SUCCESS != ierr) {
    iMesh_createTag( meshInstance, 
                     VERTEX_BYTE_TAG_NAME,
                     1, iBase_INTEGER,
                     &byteTag, &ierr,
                     strlen(VERTEX_BYTE_TAG_NAME) );
    if (iBase_SUCCESS != ierr) {
      MSQ_SETERR(err)( MsqError::INVALID_STATE, 
		       "Tag \"%s\" could not be created", 
		       VERTEX_BYTE_TAG_NAME );
      return;
    }
    createdByteTag = true;
  }
  else {
    int size, type;
    iMesh_getTagSizeBytes( meshInstance, byteTag, &size, &ierr );
    if (iBase_SUCCESS != ierr || size != sizeof(int)) {
      MSQ_SETERR(err)( MsqError::INVALID_STATE,
		       "Tag \"%s\" exists with invalid size", 
		       VERTEX_BYTE_TAG_NAME );
      return;
    }
    iMesh_getTagType( meshInstance, byteTag, &type, &ierr );
    if (iBase_SUCCESS != ierr || type != iBase_INTEGER) {
      MSQ_SETERR(err)( MsqError::INVALID_STATE,
		       "Tag \"%s\" exists with invalid type", 
		       VERTEX_BYTE_TAG_NAME );
      return;
    }
  }
  iMesh_getGeometricDimension( meshInstance, &geometricDimension, &ierr );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return;
  }
}


void MsqIMesh::set_fixed_tag( iBase_TagHandle tag, MsqError& err )
{
  iBase_TagValueType t = check_valid_flag_tag( tag, "fixed", err );
  MSQ_ERRRTN(err);
  fixedTag = tag;
  fixedTagType = t;
  haveFixedTag = true;
}

void MsqIMesh::clear_fixed_tag()
{
  haveFixedTag = false;
}

const iBase_TagHandle* MsqIMesh::get_fixed_tag() const
{
  return haveFixedTag ? &fixedTag : 0;
}

void MsqIMesh::set_slaved_tag( iBase_TagHandle tag, MsqError& err )
{
  iBase_TagValueType t = check_valid_flag_tag( tag, "slaved", err );
  MSQ_ERRRTN(err);
  slavedTag = tag;
  slavedTagType = t;
  haveSlavedTag = true;
}

void MsqIMesh::clear_slaved_tag()
{
  haveSlavedTag = false;
}

const iBase_TagHandle* MsqIMesh::get_slaved_tag() const
{
  return haveSlavedTag ? &slavedTag : 0;
}



void MsqIMesh::set_active_set( iBase_EntitySetHandle elem_set, 
                               iBase_EntityType type_in,
                               MsqError& err )
{
  inputSetType = type_in;
  inputSet = elem_set;
  
    // clear vertex byte
  std::vector<VertexHandle> verts;
  get_all_vertices( verts, err ); MSQ_ERRRTN(err);
  if (!verts.empty()) {
    std::vector<unsigned char> zeros( verts.size(), 0 );
    vertices_set_byte( arrptr(verts), arrptr(zeros), verts.size(), err );
    MSQ_CHKERR(err);
  }
}

  

// Returns whether this mesh lies in a 2D or 3D coordinate system.
int MsqIMesh::get_geometric_dimension(Mesquite::MsqError &err)
{
  return geometricDimension;
}
    

//************ Vertex Properties ********************

void MsqIMesh::get_flag_data( iBase_TagHandle tag,
                              bool have_tag,
                              iBase_TagValueType type,
                              const VertexHandle vert_array[],
                              std::vector<bool>& flag_array,
                              size_t num_vtx, 
                              MsqError& err )
{
  if (!num_vtx)
    return;

  if (!have_tag) {
    flag_array.clear();
    flag_array.resize( num_vtx, false );
    return;
  }

  flag_array.resize( num_vtx );

  assert( sizeof(VertexHandle) == sizeof(iBase_EntityHandle) );
  const iBase_EntityHandle* arr = reinterpret_cast<const iBase_EntityHandle*>(vert_array);

  int ierr, alloc = num_vtx, size = 0;
  assert((size_t)alloc == num_vtx); // size_t can hold larger values than int if 64-bit

  if (type == iBase_INTEGER) {
    std::vector<int> values(num_vtx);
    int* ptr = arrptr(values);
    iMesh_getIntArrData( meshInstance, arr, num_vtx, tag, &ptr, &alloc, &size, &ierr );
    for (int i = 0; i < size; ++i)
      flag_array[i] = !!values[i];
  }
  else if (type == iBase_BYTES) {
    std::vector<char> values(num_vtx);
    void* ptr = arrptr(values);
    iMesh_getArrData( meshInstance, arr, num_vtx, tag, &ptr, &alloc, &size, &ierr );
    for (int i = 0; i < size; ++i)
      flag_array[i] = !!values[i];
  }
  else {
    MSQ_SETERR(err)("Invalid tag type for vertex flag data", MsqError::INVALID_STATE);
    return ;
  }
  
    // check if query for tag data failed
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return;
  }
  
    // make sure we got back the requested number of values
  assert( static_cast<size_t>(size) == num_vtx );
}

// Returns true or false, indicating whether the vertex
// is allowed to be repositioned.  True indicates that the vertex
// is fixed and cannot be moved.  Note that this is a read-only
// property; this flag can't be modified by users of the
// Mesquite::Mesh interface.
void MsqIMesh::vertices_get_fixed_flag(
  const VertexHandle vert_array[], 
  std::vector<bool>& bool_array,
  size_t num_vtx, MsqError &err)
{
  get_flag_data( fixedTag, haveFixedTag, fixedTagType, vert_array, bool_array, num_vtx, err );
}


void MsqIMesh::vertices_get_slaved_flag(
  const VertexHandle vert_array[], 
  std::vector<bool>& bool_array,
  size_t num_vtx, MsqError &err)
{
  get_flag_data( slavedTag, haveSlavedTag, slavedTagType, vert_array, bool_array, num_vtx, err );
}

// Get vertex coordinates 
void MsqIMesh::vertices_get_coordinates(
  const Mesquite::Mesh::VertexHandle vert_array[],
  MsqVertex* coordinates, 
  size_t num_vtx, 
  MsqError &err)
{
  if (!num_vtx)
    return;

  std::vector<double> dbl_store( 3*num_vtx );
  double* dbl_array = arrptr(dbl_store);
  
  int ierr, junk = 3*num_vtx, junk2;
  assert( sizeof(VertexHandle) == sizeof(iBase_EntityHandle) );
  const iBase_EntityHandle* arr = reinterpret_cast<const iBase_EntityHandle*>(vert_array);
  iMesh_getVtxArrCoords( meshInstance, arr, num_vtx, iBase_INTERLEAVED, &dbl_array, &junk, &junk2, &ierr );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return;
  }
  
  if (geometricDimension == 2)
  {
    double* iter = dbl_array;
    for (size_t i = 0; i < num_vtx; ++i)
    {
      coordinates[i].x(*iter); ++iter;
      coordinates[i].y(*iter); ++iter;
      coordinates[i].z(0);
    }
  }
  else 
  {
    double* iter = dbl_array;
    for (size_t i = 0; i < num_vtx; ++i)
    {
      coordinates[i].set(iter);
      iter += 3;
    }
  }
}

void MsqIMesh::vertex_set_coordinates(
  Mesquite::Mesh::VertexHandle vertex,
  const Vector3D &coords, MsqError &err)
{
  int ierr;
  iBase_EntityHandle bh = static_cast<iBase_EntityHandle>(vertex);
  iMesh_setVtxCoord( meshInstance, bh, coords[0], coords[1], coords[2], &ierr );
  if (iBase_SUCCESS != ierr) 
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
}

// Each vertex has a byte-sized flag that can be used to store
// flags.  This byte's value is neither set nor used by the mesh
// implementation.  It is intended to be used by Mesquite algorithms.
// Until a vertex's byte has been explicitly set, its value is 0.
void MsqIMesh::vertex_set_byte (
  Mesquite::Mesh::VertexHandle vertex,
  unsigned char byte, MsqError &err)
{
  int ierr, value = byte;
  iBase_EntityHandle bh = static_cast<iBase_EntityHandle>(vertex);
  iMesh_setIntData( meshInstance, bh, byteTag, value, &ierr );
  if (iBase_SUCCESS != ierr) 
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
}

void MsqIMesh::vertices_set_byte (
  const VertexHandle *vert_array,
  const unsigned char *byte_array,
  size_t array_size, MsqError &err)
{
  if (!array_size)
    return;

  std::vector<int> data(array_size);
  std::copy( byte_array, byte_array + array_size, data.begin() );
  int ierr;
  assert( sizeof(VertexHandle) == sizeof(iBase_EntityHandle) );
  const iBase_EntityHandle* arr = reinterpret_cast<const iBase_EntityHandle*>(vert_array);
  iMesh_setIntArrData( meshInstance, arr, array_size, byteTag, arrptr(data), array_size, &ierr );
  if (iBase_SUCCESS != ierr) 
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
}

// Retrieve the byte value for the specified vertex or vertices.
// The byte value is 0 if it has not yet been set via one of the
// *_set_byte() functions.
void MsqIMesh::vertex_get_byte(
  Mesquite::Mesh::VertexHandle vertex,
  unsigned char *byte, MsqError &err)
{
  int ierr, value;
  iBase_EntityHandle bh = static_cast<iBase_EntityHandle>(vertex);
  iMesh_getIntData( meshInstance, bh, byteTag, &value, &ierr );
  if (iBase_SUCCESS != ierr) 
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
  *byte = value;
}

void MsqIMesh::vertices_get_byte(
  const VertexHandle *vert_array,
  unsigned char *byte_array,
  size_t array_size, MsqError &err)
{
  if (!array_size)
    return;

  std::vector<int> data(array_size);
  int ierr;
  int* ptr = arrptr(data);
  int junk1 = data.size(), junk2;
  assert( sizeof(VertexHandle) == sizeof(iBase_EntityHandle) );
  const iBase_EntityHandle* arr = reinterpret_cast<const iBase_EntityHandle*>(vert_array);
  iMesh_getIntArrData( meshInstance, arr, array_size, byteTag, &ptr, &junk1, &junk2, &ierr );
  if (iBase_SUCCESS != ierr) 
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
  std::copy( data.begin(), data.end(), byte_array );
}


//**************** Topology *****************

void MsqIMesh::get_adjacent_entities( const iBase_EntityHandle* source,
                                          size_t num_source,
                                          iBase_EntityType target_type,
                                          std::vector<EntityHandle>& target,
                                          std::vector<size_t>& offsets,
                                          MsqError& err )
{
  if (num_source == 0) {
    target.clear();
    offsets.clear();
    offsets.reserve(1);
    offsets.push_back(0);
    return;
  }
  
  int ierr, num_adj = 0, num_offset;
  
  assert( sizeof(size_t) >= sizeof(int) );
  offsets.resize( num_source + 1 );
  int* ptr2 = (int*)arrptr(offsets);
  bool expand = false;
  if (sizeof(size_t) > sizeof(int))
    expand = true;
  
  assert( sizeof(iBase_EntityHandle) == sizeof(EntityHandle) );
  bool have_adj = false;
    // If passed vector has allocated storage, try to use existing space
  if (target.capacity() >= num_source)
  {
    target.resize( target.capacity() );
    int junk1 = target.capacity(), junk3 = offsets.size();
    iBase_EntityHandle* ptr = reinterpret_cast<iBase_EntityHandle*>(arrptr(target));
    iMesh_getEntArrAdj( meshInstance, source, num_source,
                        target_type, 
                        &ptr, &junk1, &num_adj, 
                        &ptr2, &junk3, &num_offset, 
                        &ierr );
    if (iBase_SUCCESS == ierr) {
      have_adj = true;
      target.resize( num_adj );
    }
  }
  
    // If implementation passed back a size, try that
  if (!have_adj && num_adj && (unsigned)num_adj > target.capacity())
  {
    target.resize( num_adj );
    int junk1 = target.capacity(), junk3 = offsets.size();
    iBase_EntityHandle* ptr = reinterpret_cast<iBase_EntityHandle*>(arrptr(target));
    iMesh_getEntArrAdj( meshInstance, source, num_source,
                        target_type, 
                        &ptr, &junk1, &num_adj, 
                        &ptr2, &junk3, &num_offset, 
                        &ierr );
    if (iBase_SUCCESS == ierr)
      have_adj = true;
  }

    // Try with empty sidl array, and copy into elements vector
  if (!have_adj)
  {
    iBase_EntityHandle* mArray = 0;
    int junk1 = 0, junk3 = offsets.size();
    iMesh_getEntArrAdj( meshInstance, source, num_source,
                        target_type, 
                        &mArray, &junk1, &num_adj, 
                        &ptr2, &junk3, &num_offset, 
                        &ierr );
    if (iBase_SUCCESS != ierr) {
      MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
      return;
    }
    
    target.resize( num_adj );
    std::copy( mArray, mArray + num_adj, reinterpret_cast<iBase_EntityHandle*>(arrptr(target)) );
    free( mArray );
  }
  
  if (expand) {
    for (size_t i = num_offset; i > 0; --i)
      offsets[i-1] = ptr2[i-1];
  }
  
  // iMesh implementations seem to be inconsistent with regard to 
  // placing the last value on this list.
  if (offsets.size() - num_offset == 1)
    offsets[num_offset++] = num_adj;
  assert( (unsigned)num_offset == offsets.size() );
}


void MsqIMesh::vertices_get_attached_elements( 
                                     const VertexHandle* vertices,
                                     size_t num_vertex,
                                     std::vector<ElementHandle>& elements,
                                     std::vector<size_t>& offsets,
                                     MsqError& err )
{
  int ierr, cont;
  assert( sizeof(EntityHandle) == sizeof(iBase_EntityHandle) );
  const iBase_EntityHandle* verts = reinterpret_cast<const iBase_EntityHandle*>(vertices);
  get_adjacent_entities( verts, num_vertex, inputSetType, elements, offsets, err ); 
  MSQ_ERRRTN(err);
  
  iBase_EntitySetHandle root_set;
  iMesh_getRootSet( meshInstance, &root_set, &ierr );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return;
  }
  
    // Remove all elements not in inputSet
  if (root_set != inputSet) {
    std::vector<size_t>::iterator offset_iter = offsets.begin();
    size_t read_idx, write_idx;
    for (read_idx = write_idx = 0; read_idx < elements.size(); ++read_idx)
    {
      if (*offset_iter == read_idx)
      {
        *offset_iter = write_idx;
        ++offset_iter;
      }

      iBase_EntityHandle bh = static_cast<iBase_EntityHandle>(elements[read_idx]);
      iMesh_isEntContained( meshInstance, inputSet, bh, &cont, &ierr );
      if (iBase_SUCCESS != ierr) {
        MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
        return;
      }

      if (cont)
        elements[write_idx++] = elements[read_idx];
    }
    *offset_iter = write_idx;
    elements.resize(write_idx);    
  }
}


//**************** Element Topology *****************


/** Get connectivity
 *\param elements - Array of length num_elems containing elements
 *                  handles of elements for which connectivity is to
 *                  be queried.
 *\param vertices - Array of vertex handles in connectivity list.
 *\param offsets  - Indices into \ref vertex_handles, one per element
 */
void MsqIMesh::elements_get_attached_vertices(
  const ElementHandle *elements,
  size_t num_elems,
  std::vector<VertexHandle>& vertices,
  std::vector<size_t>& offsets,
  Mesquite::MsqError &err)
{
  assert( sizeof(iBase_EntityHandle) == sizeof(EntityHandle) );
  const iBase_EntityHandle* elems = reinterpret_cast<const iBase_EntityHandle*>(elements);
  get_adjacent_entities( elems, num_elems, iBase_VERTEX, vertices, offsets, err );
  MSQ_CHKERR(err);
}


void MsqIMesh::get_all_elements( std::vector<ElementHandle>& elements,
                                     MsqError& err )
{
  int ierr, count_in, count_out;

  if (inputSetType == iBase_ALL_TYPES) {
    int num_vol, num_face;

    iMesh_getNumOfType( meshInstance, inputSet, iBase_FACE, &num_face, &ierr );
    if (iBase_SUCCESS != ierr) {
      MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
      return;
    }
    iMesh_getNumOfType( meshInstance, inputSet, iBase_REGION, &num_vol, &ierr );
    if (iBase_SUCCESS != ierr) {
      MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
      return;
    }
    elements.resize( num_face + num_vol );
    if (elements.empty())
      return;
    
    iBase_EntityHandle* ptr = reinterpret_cast<iBase_EntityHandle*>(arrptr(elements));
    if (num_face) {
      count_in = num_face+num_vol;
      iMesh_getEntities( meshInstance, inputSet, 
                         iBase_FACE, iMesh_ALL_TOPOLOGIES,
                         &ptr, &count_in, &count_out, &ierr );
      if (iBase_SUCCESS != ierr) {
        MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
        return;
      }
      assert (count_out == num_face);
    }
    
    if (num_vol) {
      ptr += num_face;
      count_in = num_vol;
      iMesh_getEntities( meshInstance, inputSet, 
                         iBase_REGION, iMesh_ALL_TOPOLOGIES,
                         &ptr, &count_in, &count_out, &ierr );
      if (iBase_SUCCESS != ierr) {
        MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
        return;
      }
      assert (count_out == num_vol);
    }
  }
  else {
    int count;
    iMesh_getNumOfType( meshInstance, inputSet, inputSetType, &count, &ierr );
    if (iBase_SUCCESS != ierr) {
      MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
      return;
    }
    
    if (!count)
      return;
    elements.resize( count );
    
    iBase_EntityHandle* ptr = reinterpret_cast<iBase_EntityHandle*>(arrptr(elements));
    count_in = count;
    iMesh_getEntities( meshInstance, inputSet, 
                       inputSetType, iMesh_ALL_TOPOLOGIES,
                       &ptr, &count_in, &count_out, &ierr );
    if (iBase_SUCCESS != ierr) {
      MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
      return;
    }
    assert (count_out == count);
  }
}

void MsqIMesh::get_all_vertices( std::vector<VertexHandle>& vertices,
                                     MsqError& err )
{
  std::vector<ElementHandle> elems;
  get_all_elements( elems, err ); MSQ_CHKERR(err);
  if (elems.empty())
    return;  
  
  std::vector<size_t> offsets;
  elements_get_attached_vertices( arrptr(elems), elems.size(), vertices, offsets, err );
  MSQ_CHKERR(err);
  
  std::sort( vertices.begin(), vertices.end() );
  vertices.erase( std::unique( vertices.begin(), vertices.end() ), vertices.end() );
}
      

// Returns the topologies of the given entities.  The "entity_topologies"
// array must be at least "num_elements" in size.
void MsqIMesh::elements_get_topologies(
  const ElementHandle *element_handle_array,
  EntityTopology *element_topologies,
  size_t num_elements, MsqError &err)
{
  if (!num_elements)
    return;

    // don't copy unless we have to
  std::vector<int> topo_store;
  int* topo_array;
  if (sizeof(EntityTopology) == sizeof(int))
    topo_array = (int*)element_topologies;
  else {
    topo_store.resize(num_elements);
    topo_array = arrptr(topo_store);
  }
  
  int ierr, junk1 = num_elements, junk2;
  assert( sizeof(ElementHandle) == sizeof(iBase_EntityHandle) );
  const iBase_EntityHandle* arr = reinterpret_cast<const iBase_EntityHandle*>(element_handle_array);
  iMesh_getEntArrTopo( meshInstance, arr, num_elements, &topo_array, &junk1, &junk2, &ierr );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return;
  }
  
  for (size_t i = 0; i < num_elements; ++i)
    element_topologies[i] = topologyMap[topo_array[i]];
}

//**************** Memory Management ****************
// Tells the mesh that the client is finished with a given
// entity handle.  
void MsqIMesh::release_entity_handles(
  const Mesquite::Mesh::EntityHandle */*handle_array*/,
  size_t /*num_handles*/, MsqError &/*err*/)
{
    // Do nothing
}

// Instead of deleting a Mesh when you think you are done,
// call release().  In simple cases, the implementation could
// just call the destructor.  More sophisticated implementations
// may want to keep the Mesh object to live longer than Mesquite
// is using it.
void MsqIMesh::release()
{
}

//**************** Tags ****************
TagHandle MsqIMesh::tag_create( const std::string& name, 
                                    TagType type, unsigned length,
                                    const void* ,
                                    MsqError& err )
{
  int itaps_type;
  switch (type) {
    case Mesquite::Mesh::BYTE:   itaps_type = iBase_BYTES;         break;
    case Mesquite::Mesh::INT:    itaps_type = iBase_INTEGER;       break;
    case Mesquite::Mesh::DOUBLE: itaps_type = iBase_DOUBLE;        break;
    case Mesquite::Mesh::HANDLE: itaps_type = iBase_ENTITY_HANDLE; break;
    default:
      MSQ_SETERR(err)("Invalid tag type", MsqError::INVALID_ARG );
      return 0;
  }
  
  int ierr;
  iBase_TagHandle result;
  iMesh_createTag( meshInstance, 
                   name.c_str(), 
                   length, itaps_type,
                   &result, &ierr,
                   name.size() );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return 0;
  }
  
  return static_cast<TagHandle>(result);
}

void MsqIMesh::tag_delete( TagHandle handle, MsqError& err )
{
  int ierr;
  iMesh_destroyTag( meshInstance, static_cast<iBase_TagHandle>(handle), true, &ierr );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return;
  }
}

TagHandle MsqIMesh::tag_get( const std::string& name, MsqError& err )
{
  iBase_TagHandle handle = 0;
  int ierr;
  iMesh_getTagHandle( meshInstance, name.c_str(), &handle, &ierr, name.length() );
  if (iBase_TAG_NOT_FOUND == ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::TAG_NOT_FOUND );
    return 0;
  }
  else if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return 0;
  }
  return static_cast<TagHandle>(handle);
}


void MsqIMesh::tag_properties( TagHandle handle,
                                   std::string& name_out,
                                   TagType& type_out,
                                   unsigned& length_out,
                                   MsqError& err )
{
  char buffer[256];
  int ierr1, ierr2, ierr3, itype;
  
  iBase_TagHandle th = static_cast<iBase_TagHandle>(handle);
  iMesh_getTagName( meshInstance, th, buffer, &ierr1, sizeof(buffer) );
  iMesh_getTagSizeValues( meshInstance, th, (int*)&length_out, &ierr2 );
  iMesh_getTagType( meshInstance, th, &itype, &ierr3 );
  
  int ierr = iBase_SUCCESS;
  if (ierr1 != iBase_SUCCESS)
    ierr = ierr1;
  else if (ierr2 != iBase_SUCCESS)
    ierr = ierr2;
  else if (ierr3 != iBase_SUCCESS)
    ierr = ierr3;
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return;
  }
  
  buffer[255] = '\0';
  name_out = buffer;
  switch (itype) {
    case iBase_BYTES        : type_out = Mesquite::Mesh::BYTE  ; break;
    case iBase_INTEGER      : type_out = Mesquite::Mesh::INT   ; break;
    case iBase_DOUBLE       : type_out = Mesquite::Mesh::DOUBLE; break;
    case iBase_ENTITY_HANDLE: type_out = Mesquite::Mesh::HANDLE; break;
    default:
      MSQ_SETERR(err)("Unsupported iMesh tag type", MsqError::NOT_IMPLEMENTED );
      return;
  }
}

void MsqIMesh::tag_set_element_data( TagHandle tag, 
                                         size_t num_elems,
                                         const ElementHandle* array,
                                         const void* data,
                                         MsqError& err )
{
  tag_set_data( tag, num_elems, array, data, err );
}

void MsqIMesh::tag_set_vertex_data( TagHandle tag, 
                                        size_t num_elems,
                                        const VertexHandle* array,
                                        const void* data,
                                        MsqError& err )
{
  tag_set_data( tag, num_elems, array, data, err );
}
    
void MsqIMesh::tag_set_data( TagHandle tag, 
                                 size_t num_elems,
                                 const EntityHandle* array,
                                 const void* data,
                                 MsqError& err )
{
  int ierr, size;
  iMesh_getTagSizeBytes( meshInstance, static_cast<iBase_TagHandle>(tag), &size, &ierr );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return;
  }
  assert( sizeof(EntityHandle) == sizeof(iBase_EntityHandle) );
  const iBase_EntityHandle* arr = reinterpret_cast<const iBase_EntityHandle*>(array);
  iMesh_setArrData( meshInstance, arr, num_elems, 
                    static_cast<iBase_TagHandle>(tag), 
                    static_cast<const char*>(data), 
                    size*num_elems, &ierr );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return;
  }
}  


void MsqIMesh::tag_get_element_data( TagHandle tag, 
                                         size_t num_elems,
                                         const ElementHandle* array,
                                         void* data,
                                         MsqError& err )
{
  tag_get_data( tag, num_elems, array, data, err );
}

void MsqIMesh::tag_get_vertex_data( TagHandle tag, 
                                        size_t num_elems,
                                        const VertexHandle* array,
                                        void* data,
                                        MsqError& err )
{
  tag_get_data( tag, num_elems, array, data, err );
}
    
void MsqIMesh::tag_get_data( TagHandle tag, 
                                 size_t num_elems,
                                 const EntityHandle* array,
                                 void* data,
                                 MsqError& err )
{
  int ierr, size;
  iMesh_getTagSizeBytes( meshInstance, static_cast<iBase_TagHandle>(tag), &size, &ierr );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return;
  }
#if IMESH_VERSION_ATLEAST(1,1)
  void* ptr = data;
#else
  char* ptr = static_cast<char*>(data);
#endif
  int junk1 = size*num_elems, junk2;
  assert( sizeof(EntityHandle) == sizeof(iBase_EntityHandle) );
  const iBase_EntityHandle* arr = reinterpret_cast<const iBase_EntityHandle*>(array);
  iMesh_getArrData( meshInstance, arr, num_elems, 
                    static_cast<iBase_TagHandle>(tag), 
                    &ptr, &junk1, &junk2, &ierr );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)( process_itaps_error( ierr ), MsqError::INTERNAL_ERROR );
    return;
  }
}

} // namespace Mesquite
