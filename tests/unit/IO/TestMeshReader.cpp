// *****************************************************************************
/*!
  \file      tests/unit/IO/TestMeshReader.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Unit tests for MeshReader polymorphic interface
  \details   Unit tests for MeshReader polymorphic interface
*/
// *****************************************************************************

#include <iterator>
#include <algorithm>

#include "NoWarning/tut.hpp"

#include "TUTConfig.hpp"
#include "QuinoaConfig.hpp"
#include "MeshReader.hpp"
#include "ContainerUtil.hpp"

#ifndef DOXYGEN_GENERATING_OUTPUT

namespace tut {

//! All tests in group inherited from this base
struct MeshReader_common {

  // Coordinates for simple tetrahedron-mesh
  std::vector< tk::real > box24_coord { 0,   0,   0,
                                        1,   0,   0,
                                        1,   1,   0,
                                        0,   1,   0,
                                        0,   0,   1,
                                        1,   0,   1,
                                        1,   1,   1,
                                        0,   1,   1,
                                        0.5, 0.5, 0,
                                        0.5, 0.5, 1,
                                        0.5, 0,   0.5,
                                        1,   0.5, 0.5,
                                        0.5, 1,   0.5,
                                        0,   0.5, 0.5 };

  // Element connectivity for simple tetrahedron-mesh
  std::vector< std::size_t > box24_inpoel { 11, 13,  8, 10,
                                             9, 13, 12, 11,
                                            13, 12, 11,  8,
                                             9, 13, 11, 10,
                                             0, 13,  4, 10,
                                             6,  5,  9, 11,
                                            13,  7,  4,  9,
                                             7,  6,  9, 12,
                                             6, 12,  2, 11,
                                             0,  3, 13,  8,
                                            12,  3,  2,  8,
                                             2,  1, 11,  8,
                                             3,  7, 13, 12,
                                             5,  4,  9, 10,
                                             0,  1,  8, 10,
                                             1,  5, 11, 10,
                                             5,  9, 11, 10,
                                             1, 11,  8, 10,
                                             4, 13,  9, 10,
                                            13,  7,  9, 12,
                                            12,  2, 11,  8,
                                             6,  9, 12, 11,
                                            13,  3, 12,  8,
                                            13,  0,  8, 10 };

  //! Verify MeshReader by invoking various member functions yielding tetrahedra
  //! \param[in,out] mr MeshReader object to verify
  void verifyTets( tk::MeshReader& mr ) {
    // Read mesh graph in serial
    std::vector< std::size_t > ginpoel, inpoel, triinpoel;
    std::unordered_map< std::size_t, std::size_t > lid;
    tk::UnsMesh::Coords coord;
    std::unordered_map< std::size_t, std::set< std::size_t > > elemBlockId;
    mr.readMeshPart( ginpoel, inpoel, triinpoel, lid, coord, elemBlockId );

    // Test if the number of elements is correct
    ensure_equals( "number of elements incorrect",
                   inpoel.size()/4, box24_inpoel.size()/4 );

    // Test if the mesh element connectivity is correct
    ensure( "element connectivity incorrect", inpoel == box24_inpoel );

    // Test if the number of coordinates is correct
    ensure_equals( "number of x coordinates incorrect",
                   coord[0].size(), box24_coord.size()/3 );
    ensure_equals( "number of y coordinates incorrect",
                   coord[1].size(), box24_coord.size()/3 );
    ensure_equals( "number of z coordinates incorrect",
                   coord[2].size(), box24_coord.size()/3 );

    // Test if the mesh node coordinates are correct
    std::vector< tk::real > x, y, z;
    for (std::size_t p=0; p<box24_coord.size()/3; ++p) {
      x.push_back( box24_coord[p*3] );
      y.push_back( box24_coord[p*3+1] );
      z.push_back( box24_coord[p*3+2] );
    }
    ensure( "nodes' x coordinates incorrect", coord[0] == x );
    ensure( "nodes' y coordinates incorrect", coord[1] == y );
    ensure( "nodes' z coordinates incorrect", coord[2] == z );
  }

  //! Verify MeshReader by invoking various member functions yielding triangles
  //! \param[in,out] mr MeshReader object to verify
  void verifyTris( tk::MeshReader& mr ) {
    // Read side set faces
    std::map< int, std::vector< std::size_t > > bface;
    std::map< int, std::vector< std::size_t > > faceid;
    mr.readSidesetFaces( bface, faceid );
    // Test if the number of boundary face element ids is correct
    ensure_equals( "number of boundary face elements incorrect",
                   tk::sumvalsize(bface), 2398 );
    // Test if the number of side set sides is correct
    ensure_equals( "number of side set sides incorrect",
                   tk::sumvalsize(faceid), 2398 );

    // Read face connectivity
    std::vector< std::size_t > triinpoel;
    mr.readFaces( triinpoel );
    // Test if the number of faces is correct
    ensure_equals( "number of faces in triangle connectivity incorrect",
                   triinpoel.size()/3, 2398 );

    // Read node lists associated to side sets
    auto bnode = mr.readSidesetNodes();
    // Test if the number of nodes is correct
    ensure_equals( "number of nodes of sidesets incorrect",
                   tk::sumvalsize(bnode), 1365 );
  }
};

//! Test group shortcuts
using MeshReader_group = test_group< MeshReader_common, MAX_TESTS_IN_GROUP >;
using MeshReader_object = MeshReader_group::object;

//! Define test group
//! \note Those test groups whose name contains "MPISingle" will be started as
//!    MPI tests (from a Charm++ nodegroup) and from only a single MPI rank.
static MeshReader_group MeshReader( "IO/MeshReader_MPISingle" );

//! Test definitions for group

//! Test mesh reader contructor dispatching to ExodusII mesh reader
template<> template<>
void MeshReader_object::test< 1 >() {
  set_test_name( "ctor dispatching to ExodusII reader" );

  //! Mesh reader configured for ExodusIIMesReader for a simple mesh w/o faces
  tk::MeshReader er( tk::regression_dir() + "/meshconv/gmsh_output/box_24.exo" );
  //! Mesh reader configured for ExodusIIMesReader for mesh w/ faces
  tk::MeshReader erf( tk::regression_dir() +
         "/inciter/transport/GaussHump/unitsquare_01_3.6k.exo" );

  // Verify the output of mesh reader dispatching to ExodusII reader using two
  // different meshes, one with only tets, one with faces/sidesets.
  verifyTets( er );
  verifyTris( erf );
}

//! Test mesh reader copy constructor dispatching to ExodusII mesh reader
template<> template<>
void MeshReader_object::test< 2 >() {
  set_test_name( "copy ctor dispatching to ExodusII reader" );

  //! Mesh reader configured for ExodusIIMesReader for a simple mesh w/o faces
  tk::MeshReader er( tk::regression_dir() + "/meshconv/gmsh_output/box_24.exo" );
  //! Mesh reader configured for ExodusIIMesReader for mesh w/ faces
  tk::MeshReader erf( tk::regression_dir() +
         "/inciter/transport/GaussHump/unitsquare_01_3.6k.exo" );

  std::vector< tk::MeshReader > v;

  // Invoke copy constructor
  v.push_back( er );
  // Verify that source of copy still works
  verifyTets( er );
  // Verify that the copy works
  verifyTets( v[0] );

  // Invoke copy constructor
  v.push_back( erf );
  // Verify that source of copy still works
  verifyTris( erf );
  // Verify that the copy works
  verifyTris( v[1] );
}

//! Test mesh reader move constructor dispatching to ExodusII mesh reader
template<> template<>
void MeshReader_object::test< 3 >() {
  set_test_name( "move ctor dispatching to ExodusII reader" );

  //! Mesh reader configured for ExodusIIMesReader for a simple mesh w/o faces
  tk::MeshReader er( tk::regression_dir() + "/meshconv/gmsh_output/box_24.exo" );
  //! Mesh reader configured for ExodusIIMesReader for mesh w/ faces
  tk::MeshReader erf( tk::regression_dir() +
         "/inciter/transport/GaussHump/unitsquare_01_3.6k.exo" );

  std::vector< tk::MeshReader > v;

  // Invoke move constructor
  auto p = er;
  v.emplace_back( std::move(p) );
  // Verify that the newly moved mesh reader works
  verifyTets( v[0] );

  // Invoke move constructor
  auto q = erf;
  v.emplace_back( std::move(q) );
  // Verify that the newly moved mesh reader works
  verifyTris( v[1] );
}

//! Test mesh reader copy assignment dispatching to ExodusII mesh reader
template<> template<>
void MeshReader_object::test< 4 >() {
  set_test_name( "copy assignment dispatching to ExoII reader" );

  //! Mesh reader configured for ExodusIIMesReader for a simple mesh w/o faces
  tk::MeshReader er( tk::regression_dir() + "/meshconv/gmsh_output/box_24.exo" );
  //! Mesh reader configured for ExodusIIMesReader for mesh w/ faces
  tk::MeshReader erf( tk::regression_dir() +
         "/inciter/transport/GaussHump/unitsquare_01_3.6k.exo" );

  // Invoke constructor
  tk::MeshReader q( er );
  // Invoke copy assignment
  q = er;
  // Verify that source of copy still works
  verifyTets( er );
  // Verify that the copy works
  verifyTets( q );

  // Invoke constructor
  tk::MeshReader p( erf );
  // Invoke copy assignment
  p = erf;
  // Verify that source of copy still works
  verifyTris( erf );
  // Verify that the copy works
  verifyTris( p );
}

//! Test mesh reader move assignment dispatching to ExodusII mesh reader
template<> template<>
void MeshReader_object::test< 5 >() {
  set_test_name( "move assignment dispatching to ExoII reader" );

  //! Mesh reader configured for ExodusIIMesReader for a simple mesh w/o faces
  tk::MeshReader er( tk::regression_dir() + "/meshconv/gmsh_output/box_24.exo" );
  //! Mesh reader configured for ExodusIIMesReader for mesh w/ faces
  tk::MeshReader erf( tk::regression_dir() +
         "/inciter/transport/GaussHump/unitsquare_01_3.6k.exo" );

  // Invoke move assignment
  auto c = er;
  auto p = std::move(c);
  // Verify that the newly moved mesh reader works
  verifyTets( p );

  // Invoke move assignment
  auto d = erf;
  auto q = std::move(d);
  // Verify that the newly moved mesh reader works
  verifyTris( q );
}

//! Test mesh reader contructor throws on undetectable file
template<> template<>
void MeshReader_object::test< 6 >() {
  set_test_name( "ctor throws on input file not detected" );

  try {
    // Pass a CMakeLists.txt's name as mesh file, should throw
    tk::MeshReader e( tk::regression_dir()+"/CMakeLists.txt" );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }
}

} // tut::

#endif  // DOXYGEN_GENERATING_OUTPUT
