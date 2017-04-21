// *****************************************************************************
/*!
  \file      src/UnitTest/tests/IO/TestMesh.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Unit tests for unstructured-mesh reader and writers in IO
  \details   Unit tests for unstructured-mesh reader and writers in IO
*/
// *****************************************************************************
#ifndef test_Mesh_h
#define test_Mesh_h

#include "NoWarning/tut.h"

#include "MeshFactory.h"
#include "Reorder.h"
#include "DerivedData.h"
#include "ProcessControl.h"
#include "GmshMeshWriter.h"
#include "GmshMeshReader.h"
#include "ExodusIIMeshWriter.h"
#include "ExodusIIMeshReader.h"
#include "NetgenMeshWriter.h"
#include "NetgenMeshReader.h"

namespace tut {

//! All tests in group inherited from this base
struct Mesh_common {

  //! Generic test function for testing writing and reading a tetrahedron mesh
  //! \param[in] reader Reader type
  //! \param[in] ascii Boolean selecting ASCII (TEXT) or binary mesh type
  //! \author J. Bakosi
  void testPureTetMesh( tk::MeshReader reader, bool ascii = false ) {
    // Coordinates for simple tetrahedron-mesh
    std::vector< tk::real > coord { 0,   0,   0,
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
    std::vector< std::size_t > inpoel { 12, 14,  9, 11,
                                        10, 14, 13, 12,
                                        14, 13, 12,  9,
                                        10, 14, 12, 11,
                                        1,  14,  5, 11,
                                        7,   6, 10, 12,
                                        14,  8,  5, 10,
                                        8,   7, 10, 13,
                                        7,  13,  3, 12,
                                        1,   4, 14,  9,
                                        13,  4,  3,  9,
                                        3,   2, 12,  9,
                                        4,   8, 14, 13,
                                        6,   5, 10, 11,
                                        1,   2,  9, 11,
                                        2,   6, 12, 11,
                                        6,  10, 12, 11,
                                        2,  12,  9, 11,
                                        5,  14, 10, 11,
                                        14,  8, 10, 13,
                                        13,  3, 12,  9,
                                        7,  10, 13, 12,
                                        14,  4, 13,  9,
                                        14,  1,  9, 11 };

    // Shift node IDs to start from zero
    tk::shiftToZero( inpoel );

    // Create unstructured-mesh object initializing element connectivity
    tk::UnsMesh outmesh( std::move(inpoel) );

    // Fill output mesh point coordinates
    for (std::size_t p=0; p<coord.size()/3; ++p) {
      outmesh.x().push_back( coord[p*3] );
      outmesh.y().push_back( coord[p*3+1] );
      outmesh.z().push_back( coord[p*3+2] );
    }

    std::string filename;

    // Write out mesh to file in format selected. The writer must be in its own
    // scope so the destructor closes the file immediately, otherwise we risk of
    // attempting to read from an unwritten file below. Also the filenames
    // should differ among the various readers, since the tests can be executed
    // concurrently.
    if (reader == tk::MeshReader::GMSH && ascii) {

      filename = "out_gmsh_asc.msh";
      tk::GmshMeshWriter( filename, tk::GmshFileType::ASCII ).
        writeMesh( outmesh );
 
    } else if (reader == tk::MeshReader::GMSH && !ascii) {

      filename = "out_gmsh_bin.msh";
      tk::GmshMeshWriter( filename, tk::GmshFileType::BINARY ).
        writeMesh( outmesh );

    } else if (reader == tk::MeshReader::NETGEN) {

      filename = "out.mesh";
      tk::NetgenMeshWriter( filename ).writeMesh( outmesh );

    } else if (reader == tk::MeshReader::EXODUSII) {

      filename = "out.exo";
      tk::ExodusIIMeshWriter( filename, tk::ExoWriter::CREATE ).
        writeMesh( outmesh );

    }

    // Create unstructured-mesh object to read into
    tk::UnsMesh inmesh;

    // Read in mesh just written out
    if (reader == tk::MeshReader::GMSH && ascii) {

      tk::GmshMeshReader( filename ).readMesh( inmesh );
 
    } else if (reader == tk::MeshReader::GMSH && !ascii) {

      tk::GmshMeshReader( filename ).readMesh( inmesh );

    } else if (reader == tk::MeshReader::NETGEN) {

      tk::NetgenMeshReader( filename ).readMesh( inmesh );

    } else if (reader == tk::MeshReader::EXODUSII) {

      tk::ExodusIIMeshReader( filename ).readMesh( inmesh );

    }

    // Test if mesh extents are the same as was written out
    ensure_equals( "number of nodes incorrect",
                   coord.size()/3, inmesh.x().size() );
    ensure_equals( "number of nodes incorrect",
                   coord.size()/3, inmesh.y().size() );
    ensure_equals( "number of nodes incorrect",
                   coord.size()/3, inmesh.z().size() );
    ensure_equals( "number of elements incorrect",
                   outmesh.tetinpoel().size(), inmesh.tetinpoel().size() );

    // Test if mesh nodes' coordinates are the same as was written out
    std::vector< tk::real > x, y, z;
    for (std::size_t p=0; p<coord.size()/3; ++p) {
       x.push_back( coord[p*3] );
       y.push_back( coord[p*3+1] );
       z.push_back( coord[p*3+2] );
    }
    ensure( "nodes' x coordinates incorrect", x == inmesh.x() );
    ensure( "nodes' y coordinates incorrect", y == inmesh.y() );
    ensure( "nodes' z coordinates incorrect", z == inmesh.z() );

    // Test if mesh element connectivity is the same as was written out
    ensure( "element connectivity incorrect",
            outmesh.tetinpoel() == inmesh.tetinpoel() );

    // remove mesh file from disk
    tk::rm( filename );
  }

};

//! Test group shortcuts
using Mesh_group = test_group< Mesh_common, MAX_TESTS_IN_GROUP >;
using Mesh_object = Mesh_group::object;

//! Define test group
static Mesh_group Mesh( "IO/Mesh" );

//! Test definitions for group

//! Write and read Gmsh ascii mesh
//! \author J. Bakosi
template<> template<>
void Mesh_object::test< 1 >() {
  set_test_name( "write/read Gmsh ASCII tet-mesh" );
  testPureTetMesh( tk::MeshReader::GMSH, true );
}

//! Write and read Gmsh binary mesh
//! \author J. Bakosi
template<> template<>
void Mesh_object::test< 2 >() {
  set_test_name( "write/read Gmsh binary tet-mesh" );
  testPureTetMesh( tk::MeshReader::GMSH );
}

//! Write and read ExodusII mesh
//! \author J. Bakosi
template<> template<>
void Mesh_object::test< 3 >() {
  set_test_name( "write/read ExodusII tet-mesh" );
  testPureTetMesh( tk::MeshReader::EXODUSII );
}

//! Write and read Netgen mesh
//! \author J. Bakosi
template<> template<>
void Mesh_object::test< 4 >() {
  set_test_name( "write/read Netgen tet-mesh" );
  testPureTetMesh( tk::MeshReader::NETGEN );
}

} // tut::

#endif // test_Mesh_h
