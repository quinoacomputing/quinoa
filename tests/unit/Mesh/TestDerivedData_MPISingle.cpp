// *****************************************************************************
/*!
  \file      tests/unit/Mesh/TestDerivedData_MPISingle.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Unit tests for Mesh/DerivedData that require to be invoked from a
             single thread
  \details   Unit tests for Mesh/DerivedData that require to be invoked from a
             single thread.
*/
// *****************************************************************************

#include "NoWarning/tut.hpp"

#include "TUTConfig.hpp"
#include "DerivedData.hpp"
#include "Reorder.hpp"
#include "ExodusIIMeshReader.hpp"

#ifndef DOXYGEN_GENERATING_OUTPUT

namespace tut {

#if defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wsuggest-attribute=noreturn"
#endif

//! All tests in group inherited from this base
struct DerivedData_MPISingle_common {};

// Test group shortcuts
// The 2nd template argument is the max number of tests in this group. If
// omitted, the default is 50, specified in tut/tut.hpp.
using DerivedData_MPISingle_group =
  test_group< DerivedData_MPISingle_common, MAX_TESTS_IN_GROUP >;
using DerivedData_MPISingle_object = DerivedData_MPISingle_group::object;

//! Define test group
//! \note Those test groups whose name contains "MPISingle" will be started as
//!    MPI tests (from a Charm++ nodegroup) and from only a single MPI rank.
static DerivedData_MPISingle_group
  DerivedData_MPISingle( "Mesh/DerivedData_MPISingle" );

//! Test definitions for group

//! Generate and test number of boundary-faces for tet-only mesh
template<> template<>
void DerivedData_MPISingle_object::test< 1 >() {
  set_test_name( "Boundary faces for tetrahedra" );

  std::map< int, std::vector< std::size_t > > bface, t_bface;

  // Create unstructured-mesh object to read into
  tk::UnsMesh inmesh;
  // Read in mesh from file
  std::string infile( tk::regression_dir() +
                      "/meshconv/gmsh_output/box_24_ss1.exo" );
  tk::ExodusIIMeshReader er( infile );
  std::map< int, std::vector< std::size_t > > faceid;
  er.readSidesetFaces( t_bface, faceid );

  auto tnbfac = tk::sumvalsize( t_bface );

  // Test if the number of boundary faces is correct
  ensure_equals( "total number of boundary faces incorrect",
                 tnbfac, 24 );

  // Mesh connectivity for simple tetrahedron-only mesh
  std::vector< std::size_t > inpoel { 12, 14,  9, 11,
                                      10, 14, 13, 12,
                                      14, 13, 12,  9,
                                      10, 14, 12, 11,
                                       1, 14,  5, 11,
                                       7,  6, 10, 12,
                                      14,  8,  5, 10,
                                       8,  7, 10, 13,
                                       7, 13,  3, 12,
                                       1,  4, 14,  9,
                                      13,  4,  3,  9,
                                       3,  2, 12,  9,
                                       4,  8, 14, 13,
                                       6,  5, 10, 11,
                                       1,  2,  9, 11,
                                       2,  6, 12, 11,
                                       6, 10, 12, 11,
                                       2, 12,  9, 11,
                                       5, 14, 10, 11,
                                      14,  8, 10, 13,
                                      13,  3, 12,  9,
                                       7, 10, 13, 12,
                                      14,  4, 13,  9,
                                      14,  1,  9, 11 };

  // Boundary-face node connectivity for the entire mesh
  std::vector< std::size_t > t_triinpoel {  2,  9,  3,
                                            1,  9,  2,
                                            3,  9,  4,
                                            1,  4,  9,
                                            7, 10,  6,
                                            6, 10,  5,
                                            8, 10,  7,
                                           10,  8,  5,
                                            6, 11,  2,
                                            2, 11,  1,
                                           11,  6,  5,
                                           11,  5,  1,
                                            3, 12,  2,
                                           12,  6,  2,
                                            7, 12,  3,
                                           12,  7,  6,
                                           13,  7,  3,
                                            4, 13,  3,
                                           13,  8,  7,
                                            8, 13,  4,
                                            5,  8, 14,
                                            1,  5, 14,
                                            4, 14,  8,
                                            1, 14,  4 };

  // Correct Boundary-face node connectivity
  std::vector< std::size_t > correct_triinpoel {  7, 10,  6,
                                                  6, 10,  5,
                                                  8, 10,  7,
                                                 10,  8,  5,
                                                  5,  8, 14,
                                                  1,  5, 14,
                                                  4, 14,  8,
                                                  1, 14,  4,
                                                  2,  9,  3,
                                                  1,  9,  2,
                                                  3,  9,  4,
                                                  1,  4,  9,
                                                  3, 12,  2,
                                                 12,  6,  2,
                                                  7, 12,  3,
                                                 12,  7,  6,
                                                 13,  7,  3,
                                                  4, 13,  3,
                                                 13,  8,  7,
                                                  8, 13,  4,
                                                  6, 11,  2,
                                                  2, 11,  1,
                                                 11,  6,  5,
                                                 11,  5,  1 };

  std::unordered_map< std::size_t, std::size_t > lid {
          { {0}, {0} },
          { {1}, {1} },
          { {2}, {2} },
          { {3}, {3} },
          { {4}, {4} },
          { {5}, {5} },
          { {6}, {6} },
          { {7}, {7} },
          { {8}, {8} },
          { {9}, {9} },
          { {10}, {10} },
          { {11}, {11} },
          { {12}, {12} },
          { {13}, {13} } };

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );
  tk::shiftToZero( t_triinpoel );
  tk::shiftToZero( correct_triinpoel );

  std::vector< std::size_t > triinpoel;

  auto nbfac = tk::genNbfacTet( tnbfac, inpoel, t_triinpoel, t_bface, lid,
                                triinpoel, bface );

  ensure_equals( "number of boundary-faces is incorrect",
                 nbfac, tnbfac );

  ensure_equals( "total number of entries in triinpoel is incorrect",
                 triinpoel.size(), correct_triinpoel.size() );

  for(std::size_t i=0 ; i<triinpoel.size(); ++i)
  {
    ensure_equals("incorrect entry " + std::to_string(i) + " in triinpoel",
                    triinpoel[i], correct_triinpoel[i]);
  }
}

//! Generate and test boundary-element vector for tet-only mesh
template<> template<>
void DerivedData_MPISingle_object::test< 2 >() {
  set_test_name( "genBelemTet for tetrahedra" );

  std::map< int, std::vector< std::size_t > > bface;

  // Create unstructured-mesh object to read into
  tk::UnsMesh inmesh;
  // Read in mesh from file
  std::string infile( tk::regression_dir() +
                      "/meshconv/gmsh_output/box_24_ss1.exo" );
  tk::ExodusIIMeshReader er( infile );
  std::map< int, std::vector< std::size_t > > faceid;
  er.readSidesetFaces( bface, faceid );

  auto nbfac = tk::sumvalsize(bface);

  // Test if the number of boundary faces is correct
  ensure_equals( "total number of boundary faces incorrect",
                 nbfac, 24 );

  // Mesh connectivity for simple tetrahedron-only mesh
  std::vector< std::size_t > inpoel { 12, 14,  9, 11,
                                      10, 14, 13, 12,
                                      14, 13, 12,  9,
                                      10, 14, 12, 11,
                                       1, 14,  5, 11,
                                       7,  6, 10, 12,
                                      14,  8,  5, 10,
                                       8,  7, 10, 13,
                                       7, 13,  3, 12,
                                       1,  4, 14,  9,
                                      13,  4,  3,  9,
                                       3,  2, 12,  9,
                                       4,  8, 14, 13,
                                       6,  5, 10, 11,
                                       1,  2,  9, 11,
                                       2,  6, 12, 11,
                                       6, 10, 12, 11,
                                       2, 12,  9, 11,
                                       5, 14, 10, 11,
                                      14,  8, 10, 13,
                                      13,  3, 12,  9,
                                       7, 10, 13, 12,
                                      14,  4, 13,  9,
                                      14,  1,  9, 11 };

  // Boundary-face node connectivity
  std::vector< std::size_t > triinpoel {  2,  9,  3,
                                          1,  9,  2,
                                          3,  9,  4,
                                          1,  4,  9,
                                          7, 10,  6,
                                          6, 10,  5,
                                          8, 10,  7,
                                         10,  8,  5,
                                          6, 11,  2,
                                          2, 11,  1,
                                         11,  6,  5,
                                         11,  5,  1,
                                          3, 12,  2,
                                         12,  6,  2,
                                          7, 12,  3,
                                         12,  7,  6,
                                         13,  7,  3,
                                          4, 13,  3,
                                         13,  8,  7,
                                          8, 13,  4,
                                          5,  8, 14,
                                          1,  5, 14,
                                          4, 14,  8,
                                          1, 14,  4 };

  // Correct boundary elements
  std::vector< std::size_t > correct_belem{ 12,
                                            15,
                                            11,
                                            10,
                                             6,
                                            14,
                                             8,
                                             7,
                                            16,
                                            15,
                                            14,
                                             5,
                                            12,
                                            16,
                                             9,
                                             6,
                                             9,
                                            11,
                                             8,
                                            13,
                                             7,
                                             5,
                                            13,
                                            10 };

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );
  tk::shiftToZero( triinpoel );

  auto esup = tk::genEsup( inpoel, 4 );
  auto esuel = tk::genEsuelTet( inpoel, esup );
  auto ntfac = tk::genNipfac( 4, nbfac, esuel );
  auto inpofa = tk::genInpofaTet( ntfac, nbfac, inpoel, triinpoel, esuel );
  auto belem = tk::genBelemTet( nbfac, inpofa, esup );

  ensure_equals( "total number of entries in belem is incorrect",
                 belem.size(), correct_belem.size() );

  for(std::size_t i=0 ; i<belem.size(); ++i)
  {
    ensure_equals("incorrect entry " + std::to_string(i) + " in belem",
                    belem[i], correct_belem[i]-1);
  }
}

//! Generate and test boundary-element vector for tet-only mesh
template<> template<>
void DerivedData_MPISingle_object::test< 3 >() {
  set_test_name( "genBelemTet for tetrahedra" );

  std::map< int, std::vector< std::size_t > > bface;

  // Create unstructured-mesh object to read into
  tk::UnsMesh inmesh;
  // Read in mesh from file
  std::string infile( tk::regression_dir() +
                      "/meshconv/gmsh_output/box_24_ss1.exo" );
  tk::ExodusIIMeshReader er( infile );
  std::map< int, std::vector< std::size_t > > faceid;
  er.readSidesetFaces( bface, faceid );

  auto nbfac = tk::sumvalsize(bface);

  // Test if the number of boundary faces is correct
  ensure_equals( "total number of boundary faces incorrect",
                 nbfac, 24 );

  // Mesh connectivity for simple tetrahedron-only mesh
  std::vector< std::size_t > inpoel { 12, 14,  9, 11,
                                      10, 14, 13, 12,
                                      14, 13, 12,  9,
                                      10, 14, 12, 11,
                                       1, 14,  5, 11,
                                       7,  6, 10, 12,
                                      14,  8,  5, 10,
                                       8,  7, 10, 13,
                                       7, 13,  3, 12,
                                       1,  4, 14,  9,
                                      13,  4,  3,  9,
                                       3,  2, 12,  9,
                                       4,  8, 14, 13,
                                       6,  5, 10, 11,
                                       1,  2,  9, 11,
                                       2,  6, 12, 11,
                                       6, 10, 12, 11,
                                       2, 12,  9, 11,
                                       5, 14, 10, 11,
                                      14,  8, 10, 13,
                                      13,  3, 12,  9,
                                       7, 10, 13, 12,
                                      14,  4, 13,  9,
                                      14,  1,  9, 11 };

  // Boundary-face node connectivity
  std::vector< std::size_t > triinpoel {  2,  9,  3,
                                          1,  9,  2,
                                          3,  9,  4,
                                          1,  4,  9,
                                          7, 10,  6,
                                          6, 10,  5,
                                          8, 10,  7,
                                         10,  8,  5,
                                          6, 11,  2,
                                          2, 11,  1,
                                         11,  6,  5,
                                         11,  5,  1,
                                          3, 12,  2,
                                         12,  6,  2,
                                          7, 12,  3,
                                         12,  7,  6,
                                         13,  7,  3,
                                          4, 13,  3,
                                         13,  8,  7,
                                          8, 13,  4,
                                          5,  8, 14,
                                          1,  5, 14,
                                          4, 14,  8,
                                          1, 14,  4 };

  // Correct boundary elements
  std::vector< std::size_t > correct_belem{ 12,
                                            15,
                                            11,
                                            10,
                                             6,
                                            14,
                                             8,
                                             7,
                                            16,
                                            15,
                                            14,
                                             5,
                                            12,
                                            16,
                                             9,
                                             6,
                                             9,
                                            11,
                                             8,
                                            13,
                                             7,
                                             5,
                                            13,
                                            10 };

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );
  tk::shiftToZero( triinpoel );

  auto esup = tk::genEsup( inpoel, 4 );
  auto esuel = tk::genEsuelTet( inpoel, esup );
  auto nipfac = tk::genNipfac( 4, nbfac, esuel );
  auto inpofa = tk::genInpofaTet( nipfac, nbfac, inpoel, triinpoel, esuel );
  auto belem = tk::genBelemTet( nbfac, inpofa, esup );

  ensure_equals( "total number of entries in belem is incorrect",
                 belem.size(), correct_belem.size() );

  for(std::size_t i=0 ; i<belem.size(); ++i)
  {
    ensure_equals("incorrect entry " + std::to_string(i) + " in belem",
                    belem[i], correct_belem[i]-1);
  }
}

#if defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

} // tut::

#endif  // DOXYGEN_GENERATING_OUTPUT
