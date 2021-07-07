// *****************************************************************************
/*!
  \file      tests/unit/Inciter/ALE/TestVorticity.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Unit tests for ALE and computing vorticity
  \details   Unit tests for ALE and computing vorticity.
*/
// *****************************************************************************

#include <limits>

#include "NoWarning/tut.hpp"

#include "TUTConfig.hpp"
#include "Types.hpp"
#include "Fields.hpp"
#include "DerivedData.hpp"
#include "Reorder.hpp"
#include "Vector.hpp"
#include "ALE/MeshMotion.hpp"
#include "GmshMeshReader.hpp"

#ifndef DOXYGEN_GENERATING_OUTPUT

namespace tut {

//! All tests in group inherited from this base
struct ALE_common {
  // floating point precision tolerance
  const tk::real pr = std::numeric_limits< tk::real >::epsilon();

  // mesh node coordinates
  std::array< std::vector< tk::real >, 3 > coord {{
    //{{ 0, 1, 1, 0, 0, 1, 1, 0, 0.5, 0.5, 0.5, 1,   0.5, 0 }},
    //{{ 0, 0, 1, 1, 0, 0, 1, 1, 0.5, 0.5, 0,   0.5, 1,   0.5 }},
    //{{ 0, 0, 0, 0, 1, 1, 1, 1, 0,   1,   0.5, 0.5, 0.5, 0.5 }} }};
    {{ 0.5, 0,   0.5, 1   }},
    {{ 0.5, 0.5, 1,   0.5 }},
    {{ 1,   0.5, 0.5, 0.5 }} }};

  // mesh connectivity for simple tetrahedron-only mesh
  std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
                                      //12, 14,  9, 11,
                                      //10, 14, 13, 12,
                                      //14, 13, 12,  9,
                                      //10, 14, 12, 11,
                                      //1,  14,  5, 11,
                                      //7,   6, 10, 12,
                                      //14,  8,  5, 10,
                                      //8,   7, 10, 13,
                                      //7,  13,  3, 12,
                                      //1,   4, 14,  9,
                                      //13,  4,  3,  9,
                                      //3,   2, 12,  9,
                                      //4,   8, 14, 13,
                                      //6,   5, 10, 11,
                                      //1,   2,  9, 11,
                                      //2,   6, 12, 11,
                                      //6,  10, 12, 11,
                                      //2,  12,  9, 11,
                                      //5,  14, 10, 11,
                                      //14,  8, 10, 13,
                                      //13,  3, 12,  9,
                                      //7,  10, 13, 12,
                                      //14,  4, 13,  9,
                                      //14,  1,  9, 11 };
};

//! Test group shortcuts
using ALE_group = test_group< ALE_common, MAX_TESTS_IN_GROUP >;
using ALE_object = ALE_group::object;

//! Define test group
static ALE_group ALE( "Inciter/ALE" );

//! Test definitions for group

} // tut::

#endif  // DOXYGEN_GENERATING_OUTPUT
