// *****************************************************************************
/*!
  \file      src/LoadBalance/ZoltanInterOp.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Interoperation with the Zoltan library
  \details   Interoperation with the Zoltan library, used for static mesh
    partitioning.
*/
// *****************************************************************************
#ifndef ZoltanInterOp_h
#define ZoltanInterOp_h

#include <vector>
#include <utility>
#include <cstddef>
#include <cstdint>

#include "Types.hpp"
#include "Options/PartitioningAlgorithm.hpp"

namespace tk {

class Print;
class UnsMesh;

//! Interoperation with the Zoltan library, used for static mesh partitioning
namespace zoltan {

//! Partition mesh using Zoltan2 with a geometric partitioner, such as RCB, RIB
std::vector< std::size_t >
geomPartMesh( tk::ctr::PartitioningAlgorithmType algorithm,
              const std::array< std::vector< tk::real >, 3 >& elemcoord,
              const std::vector< long long >& elemid,
              int npart );

} // zoltan::
} // tk::

#endif // ZoltanInterOp_h
