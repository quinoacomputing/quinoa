//******************************************************************************
/*!
  \file      src/LinSys/ZoltanInterOp.h
  \author    J. Bakosi
  \date      Fri 20 Nov 2015 10:57:58 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Interoperation with the Zoltan library
  \details   Interoperation with the Zoltan library, used for static mesh
    partitioning.
*/
//******************************************************************************
#ifndef ZoltanInterOp_h
#define ZoltanInterOp_h

#include <vector>
#include <utility>
#include <cstddef>
#include <cstdint>

#include "Options/PartitioningAlgorithm.h"

namespace tk {

class Print;
class UnsMesh;

//! Interoperation with the Zoltan library, used for static mesh partitioning
namespace zoltan {

//! Partition mesh using Zoltan2 with a geometric partitioner, such as RCB, RIB
std::vector< std::size_t >
geomPartMesh( tk::ctr::PartitioningAlgorithmType algorithm,
              const std::array< std::vector< tk::real >, 3 >& elemcoord,
              const std::vector< std::size_t >& elemid,
              std::size_t nelem,
              uint64_t npart );

} // zoltan::
} // tk::

#endif // ZoltanInterOp_h
