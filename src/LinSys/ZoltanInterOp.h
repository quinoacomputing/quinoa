//******************************************************************************
/*!
  \file      src/LinSys/ZoltanInterOp.h
  \author    J. Bakosi
  \date      Sat 30 May 2015 11:47:56 AM MDT
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

namespace tk {

class Print;
class UnsMesh;

//! Interoperation with the Zoltan library, used for static mesh partitioning
namespace zoltan {

//! Partition mesh using Zoltan's hypergraph algorithm in serial
std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
partitionMesh( tk::UnsMesh& graph,
               uint64_t npart,
               const tk::Print& print );

} // zoltan::
} // tk::

#endif // ZoltanInterOp_h
