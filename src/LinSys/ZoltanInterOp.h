//******************************************************************************
/*!
  \file      src/LinSys/ZoltanInterOp.h
  \author    J. Bakosi
  \date      Mon 16 Nov 2015 02:07:59 PM MST
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

//! Partition mesh using Zoltan's hypergraph algorithm in parallel
std::vector< std::size_t >
partitionMesh( const tk::UnsMesh& graph,
               const std::array< std::vector< tk::real >, 3 >& elemcoord,
               const std::vector< std::size_t > elemid,
               std::size_t nelem,
               uint64_t npart,
               const tk::Print& print );

} // zoltan::
} // tk::

#endif // ZoltanInterOp_h
