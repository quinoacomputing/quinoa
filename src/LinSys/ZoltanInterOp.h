//******************************************************************************
/*!
  \file      src/LinSys/ZoltanInterOp.h
  \author    J. Bakosi
  \date      Wed 15 Apr 2015 09:22:35 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Interoperation with the Zoltan library
  \details   Interoperation with the Zoltan library, used for static mesh
    partitioning.
*/
//******************************************************************************
#ifndef ZoltanInterOp_h
#define ZoltanInterOp_h

#include <UnsMesh.h>
#include <Tags.h>
#include <TaggedTuple.h>

namespace tk {
//! Interoperation with the Zoltan library, used for static mesh partitioning
namespace zoltan {

//! Partition mesh using Zoltan's hypergraph algorithm in serial
std::vector< std::size_t >
partitionMesh( const tk::UnsMesh& mesh,
               uint64_t npart,
               const tk::Print& print );

} // zoltan::
} // tk::

#endif // ZoltanInterOp_h
