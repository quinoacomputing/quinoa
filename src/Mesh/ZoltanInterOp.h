//******************************************************************************
/*!
  \file      src/Mesh/ZoltanInterOp.h
  \author    J. Bakosi
  \date      Mon 06 Apr 2015 01:48:49 PM MDT
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
tk::tuple::tagged_tuple< tag::psup,  std::pair< std::vector< std::size_t >,
                                                std::vector< std::size_t > >,
                         tag::owner, std::vector< int > >
partitionMesh( const tk::UnsMesh& mesh, uint64_t npart );

} // zoltan::
} // tk::

#endif // ZoltanInterOp_h
