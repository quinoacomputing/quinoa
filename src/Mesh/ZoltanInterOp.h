//******************************************************************************
/*!
  \file      src/Mesh/ZoltanInterOp.h
  \author    J. Bakosi
  \date      Thu 02 Apr 2015 02:34:31 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Interoperation with the Zoltan library
  \details   Interoperation with the Zoltan library, used for static mesh
    partitioning.
*/
//******************************************************************************
#ifndef ZoltanInterOp_h
#define ZoltanInterOp_h

#include <UnsMesh.h>

namespace tk {
//! Interoperation with the Zoltan library, used for static mesh partitioning
namespace zoltan {

//! Partition mesh using Zoltan
void partitionMesh( const tk::UnsMesh& mesh );

} // zoltan::
} // tk::

#endif // ZoltanInterOp_h
