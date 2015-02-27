//******************************************************************************
/*!
  \file      src/Mesh/ZoltanInterOp.h
  \author    J. Bakosi
  \date      Mon 23 Feb 2015 03:31:01 PM MST
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

void partitionMesh( const tk::UnsMesh& mesh );

} // zoltan::
} // tk::

#endif // ZoltanInterOp_h
