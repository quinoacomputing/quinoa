//******************************************************************************
/*!
  \file      src/Mesh/HypreInterOp.h
  \author    J. Bakosi
  \date      Thu 16 Apr 2015 05:36:06 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Interoperation with the Hypre library
  \details   Interoperation with the Hypre library, used for linear solvers.
*/
//******************************************************************************
#ifndef HypreInterOp_h
#define HypreInterOp_h

#include <UnsMesh.h>
#include <Tags.h>
#include <TaggedTuple.h>

namespace tk {
//! Interoperation with the Hypre library, used for linear solvers
namespace hypre {

void test( /*const std::map< int, std::vector< std::size_t > >& contribute*/ );

} // hypre::
} // tk::

#endif // HypreInterOp_h
