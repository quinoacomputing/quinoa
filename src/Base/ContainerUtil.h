//******************************************************************************
/*!
  \file      src/Base/ContainerUtil.h
  \author    J. Bakosi
  \date      Fri 24 Apr 2015 06:05:25 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Container utilities
  \details   Various STL container utilities.
*/
//******************************************************************************
#ifndef ContainerUtil_h
#define ContainerUtil_h

#include <algorithm>

#include <Exception.h>

namespace tk {

//! \brief Make elements of container unique
//! \param[inout] c Container
//! \author  J. Bakosi
template< class Container >
void unique( Container& c ) {
  std::sort( begin(c), end(c) );
  auto it = std::unique( begin(c), end(c) );
  auto d = std::distance( begin(c), it );
  Assert( d >= 0, "Distance must be non-negative in tk::unique()" );
  c.resize( static_cast< std::size_t >( d ) );
}

} // tk::

#endif // ContainerUtil_h
