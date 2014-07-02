//******************************************************************************
/*!
  \file      src/Base/flip_map.h
  \author    J. Bakosi
  \date      Wed 02 Jul 2014 08:19:18 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Flip a std::map yielding a multimap sorted by std::map::value_type
  \details   Flip a std::map yielding a multimap sorted by std::map::value_type
*/
//******************************************************************************
#ifndef flip_map_h
#define flip_map_h

namespace tk {

// Credit goes to Oli Charlesworth: http://stackoverflow.com/a/5056797

template< typename A, typename B >
std::pair< B, A > flip_pair( const std::pair< A ,B >& p )
{ return std::pair< B, A >( p.second, p.first ); }

template< typename A, typename B >
std::multimap< B, A > flip_map( const std::map< A, B >& src ) {
  std::multimap< B, A > dst;
  std::transform( src.begin(), src.end(), std::inserter(dst, dst.begin()), 
                  flip_pair< A ,B > );
  return dst;
}

} // tk::

#endif // flip_map_h
