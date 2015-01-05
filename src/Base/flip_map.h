//******************************************************************************
/*!
  \file      src/Base/flip_map.h
  \author    J. Bakosi
  \date      Sun 21 Dec 2014 10:04:01 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Flip a std::map yielding a multimap sorted by std::map::value_type
  \details   Flip a std::map yielding a multimap sorted by std::map::value_type.
    Credit goes to Oli Charlesworth: http://stackoverflow.com/a/5056797
*/
//******************************************************************************
#ifndef flip_map_h
#define flip_map_h

namespace tk {

//! Flip a std::pair of arbitrary types
//! \param[in] p std::pair of arbitrary types, A and B
//! \return std::pair of arbitrary types, B and A
//! \author J. Bakosi
template< typename A, typename B >
std::pair< B, A > flip_pair( const std::pair< A ,B >& p )
{ return std::pair< B, A >( p.second, p.first ); }

//! Flip a std::map of arbitrary types, yielding a std::multimap sorted by
//! std::map::value_type.
//! \param[in] src std::map of arbitrary key and value pairs of types A and B
//! \return std::multimap of arbitrary key and value pairs of types B and A
//! \author J. Bakosi
template< typename A, typename B >
std::multimap< B, A > flip_map( const std::map< A, B >& src ) {
  std::multimap< B, A > dst;
  std::transform( src.begin(), src.end(), std::inserter(dst, dst.begin()), 
                  flip_pair< A ,B > );
  return dst;
}

} // tk::

#endif // flip_map_h
