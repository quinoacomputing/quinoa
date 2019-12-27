// *****************************************************************************
/*!
  \file      src/Base/TaggedTupleiPrint.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Tagged tuple printer
  \details   Tagged tuple printer.
*/
// *****************************************************************************
#ifndef TaggedTuplePrint_h
#define TaggedTuplePrint_h

#include <ostream>

#include <brigand/algorithms/for_each.hpp>
#include <brigand/sequences/has_key.hpp>

#include "NoWarning/set.hpp"
#include "TaggedTuple.hpp"

namespace tk {

//! Function object type to print contents of a TaggedTuple
//! \tparam List brigand::list of types in the tagged tuple
//! \tparam Ignore brigand::list of types to not print
template< class List, class Ignore = brigand::set<> >
struct TuplePrinter {
  std::ostream& os;
  const tk::TaggedTuple< List >& tuple;
  //! Constructor
  TuplePrinter( std::ostream& s, const tk::TaggedTuple< List >& t ) :
    os(s), tuple(t) {}
  //! Function call operator templated on the type being output
  template< typename Key > void operator()( brigand::type_<Key> ) {
    using ignored = brigand::has_key< Ignore, Key >;
    if constexpr( std::is_same_v< ignored, brigand::false_type > ) {
      using Tuple = tk::TaggedTuple< List >;
      const auto& key = Key::name();
      const auto& value = tuple.template get< Key >();
      if constexpr( Tuple::template is_tagged_tuple< Key >::value )
        os << "{ " << key << ": " << value << "} ";
      else
        os << key << ": '" << value << "' ";
    }
  }
};

//! Simple (unformatted) output of a TaggedTuple to an output stream
//! \tparam List brigand::list of types in the tagged tuple
//! \param os Output stream to output to
//! \param t TaggedTuple to print
//! \return Output stream
template< class List >
inline std::ostream&
operator<< ( std::ostream& os, const tk::TaggedTuple< List >& t ) {
  using keys = typename tk::TaggedTuple< List >::Keys;
  os << std::boolalpha;
  brigand::for_each< keys >( TuplePrinter< List >( os, t ) );
  return os;
}

//! Simple (unformatted) output of a TaggedTuple to an output stream
//! \tparam List brigand::list of types in the tagged tuple
//! \tparam Ignore brigand::list of types to not print
//! \param os Output stream to output to
//! \param t TaggedTuple to print
template< class List, class Ignore >
inline void
print( std::ostream& os, const tk::TaggedTuple< List >& t ) {
  using keys = typename tk::TaggedTuple< List >::Keys;
  os << std::boolalpha;
  brigand::for_each< keys >( TuplePrinter< List, Ignore >( os, t ) );
}

} // tk::

#endif // TaggedTuplePrint_h
