// *****************************************************************************
/*!
  \file      src/Base/TaggedTupleDeepPrint.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Structured TaggedTuple printer with depth/indentation
  \details   Structured TaggedTuple printer with depth/indentation.
*/
// *****************************************************************************
#ifndef TaggedTupleDeepPrint_h
#define TaggedTupleDeepPrint_h

#include <ostream>

#include <brigand/algorithms/for_each.hpp>
#include <brigand/sequences/has_key.hpp>

#include "NoWarning/format.hpp"
#include "NoWarning/set.hpp"

#include "Has.hpp"

// The include order here is important: it populates the overloads of
// operator<< for various types, followed by TaggedTuple, the (simple)
// TaggedTuplePrint (which will be accessible by the upstream, simpler
// operator<< for vector, map, etc.) and finally, the most complex TaggedTuple
// printer with depth, defined below.
#include "PrintUtil.hpp"
#include "TaggedTuple.hpp"
#include "TaggedTuplePrint.hpp"

namespace tk {

//! Function object type to print contents of a TaggedTuple at depth
//! \details Compared to tk::TuplePrinter this prints every key and value
//!   in a new line and nested tagged tuples starts at increasing depths
//!   (indents).
//! \tparam List brigand::list of types in the tagged tuple
//! \tparam Ignore brigand::list of types to not print
template< class List, class Ignore = brigand::set<> >
struct DeepTuplePrinter {
  std::ostream& os;
  const tk::TaggedTuple< List >& tuple;
  std::size_t& depth;
  //! Constructor
  DeepTuplePrinter( std::ostream& s, const tk::TaggedTuple< List >& t,
                    std::size_t& d ) : os(s), tuple(t), depth(d) {}
  //! Function call operator templated on the type being output
  template< typename Key > void operator()( brigand::type_<Key> ) {
    using ignored = brigand::has_key< Ignore, Key >;
    if constexpr( std::is_same_v< ignored, brigand::false_type > ) {
      using Tuple = tk::TaggedTuple< List >;
      const auto& key = Key::name();
      const auto& value = tuple.template get< Key >();
      if constexpr( Tuple::template is_tagged_tuple< Key >::value ) {
        std::string indent( depth * 2, ' ' );
        os << '\n' << indent << '\'' << key << "' {";
        using ituple = typename Tuple::template TupleElement< Key >;
        using ikeys = typename ituple::Keys;
        using ilist = typename ituple::PairList;
        brigand::for_each< ikeys >(
          DeepTuplePrinter< ilist, Ignore >( os, value, ++depth ) );
        os << '\n' << indent << '}';
        --depth;
      } else {
        std::string indent( depth * 2, ' ' );
        os << boost::format("\n%s%-15s : ") % indent % key;
        os << std::boolalpha << value;
      }
    }
  }
};

//! Output command line object (a TaggedTuple) to file
//! \tparam Tuple Tuple object type
//! \param[in,out] os Output stream to print to
//! \param[in] name Name of (root) Tuple
//! \param[in] c Command line object to output to file
template< class Tuple >
void print( std::ostream& os, const std::string& name, const Tuple& c ) {
  static_assert( tk::HasTypedef_i_am_tagged_tuple_v< Tuple > );
  using Keys = typename Tuple::Keys;
  using Ignore = typename Tuple::ignore;
  using List = typename Tuple::PairList;
  os << "# vim: filetype=sh:\n#\n"
        "# Contents of a tagged tuple.\n#\n"
        "# A string in single quotes denotes the name/tag of a (nested)\n"
        "# tagged tuple. The contents of tuples are enclosed within braces.\n"
        "# Vectors are enclosed within square brackets. Keys of associative\n"
        "# containers are in paretheses.\n\n";
  os << '\'' << name << "' {";
  std::size_t depth = 1;
  brigand::for_each< Keys >( DeepTuplePrinter< List, Ignore >( os, c, depth ) );
  os << "\n}";
}

} // tk::

#endif // TaggedTupleDeepPrint_h
