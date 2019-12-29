// *****************************************************************************
/*!
  \file      src/Control/CmdLinePrint.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     CmdLine-specific pretty printer functionality
  \details   CmdLine-specific pretty printer functionality.
*/
// *****************************************************************************
#ifndef CmdLinePrint_h
#define CmdLinePrint_h

#include <ostream>

#include <brigand/algorithms/for_each.hpp>
#include <brigand/sequences/has_key.hpp>

#include "NoWarning/format.hpp"
#include "NoWarning/set.hpp"

#include "TaggedTuple.hpp"

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
        os << " }";
        --depth;
      } else {
        std::string indent( depth * 2, ' ' );
        os << boost::format("\n%s%-10s : %b") % indent % key % value;
      }
    }
  }
};

//! Output command line object (a TaggedTuple) to file
//! \tparam CmdLine Command line object type
//! \param[in,out] os Output stream to print to
//! \param[in] name Name of the command line, e.g., executable name
//! \param[in] c Command line object to output to file
template< class CmdLine >
void print( std::ostream& os,
            const std::string& name,
            const CmdLine& c )
{
  using Keys = typename CmdLine::Keys;
  using Ignore = typename CmdLine::ignore;
  using List = typename CmdLine::PairList;
  os << "# vim: filetype=sh:\n#\n"
        "# Contents of a tagged tuple.\n#\n"
        "# The first string is the name of the tuple followed by its type in\n"
        "# double quotes. A string in single quotes denote the name/tag of a\n"
        "# (nested) tagged tuple. The contents of tuples are enclosed within\n"
        "# braces, indented, and aligned to the same column, compared to the\n"
        "# parent tuple.\n\n";
  os << name << " \"cmdline\" {";
  std::size_t depth = 1;
  brigand::for_each< Keys >( DeepTuplePrinter< List, Ignore >( os, c, depth ) );
  os << " }";
}

} // tk::

#endif // CmdLinePrint_h
