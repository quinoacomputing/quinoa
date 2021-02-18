// *****************************************************************************
/*!
  \file      tests/unit/Base/MigratedTypes.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Type definitions for unit tests of Charm++ migrated data in Base
  \details   Type definitions for unit tests of Charm++ migrated data in Base.
*/
// *****************************************************************************
#ifndef test_MigratedTypes_h
#define test_MigratedTypes_h

#include <variant>
#include <optional>

#include "PUPUtil.hpp"
#include "TaggedTuple.hpp"

//! Unit test declarations and definitions
namespace tut {
//! Unit test declarations and definitions for testing Charm++ chares
namespace charm {

//! Tags for migrated tagged tuple test(s)
namespace tag {
  struct name {};
  struct age {};
  struct email {};
} // tag::

// Tested types migrated
enum class Enum_default { F1=32, F2 };
enum class Enum_uint8_t : uint8_t { F1=14, F2 };
enum Enum_cstyle { F1=54, F2 };
using Pair = std::pair< int, double >;
using Vector = std::vector< std::string >;
using Tuple = std::tuple< int, double, std::vector< std::string >,
                          Enum_default, std::map< Enum_uint8_t, std::string > >;
using Array = std::array< int, 2 >;
using UnorderedMap = std::unordered_map< int, std::string >;
using UnorderedSet = std::unordered_set< int >;
using OptionalStr = std::optional< std::string >;
using OptionalInt = std::optional< int >;
using TaggedTuple = tk::TaggedTuple< brigand::list< tag::name,  std::string,
                                                    tag::age,   int,
                                                    tag::email, std::string > >;
using Variant = std::variant< int, double >;

//! Pack/Unpack: delegate to PUP::
inline void operator|( PUP::er& p, Enum_default& e ) { PUP::pup( p, e ); }
inline void operator|( PUP::er& p, Enum_uint8_t& e ) { PUP::pup( p, e ); }
inline void operator|( PUP::er& p, Enum_cstyle& e ) { PUP::pup( p, e ); }
inline void operator|( PUP::er& p, Tuple& t ) { PUP::pup( p, t ); }
inline void operator|( PUP::er& p, Array& a ) { PUP::pup( p, a ); }
inline void operator|( PUP::er& p, UnorderedMap& m ) { PUP::pup( p, m ); }
inline void operator|( PUP::er& p, UnorderedSet& s ) { PUP::pup( p, s ); }
inline void operator|( PUP::er& p, OptionalStr& o ) { PUP::pup( p, o ); }
inline void operator|( PUP::er& p, OptionalInt& o ) { PUP::pup( p, o ); }

} // charm::
} // tut::

#endif // test_MigratedTypes_h
