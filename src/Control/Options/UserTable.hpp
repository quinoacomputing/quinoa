// *****************************************************************************
/*!
  \file      src/Control/Options/UserTable.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     User-defined table (discrete y=f(x) function) options
  \details   User-defined table (discrete y=f(x) function) options.
*/
// *****************************************************************************
#ifndef UserTableOptions_h
#define UserTableOptions_h

#include <brigand/sequences/list.hpp>

#include "Toggle.hpp"
#include "Keywords.hpp"
#include "PUPUtil.hpp"

namespace tk {
namespace ctr {

//! Table types types
enum class UserTableType : uint8_t { POSITION=0,
                                     VELOCITY,
                                     ACCELERATION };

//! \brief Pack/Unpack UserTableType: forward overload to generic enum class
//!   packer
inline void operator|( PUP::er& p, UserTableType& e ) { PUP::pup( p, e ); }

//! \brief UserTable options: outsource searches to base templated on enum
//!   type
class UserTable : public tk::Toggle< UserTableType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    using keywords = brigand::list< kw::position
                                  , kw::velocity
                                  , kw::acceleration
                                  >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit UserTable() :
      tk::Toggle< UserTableType >(
        //! Group, i.e., options, name
        "User-defined table",
        //! Enums -> names
        { { UserTableType::POSITION, kw::position::name() },
          { UserTableType::VELOCITY, kw::velocity::name() },
          { UserTableType::ACCELERATION, kw::acceleration::name() } },
        //! keywords -> Enums
        { { kw::position::string(), UserTableType::POSITION },
          { kw::velocity::string(), UserTableType::VELOCITY },
          { kw::acceleration::string(), UserTableType::ACCELERATION } } ) {}
};

} // ctr::
} // tk:::

#endif // UserTableOptions_h
