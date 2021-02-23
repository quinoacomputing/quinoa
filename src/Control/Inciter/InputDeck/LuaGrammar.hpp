// *****************************************************************************
/*!
  \file      src/Control/Inciter/InputDeck/LuaGrammar.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Inciter's input deck grammar definition using Lua
  \details   Inciter's input deck grammar definition using Lua.
*/
// *****************************************************************************
#ifndef InciterInputDeckLuaGrammar_h
#define InciterInputDeckLuaGrammar_h

#include "Inciter/InputDeck/Grammar.hpp"

namespace inciter {
namespace lua {

//! Interpret inciter box ... end block using Lua
template< class eq >
void box( const sol::table& tbl, ctr::box& icbox ) {

  // get access to the 'box' table within the 'tbl' lua table
  const auto& boxlua = tbl.get_or< sol::table >( kw::box::string(), {} );

  if (boxlua != sol::lua_nil) {      // if table 'box' exists in table
    using v = std::vector< tk::real >;

    icbox.template get< tag::xmin >() =
      boxlua.get_or( kw::xmin::string(), 0.0 );
    icbox.template get< tag::xmax >() =
      boxlua.get_or( kw::xmax::string(), 0.0 );
    icbox.template get< tag::ymin >() =
      boxlua.get_or( kw::ymin::string(), 0.0 );
    icbox.template get< tag::ymax >() =
      boxlua.get_or( kw::ymax::string(), 0.0 );
    icbox.template get< tag::zmin >() =
      boxlua.get_or( kw::zmin::string(), 0.0 );
    icbox.template get< tag::zmax >() =
      boxlua.get_or( kw::zmax::string(), 0.0 );

    icbox.template get< tag::density >().back() =
      boxlua.get_or< v >( kw::density::string(), {} );
    icbox.template get< tag::velocity >().back() =
      boxlua.get_or< v >( kw::velocity::string(), {} );
    icbox.template get< tag::pressure >().back() =
      boxlua.get_or< v >( kw::pressure::string(), {} );
    icbox.template get< tag::energy >().back() =
      boxlua.get_or< v >( kw::energy::string(), {} );
    icbox.template get< tag::temperature >().back() =
      boxlua.get_or< v >( kw::temperature::string(), {} );
  }
}

//! Interpret inciter ic ... end block using Lua
template< class eq >
void ic( const sol::state& lua, ctr::InputDeck& deck ) {

  // get access to the 'ic' table in lua
  const auto& iclua = lua.get_or< sol::table >( kw::ic::string(), {} );

  if (iclua != sol::lua_nil) {      // if table 'ic' exists in lua
    using v = std::vector< tk::real >;
    auto& icblock = deck.get< tag::param, eq, tag::ic >();

    icblock.template get< tag::density >().back() =
      iclua.get_or< v >( kw::density::string(), {} );
    icblock.template get< tag::velocity >().back() =
      iclua.get_or< v >( kw::velocity::string(), {} );
    icblock.template get< tag::pressure >().back() =
      iclua.get_or< v >( kw::pressure::string(), {} );
    icblock.template get< tag::energy >().back() =
      iclua.get_or< v >( kw::energy::string(), {} );
    icblock.template get< tag::temperature >().back() =
      iclua.get_or< v >( kw::temperature::string(), {} );

    // interpret 'box' table within 'ic' table
    box< eq >( iclua, icblock.template get< tag::box >() );
  }
}

} // ::lua
} // ::inciter

#endif // InciterInputDeckLuaGrammar_h
