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

//! Interpret inciter box ... end blocks using Lua
//! \tparam eq Equation system index
//! \param[in] iclua Lua::sol table handle for IC block
//! \param[in,out] icbox Box configuration tuple in nested vectors to write to
template< class eq >
void box( const sol::table& iclua,
          std::size_t,
          std::vector< ctr::box >& icbox )
{
  const auto& boxlua = iclua[ kw::box::string() ];
  if (boxlua.valid()) {
    std::size_t i = 1;
    do {

      const auto& s = boxlua[i];
      if (not s.valid()) break;
      icbox.emplace_back();
      auto& b = icbox.back();

      b.template get< tag::xmin >() = s[ kw::xmin::string() ].get_or(0.0);
      b.template get< tag::xmax >() = s[ kw::xmax::string() ].get_or(0.0);
      b.template get< tag::ymin >() = s[ kw::ymin::string() ].get_or(0.0);
      b.template get< tag::ymax >() = s[ kw::ymax::string() ].get_or(0.0);
      b.template get< tag::zmin >() = s[ kw::zmin::string() ].get_or(0.0);
      b.template get< tag::zmax >() = s[ kw::zmax::string() ].get_or(0.0);

      b.template get< tag::density >() =
        s[ kw::density::string() ].get_or( 0.0 );
      b.template get< tag::pressure >() =
        s[ kw::pressure::string() ].get_or( 0.0 );
      b.template get< tag::temperature >() =
        s[ kw::temperature::string() ].get_or( 0.0 );
      b.template get< tag::energy >() =
        s[ kw::energy::string() ].get_or( 0.0 );
      b.template get< tag::velocity >() =
        s[ kw::velocity::string() ].get_or< std::vector< tk::real > >( {} );

      auto op = ctr::Initiate();
      auto is = s[ kw::initiate::string() ].
                  get_or( op.name(ctr::InitiateType::IMPULSE) );
      b.template get< tag::initiate >().get< tag::init >() = op.value( is );

      // TODO: The contents of table 'linear' is not yet parsed.

      ++i;

    } while (true);
  }
}

//! Interpret inciter ic ... end block using Lua
//! \tparam eq Equation system index
//! \param[in] lua Lua library handle
//! \param[in] cnt Lua block counter
//! \param[in,out] deck Input deck to write
template< class eq >
void ic( const sol::state& lua, std::size_t cnt, ctr::InputDeck& deck ) {
  // process all 'ic' tables in lua
  const auto& iclua = lua.get_or< sol::table >( kw::ic::string(), {} );

  if (iclua != sol::lua_nil) {      // if table 'ic' exists in lua
    auto& icblock = deck.get< tag::param, eq, tag::ic >();

    // Lambda to read in and store a vector of doubles from Lua
    auto vec = [&]( const std::string& token,
                    std::vector< std::vector< tk::real > >& v )
    {
      auto r = iclua.get_or< std::vector< tk::real > >( token, {} );
      if (!r.empty()) v.emplace_back( std::move(r) );
    };

    vec( kw::density::string(), icblock.template get<tag::density>() );
    vec( kw::velocity::string(), icblock.template get<tag::velocity>() );
    vec( kw::pressure::string(), icblock.template get<tag::pressure>() );
    vec( kw::energy::string(), icblock.template get<tag::energy>() );
    vec( kw::temperature::string(), icblock.template get<tag::temperature>() );

    // interpret 'box' table within 'ic' table
    box< eq >( iclua, cnt, icblock.template get< tag::box >() );
  }
}

} // ::lua
} // ::inciter

#endif // InciterInputDeckLuaGrammar_h
