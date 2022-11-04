// *****************************************************************************
/*!
  \file      src/PDE/EoS/EosVariant.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Polymorphic variant-style implementation for equations of state,
    where children implement specific EOS functions.
*/
// *****************************************************************************
#ifndef EosVariant_h
#define EosVariant_h

#include <variant>

#include "PUPUtil.hpp"
#include "Inciter/Options/Material.hpp"
#include "EoS/SGclass.hpp"
#include "EoS/JWLclass.hpp"

namespace inciter {

//! Base class for generic forwarding interface to eos types
class EOS {

  private:
    //! Variant type listing all eos types modeling the same concept
    std::variant< SGclass
                , JWLclass
                > material;

  public:
    //! Empty constructor for Charm++
    explicit EOS() {}

    //! Constructor
    explicit EOS( ctr::MaterialType mattype,
      std::size_t eqtype,
      std::size_t system,
      std::size_t k );

    //! Entry method tags for specific EOS classes to use with eosCall()
    struct eos_density {};
    struct eos_pressure {};
    struct eos_soundspeed {};
    struct eos_totalenergy {};
    struct eos_temperature {};
    struct min_eff_pressure {};
    //! Call EOS function
    //! \tparam Fn Function tag identifying the function to call
    //! \tparam Args Types of arguments to pass to function
    //! \param[in] args Arguments to member function to be called
    //! \details This function issues a call to a member function of the
    //!   EOS vector and is thus equivalent to mat_blk[imat].Fn(...).
    template< typename Fn, typename... Args >
    tk::real eosCall( Args&&... args ) const {
      return std::visit( [&]( const auto& m )-> tk::real {
          if constexpr( std::is_same_v< Fn, eos_density > )
            return m.eos_density( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, eos_pressure > )
            return m.eos_pressure( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, eos_soundspeed > )
            return m.eos_soundspeed( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, eos_totalenergy > )
            return m.eos_totalenergy( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, eos_temperature > )
            return m.eos_temperature( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, min_eff_pressure > )
            return m.min_eff_pressure( std::forward< Args >( args )... );
        }, material );
    }

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) {
      p | material;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] s EOS object reference
    friend void operator|( PUP::er& p, EOS& s ) { s.pup(p); }
    //@}
};

} // inciter::

#endif // EosVariant_h
