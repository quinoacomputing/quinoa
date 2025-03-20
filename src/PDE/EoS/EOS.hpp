// *****************************************************************************
/*!
  \file      src/PDE/EoS/EOS.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Polymorphic variant-style implementation for equations of state,
    where children implement specific EOS functions.
*/
// *****************************************************************************
#ifndef EOS_h
#define EOS_h

#include <variant>

#include "PUPUtil.hpp"
#include "Inciter/Options/Material.hpp"
#include "EoS/StiffenedGas.hpp"
#include "EoS/JWL.hpp"
#include "EoS/SmallShearSolid.hpp"
#include "EoS/WilkinsAluminum.hpp"
#include "EoS/ThermallyPerfectGas.hpp"

namespace inciter {

//! Equation types
enum class EqType : uint8_t { compflow
                            , multimat
                            , multispecies
                            };

//! Base class for generic forwarding interface to eos types
class EOS {

  private:
    //! Variant type listing all eos types modeling the same concept
    std::variant< StiffenedGas
                , JWL
                , SmallShearSolid
                , WilkinsAluminum
                , ThermallyPerfectGas
                > m_material;

  public:
    //! Empty constructor for Charm++
    explicit EOS() {}

    //! Constructor
    explicit EOS( ctr::MaterialType mattype, EqType eq, std::size_t k );

    //! Entry method tags for specific EOS classes to use with compute()
    struct density {};
    struct pressure {};
    struct soundspeed {};
    struct shearspeed {};
    struct totalenergy {};
    struct temperature {};
    struct min_eff_pressure {};
    struct refDensity {};
    struct refPressure {};
    struct rho0 {};
    //! Call EOS function
    //! \tparam Fn Function tag identifying the function to call
    //! \tparam Args Types of arguments to pass to function
    //! \param[in] args Arguments to member function to be called
    //! \details This function issues a call to a member function of the
    //!   EOS vector and is thus equivalent to mat_blk[imat].Fn(...).
    template< typename Fn, typename... Args >
    tk::real compute( Args&&... args ) const {
      return std::visit( [&]( const auto& m )-> tk::real {
          if constexpr( std::is_same_v< Fn, density > )
            return m.density( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, pressure > )
            return m.pressure( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, soundspeed > )
            return m.soundspeed( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, shearspeed > )
            return m.shearspeed( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, totalenergy > )
            return m.totalenergy( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, temperature > )
            return m.temperature( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, min_eff_pressure > )
            return m.min_eff_pressure( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, refDensity > )
            return m.refDensity( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, refPressure > )
            return m.refPressure( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, rho0 > )
            return m.rho0( std::forward< Args >( args )... );
        }, m_material );
    }

    //! Entry method tags for specific EOS classes to use with computeTensor()
    struct CauchyStress {};
    //! Call EOS function returning a tensor
    //! \tparam Fn Function tag identifying the function to call
    //! \tparam Args Types of arguments to pass to function
    //! \param[in] args Arguments to member function to be called
    //! \details This function issues a call to a member function of the
    //!   EOS vector and is thus equivalent to mat_blk[imat].Fn(...).
    template< typename Fn, typename... Args >
    std::array< std::array< tk::real, 3 >, 3 > computeTensor( Args&&... args )
    const {
      return std::visit( [&]( const auto& m )->
        std::array< std::array< tk::real, 3 >, 3 > {
          if constexpr( std::is_same_v< Fn, CauchyStress > )
            return m.CauchyStress( std::forward< Args >( args )... );

        }, m_material );
    }

    //! Entry method tags for specific EOS classes to use with set()
    struct setRho0 {};
    //! Call EOS function
    //! \tparam Fn Function tag identifying the function to call
    //! \tparam Args Types of arguments to pass to function
    //! \param[in] args Arguments to member function to be called
    //! \details This function issues a call to a member function of the
    //!   EOS vector and is thus equivalent to mat_blk[imat].Fn(...).
    template< typename Fn, typename... Args >
    void set( Args&&... args ) {
      std::visit( [&]( auto& m )-> void {
          if constexpr( std::is_same_v< Fn, setRho0 > )
            m.setRho0( std::forward< Args >( args )... );
        }, m_material );
    }

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) {
      p | m_material;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] s EOS object reference
    friend void operator|( PUP::er& p, EOS& s ) { s.pup(p); }
    //@}
};

} // inciter::

#endif // EOS_h
