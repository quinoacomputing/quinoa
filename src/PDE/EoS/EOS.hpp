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
#include "EoS/GodunovRomenski.hpp"
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
   //! Using union style instead of variant
    enum class EOSType {StiffenedGas
                , JWL
                , SmallShearSolid
                , WilkinsAluminum
                , GodunovRomenski
                , ThermallyPerfectGas};

    EOSType type;
    union {
        StiffenedGas stiffenedGas;
        JWL jwl;
        SmallShearSolid smallShearSolid;
        WilkinsAluminum wilkinsAluminum;
        GodunovRomenski godunovRomenski;
        ThermallyPerfectGas thermallyPerfectGas;
    } m_material;


  public:
    //! Empty constructor for Charm++
    explicit EOS();

    //! Constructor
    explicit EOS( ctr::MaterialType mattype, EqType eq, std::size_t k );

    //! Entry method tags for specific EOS classes to use with compute()
    struct density {};
    struct pressure {};
    struct pressure_coldcompr {};
    struct soundspeed {};
    struct shearspeed {};
    struct totalenergy {};
    struct temperature {};
    struct min_eff_pressure {};
    struct refDensity {};
    struct refPressure {};
    struct rho0 {};
    struct gas_constant {};
    struct internalenergy {};
    struct cv {};
    //! Call EOS function
    //! \tparam Fn Function tag identifying the function to call
    //! \tparam Args Types of arguments to pass to function
    //! \param[in] args Arguments to member function to be called
    //! \details This function issues a call to a member function of the
    //!   EOS vector and is thus equivalent to mat_blk[imat].Fn(...).
    template< typename Fn, typename... Args >
    tk::real compute( Args&&... args ) const {
        if (type == EOSType::StiffenedGas) {
          if constexpr( std::is_same_v< Fn, density > )
            return m_material.stiffenedGas.density( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, pressure > )
            return m_material.stiffenedGas.pressure( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, pressure_coldcompr > )
            return m_material.stiffenedGas.pressure_coldcompr( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, soundspeed > )
            return m_material.stiffenedGas.soundspeed( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, shearspeed > )
            return m_material.stiffenedGas.shearspeed( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, totalenergy > )
            return m_material.stiffenedGas.totalenergy( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, temperature > )
            return m_material.stiffenedGas.temperature( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, min_eff_pressure > )
            return m_material.stiffenedGas.min_eff_pressure( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, refDensity > )
            return m_material.stiffenedGas.refDensity( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, refPressure > )
            return m_material.stiffenedGas.refPressure( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, rho0 > )
            return m_material.stiffenedGas.rho0( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, gas_constant > )
            return m_material.stiffenedGas.gas_constant( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, internalenergy > )
            return m_material.stiffenedGas.internalenergy( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, cv > )
            return m_material.stiffenedGas.cv( std::forward< Args >( args )... );
        }
        else if (type == EOSType::JWL) {
          if constexpr( std::is_same_v< Fn, density > )
            return m_material.jwl.density( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, pressure > )
            return m_material.jwl.pressure( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, pressure_coldcompr > )
            return m_material.jwl.pressure_coldcompr( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, soundspeed > )
            return m_material.jwl.soundspeed( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, shearspeed > )
            return m_material.jwl.shearspeed( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, totalenergy > )
            return m_material.jwl.totalenergy( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, temperature > )
            return m_material.jwl.temperature( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, min_eff_pressure > )
            return m_material.jwl.min_eff_pressure( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, refDensity > )
            return m_material.jwl.refDensity( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, refPressure > )
            return m_material.jwl.refPressure( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, rho0 > )
            return m_material.jwl.rho0( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, gas_constant > )
            return m_material.jwl.gas_constant( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, internalenergy > )
            return m_material.jwl.internalenergy( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, cv > )
            return m_material.jwl.cv( std::forward< Args >( args )... );
        }
        else if (type == EOSType::SmallShearSolid) {
          if constexpr( std::is_same_v< Fn, density > )
            return m_material.smallShearSolid.density( std::forward< Args >( args )... );
          else if constexpr( std::is_same_v< Fn, pressure > )
            return m_material.smallShearSolid.pressure( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, pressure_coldcompr > )
            return m_material.smallShearSolid.pressure_coldcompr( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, soundspeed > )
            return m_material.smallShearSolid.soundspeed( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, shearspeed > )
            return m_material.smallShearSolid.shearspeed( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, totalenergy > )
            return m_material.smallShearSolid.totalenergy( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, temperature > )
            return m_material.smallShearSolid.temperature( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, min_eff_pressure > )
            return m_material.smallShearSolid.min_eff_pressure( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, refDensity > )
            return m_material.smallShearSolid.refDensity( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, refPressure > )
            return m_material.smallShearSolid.refPressure( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, rho0 > )
            return m_material.smallShearSolid.rho0( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, gas_constant > )
            return m_material.smallShearSolid.gas_constant( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, internalenergy > )
            return m_material.smallShearSolid.internalenergy( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, cv > )
            return m_material.smallShearSolid.cv( std::forward< Args >( args )... );
        } 
        else if (type == EOSType::GodunovRomenski) {
          if constexpr( std::is_same_v< Fn, density > )
            return m_material.godunovRomenski.density( std::forward< Args >( args )... );
          else if constexpr( std::is_same_v< Fn, pressure > )
            return m_material.godunovRomenski.pressure( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, pressure_coldcompr > )
            return m_material.godunovRomenski.pressure_coldcompr( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, soundspeed > )
            return m_material.godunovRomenski.soundspeed( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, shearspeed > )
            return m_material.godunovRomenski.shearspeed( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, totalenergy > )
            return m_material.godunovRomenski.totalenergy( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, temperature > )
            return m_material.godunovRomenski.temperature( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, min_eff_pressure > )
            return m_material.godunovRomenski.min_eff_pressure( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, refDensity > )
            return m_material.godunovRomenski.refDensity( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, refPressure > )
            return m_material.godunovRomenski.refPressure( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, rho0 > )
            return m_material.godunovRomenski.rho0( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, gas_constant > )
            return m_material.godunovRomenski.gas_constant( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, internalenergy > )
            return m_material.godunovRomenski.internalenergy( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, cv > )
            return m_material.godunovRomenski.cv( std::forward< Args >( args )... );
        }
        else if (type == EOSType::ThermallyPerfectGas) {
          if constexpr( std::is_same_v< Fn, density > )
            return m_material.thermallyPerfectGas.density( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, pressure > )
            return m_material.thermallyPerfectGas.pressure( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, pressure_coldcompr > )
            return m_material.thermallyPerfectGas.pressure_coldcompr( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, soundspeed > )
            return m_material.thermallyPerfectGas.soundspeed( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, shearspeed > )
            return m_material.thermallyPerfectGas.shearspeed( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, totalenergy > )
            return m_material.thermallyPerfectGas.totalenergy( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, temperature > )
            return m_material.thermallyPerfectGas.temperature( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, min_eff_pressure > )
            return m_material.thermallyPerfectGas.min_eff_pressure( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, refDensity > )
            return m_material.thermallyPerfectGas.refDensity( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, refPressure > )
            return m_material.thermallyPerfectGas.refPressure( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, rho0 > )
            return m_material.thermallyPerfectGas.rho0( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, gas_constant > )
            return m_material.thermallyPerfectGas.gas_constant( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, internalenergy > )
            return m_material.thermallyPerfectGas.internalenergy( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, cv > )
            return m_material.thermallyPerfectGas.cv( std::forward< Args >( args )... );
        }
        else if (type == EOSType::WilkinsAluminum) {
          if constexpr( std::is_same_v< Fn, density > )
            return m_material.wilkinsAluminum.density( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, pressure > )
            return m_material.wilkinsAluminum.pressure( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, pressure_coldcompr > )
            return m_material.wilkinsAluminum.pressure_coldcompr( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, soundspeed > )
            return m_material.wilkinsAluminum.soundspeed( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, shearspeed > )
            return m_material.wilkinsAluminum.shearspeed( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, totalenergy > )
            return m_material.wilkinsAluminum.totalenergy( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, temperature > )
            return m_material.wilkinsAluminum.temperature( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, min_eff_pressure > )
            return m_material.wilkinsAluminum.min_eff_pressure( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, refDensity > )
            return m_material.wilkinsAluminum.refDensity( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, refPressure > )
            return m_material.wilkinsAluminum.refPressure( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, rho0 > )
            return m_material.wilkinsAluminum.rho0( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, gas_constant > )
            return m_material.wilkinsAluminum.gas_constant( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, internalenergy > )
            return m_material.wilkinsAluminum.internalenergy( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, cv > )
            return m_material.wilkinsAluminum.cv( std::forward< Args >( args )... );
        };

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
        if (type == EOSType::StiffenedGas) {
          if constexpr(std::is_same_v<Fn, CauchyStress>) 
            return m_material.stiffenedGas.CauchyStress( std::forward< Args >( args )... );
        }
        else if (type == EOSType::GodunovRomenski) {
          if constexpr(std::is_same_v<Fn, CauchyStress>) 
            return m_material.godunovRomenski.CauchyStress( std::forward< Args >( args )... );
        }
        else if (type == EOSType::JWL) {
          if constexpr(std::is_same_v<Fn, CauchyStress>) 
            return m_material.jwl.CauchyStress( std::forward< Args >( args )... );
        }
        else if (type == EOSType::SmallShearSolid) {
          if constexpr(std::is_same_v<Fn, CauchyStress>) 
            return m_material.smallShearSolid.CauchyStress( std::forward< Args >( args )... );
        }
        else if (type == EOSType::ThermallyPerfectGas) {
          if constexpr(std::is_same_v<Fn, CauchyStress>) 
            return m_material.thermallyPerfectGas.CauchyStress( std::forward< Args >( args )... );
        }
        else if (type == EOSType::WilkinsAluminum) {
          if constexpr(std::is_same_v<Fn, CauchyStress>) 
            return m_material.wilkinsAluminum.CauchyStress( std::forward< Args >( args )... );
        }
      
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
       if (type == EOSType::StiffenedGas) {
          if constexpr(std::is_same_v<Fn, CauchyStress>) 
            return m_material.stiffenedGas.setRho0( std::forward< Args >( args )... );
        }
        else if (type == EOSType::GodunovRomenski) {
          if constexpr(std::is_same_v<Fn, CauchyStress>) 
            return m_material.godunovRomenski.setRho0( std::forward< Args >( args )... );
        }
        else if (type == EOSType::JWL) {
          if constexpr(std::is_same_v<Fn, CauchyStress>) 
            return m_material.jwl.setRho0( std::forward< Args >( args )... );
        }
        else if (type == EOSType::SmallShearSolid) {
          if constexpr(std::is_same_v<Fn, CauchyStress>) 
            return m_material.smallShearSolid.setRho0( std::forward< Args >( args )... );
        }
        else if (type == EOSType::ThermallyPerfectGas) {
          if constexpr(std::is_same_v<Fn, CauchyStress>) 
            return m_material.thermallyPerfectGas.setRho0( std::forward< Args >( args )... );
        }
        else if (type == EOSType::WilkinsAluminum) {
          if constexpr(std::is_same_v<Fn, CauchyStress>) 
            return m_material.wilkinsAluminum.setRho0( std::forward< Args >( args )... );
        }
    }

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup(PUP::er &p) {
    p | type;  // Serialize the tag first

    switch(type) {
        case EOSType::StiffenedGas:
            if (p.isUnpacking()) {
                // Construct in-place since union member needs explicit construction
                new (&m_material.stiffenedGas) StiffenedGas();
            }
            p | m_material.stiffenedGas;
            break;

        case EOSType::JWL:
            if (p.isUnpacking()) {
                new (&m_material.jwl) JWL();
            }
            p | m_material.jwl;
            break;

        case EOSType::SmallShearSolid:
            if (p.isUnpacking()) {
                new (&m_material.smallShearSolid) SmallShearSolid();
            }
            p | m_material.smallShearSolid;
            break;

        case EOSType::WilkinsAluminum:
            if (p.isUnpacking()) {
                new (&m_material.wilkinsAluminum) WilkinsAluminum();
            }
            p | m_material.wilkinsAluminum;
            break;

        case EOSType::GodunovRomenski:
            if (p.isUnpacking()) {
                new (&m_material.godunovRomenski) GodunovRomenski();
            }
            p | m_material.godunovRomenski;
            break;

        case EOSType::ThermallyPerfectGas:
            if (p.isUnpacking()) {
                new (&m_material.thermallyPerfectGas) ThermallyPerfectGas();
            }
            p | m_material.thermallyPerfectGas;
            break;
    }
}

    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] s EOS object reference
    friend void operator|( PUP::er& p, EOS& s ) { s.pup(p); }
    //@}
};

} // inciter::

#endif // EOS_h
