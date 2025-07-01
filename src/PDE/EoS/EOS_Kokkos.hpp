#ifndef EOS_KOKKOS_h
#define EOS_KOKKOS_h

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
class EOS_Kokkos {

  private:
    //! Using union style instead of variant
    enum class EOSType = {StiffenedGas
                , JWL
                , SmallShearSolid
                , WilkinsAluminum
                , GodunovRomenski
                , ThermallyPerfectGas};

    EOSType type;
    union EOSUnion {
        StiffiedGas stiffenedGas;
        JWl jwl;
        SmallShearSolid smallShearSolid;
        WilkinsAluminum wilkinsAluminum;
        GodunovRomenski godunovRomenski;
        ThermallyPerfectGas thermallyPerfectGas;

        KOKKOS_FUNCTION
        EOSUnion() {};

        KOKKOS_FUNCTION
        ~EOSUnion() {};

        KOKKOS_FUNCTION
        EOSUnion(const StiffenedGas& s) : stiffenedGas(s) {};

        KOKKOS_FUNCTION
        EOSUnion(const JWL& s) : jwl(s) {};

        KOKKOS_FUNCTION
        EOSUnion(const SmallShearSolid& s) : smallShearSolid(s) {};

        KOKKOS_FUNCTION
        EOSUnion(const WilkinsAluminum& s) : wilkinsAluminum(s) {};

        KOKKOS_FUNCTION
        EOSUnion(const GodunovRomenskis& s) : godunovRomenski(s) {};

        KOKKOS_FUNCTION
        EOSUnion(const ThermallyPerfectGas& s) : thermallyPerfectGas(s) {};

    } m_material;

  public:

    KOKKOS_FUNCTION
    EOS_Kokkos(): type(EOSType::StiffenedGas), m_material(StiffenedGas()) {};

    KOKKOS_FUNCTION
    EOS_Kokkos(const StiffenedGas& s) : type(EOSType::StiffenedGas), m_material(s) {};

    KOKKOS_FUNCTION
    EOS_Kokkos(const JWL& s) : type(EOSType::JWL), m_material(s) {};
    
    KOKKOS_FUNCTION
    EOS_Kokkos(const SmallShearSolid& s) : type(EOSType::SmallShearSolid), m_material(s) {};
    
    KOKKOS_FUNCTION
    EOS_Kokkos(const WilkinsAluminum& s) : type(EOSType::WilkinsAluminum), m_material(s) {};
    
    KOKKOS_FUNCTION
    EOS_Kokkos(const GodunovRomenski& s) : type(EOSType::GodunovRomenski), m_material(s) {};
    
    KOKKOS_FUNCTION
    EOS_Kokkos(const ThermallyPerfectGas& s) : type(EOSType::ThermallyPerfectGas), m_material(s) {};
    

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
    //! Call EOS function. Similar to regular EOS.hpp, but now
    //! utilize union style polymorphism

    template< typename Fn, typename... Args >
    KOKKOS_FUNCTION
    tk::real compute( Args&&... args ) const {
        if (type == EOSType::StiffenedGas) {
          if constexpr( std::is_same_v< Fn, density > )
            return m_material.stiffenedGas.density( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, pressure > )
            return m.pressure( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, pressure_coldcompr > )
            return m.pressure_coldcompr( std::forward< Args >( args )... );

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

          else if constexpr( std::is_same_v< Fn, gas_constant > )
            return m.gas_constant( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, internalenergy > )
            return m.internalenergy( std::forward< Args >( args )... );

          else if constexpr( std::is_same_v< Fn, cv > )
            return m.cv( std::forward< Args >( args )... );
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

};

} // inciter::

#endif // EOS_h
