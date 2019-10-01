// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/Problem.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem options for inciter
  \details   Problem options for inciter
*/
// *****************************************************************************
#ifndef ProblemOptions_h
#define ProblemOptions_h

#include <brigand/sequences/list.hpp>
#include <brigand/algorithms/for_each.hpp>

#include "Toggle.hpp"
#include "Keywords.hpp"
#include "PUPUtil.hpp"

namespace inciter {
namespace ctr {

//! Problem types
enum class ProblemType : uint8_t { USER_DEFINED,
                                   SHEAR_DIFF,
                                   VORTICAL_FLOW,
                                   NL_ENERGY_GROWTH,
                                   RAYLEIGH_TAYLOR,
                                   TAYLOR_GREEN,
                                   SLOT_CYL,
                                   GAUSS_HUMP,
                                   CYL_ADVECT,
                                   SOD_SHOCKTUBE,
                                   ROTATED_SOD_SHOCKTUBE,
                                   SEDOV_BLASTWAVE,
                                   INTERFACE_ADVECTION,
                                   GAUSS_HUMP_COMPFLOW,
                                   WATERAIR_SHOCKTUBE };

//! Pack/Unpack ProblemType: forward overload to generic enum class packer
inline void operator|( PUP::er& p, ProblemType& e ) { PUP::pup( p, e ); }

//! \brief Problem options: outsource to base templated on enum type
class Problem : public tk::Toggle< ProblemType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    using keywords = brigand::list< kw::user_defined
                                  , kw::shear_diff
                                  , kw::vortical_flow
                                  , kw::nl_energy_growth
                                  , kw::rayleigh_taylor
                                  , kw::taylor_green
                                  , kw::slot_cyl
                                  , kw::gauss_hump
                                  , kw::cyl_advect
                                  , kw::sod_shocktube
                                  , kw::rotated_sod_shocktube
                                  , kw::sedov_blastwave
                                  , kw::interface_advection
                                  , kw::gauss_hump_compflow
                                  , kw::waterair_shocktube
                                  >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit Problem() :
      tk::Toggle< ProblemType >(
        //! Group, i.e., options, name
        kw::problem::name(),
        //! Enums -> names
        { { ProblemType::USER_DEFINED, kw::user_defined::name() },
          { ProblemType::SHEAR_DIFF, kw::shear_diff::name() },
          { ProblemType::VORTICAL_FLOW, kw::vortical_flow::name() },
          { ProblemType::NL_ENERGY_GROWTH, kw::nl_energy_growth::name() },
          { ProblemType::RAYLEIGH_TAYLOR, kw::rayleigh_taylor::name() },
          { ProblemType::TAYLOR_GREEN, kw::taylor_green::name() },
          { ProblemType::SLOT_CYL, kw::slot_cyl::name() },
          { ProblemType::GAUSS_HUMP, kw::gauss_hump::name() },
          { ProblemType::CYL_ADVECT, kw::cyl_advect::name() },
          { ProblemType::SOD_SHOCKTUBE, kw::sod_shocktube::name() },
          { ProblemType::ROTATED_SOD_SHOCKTUBE,
            kw::rotated_sod_shocktube::name() },
          { ProblemType::SEDOV_BLASTWAVE, kw::sedov_blastwave::name() },
          { ProblemType::INTERFACE_ADVECTION,
            kw::interface_advection::name() },
          { ProblemType::GAUSS_HUMP_COMPFLOW,
            kw::gauss_hump_compflow::name() },
          { ProblemType::WATERAIR_SHOCKTUBE, kw::waterair_shocktube::name() } },
        //! keywords -> Enums
        { { kw::user_defined::string(), ProblemType::USER_DEFINED },
          { kw::shear_diff::string(), ProblemType::SHEAR_DIFF },
          { kw::vortical_flow::string(), ProblemType::VORTICAL_FLOW },
          { kw::nl_energy_growth::string(), ProblemType::NL_ENERGY_GROWTH },
          { kw::rayleigh_taylor::string(), ProblemType::RAYLEIGH_TAYLOR },
          { kw::taylor_green::string(), ProblemType::TAYLOR_GREEN },
          { kw::slot_cyl::string(), ProblemType::SLOT_CYL },
          { kw::gauss_hump::string(), ProblemType::GAUSS_HUMP },
          { kw::cyl_advect::string(), ProblemType::CYL_ADVECT },
          { kw::sod_shocktube::string(), ProblemType::SOD_SHOCKTUBE },
          { kw::rotated_sod_shocktube::string(),
            ProblemType::ROTATED_SOD_SHOCKTUBE },
          { kw::sod_shocktube::string(), ProblemType::SOD_SHOCKTUBE },
          { kw::sedov_blastwave::string(), ProblemType::SEDOV_BLASTWAVE },
          { kw::interface_advection::string(),
            ProblemType::INTERFACE_ADVECTION },
          { kw::gauss_hump_compflow::string(),
            ProblemType::GAUSS_HUMP_COMPFLOW },
          { kw::waterair_shocktube::string(),
            ProblemType::WATERAIR_SHOCKTUBE } } )
    {
       brigand::for_each< keywords >( assertPolicyCodes() );
    }

    //! \brief Return policy code based on Enum
    //! \param[in] p Enum value of the problem option requested
    //! \return Policy code of the option
    const std::string& code( ProblemType p ) const {
      using tk::operator<<;
      auto it = policy.find( p );
      Assert( it != end(policy),
              std::string("Cannot find policy code for problem \"") << p <<
                "\"" );
      return it->second;
    }

  private:
    //! Function object for ensuring the existence of policy codes
    struct assertPolicyCodes {
      //! \brief Function call operator templated on the type to assert the
      //!   existence of a policy code
      template< typename U > void operator()( brigand::type_<U> ) {
        static_assert( tk::HasTypedef_code_v< typename U::info >,
                       "Policy code undefined for keyword" );
      }
    };

    //! Enums -> policy code
    std::map< ProblemType, std::string > policy {
        { ProblemType::USER_DEFINED, *kw::user_defined::code() }
      , { ProblemType::SHEAR_DIFF, *kw::shear_diff::code() }
      , { ProblemType::VORTICAL_FLOW, *kw::vortical_flow::code() }
      , { ProblemType::NL_ENERGY_GROWTH, *kw::nl_energy_growth::code() }
      , { ProblemType::RAYLEIGH_TAYLOR, *kw::rayleigh_taylor::code() }      
      , { ProblemType::TAYLOR_GREEN, *kw::taylor_green::code() }      
      , { ProblemType::SLOT_CYL, *kw::slot_cyl::code() }
      , { ProblemType::GAUSS_HUMP, *kw::gauss_hump::code() }
      , { ProblemType::CYL_ADVECT, *kw::cyl_advect::code() }
      , { ProblemType::SOD_SHOCKTUBE, *kw::sod_shocktube::code() }
      , { ProblemType::ROTATED_SOD_SHOCKTUBE,
          *kw::rotated_sod_shocktube::code() }
      , { ProblemType::SEDOV_BLASTWAVE, *kw::sedov_blastwave::code() }
      , { ProblemType::INTERFACE_ADVECTION, *kw::interface_advection::code() }
      , { ProblemType::GAUSS_HUMP_COMPFLOW, *kw::gauss_hump_compflow::code() }
      , { ProblemType::WATERAIR_SHOCKTUBE, *kw::waterair_shocktube::code() }
    };
};

} // ctr::
} // inciter::

#endif // ProblemOptions_h
