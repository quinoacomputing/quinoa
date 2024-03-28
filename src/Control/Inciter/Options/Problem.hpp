// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/Problem.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
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
                                   CYL_VORTEX,
                                   SHEDDING_FLOW,
                                   SOD_SHOCKTUBE,
                                   ROTATED_SOD_SHOCKTUBE,
                                   SEDOV_BLASTWAVE,
                                   INTERFACE_ADVECTION,
                                   GAUSS_HUMP_COMPFLOW,
                                   WATERAIR_SHOCKTUBE,
                                   SHOCK_HEBUBBLE,
                                   UNDERWATER_EX,
                                   SHOCKDENSITY_WAVE,
                                   EQUILINTERFACE_ADVECT,
                                   SINEWAVE_PACKET,
                                   RICHTMYER_MESHKOV };

//! Pack/Unpack ProblemType: forward overload to generic enum class packer
inline void operator|( PUP::er& p, ProblemType& e ) { PUP::pup( p, e ); }

//! \brief Problem options: outsource to base templated on enum type
class Problem : public tk::Toggle< ProblemType > {

  public:
    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit Problem() :
      tk::Toggle< ProblemType >(
        //! Group, i.e., options, name
        "Problem",
        //! Enums -> names
        { { ProblemType::USER_DEFINED, "user_defined" },
          { ProblemType::SHEAR_DIFF, "shear_diff" },
          { ProblemType::VORTICAL_FLOW, "vortical_flow" },
          { ProblemType::NL_ENERGY_GROWTH, "nl_energy_growth" },
          { ProblemType::RAYLEIGH_TAYLOR, "rayleigh_taylor" },
          { ProblemType::TAYLOR_GREEN, "taylor_green" },
          { ProblemType::SLOT_CYL, "slot_cyl" },
          { ProblemType::GAUSS_HUMP, "gauss_hump" },
          { ProblemType::CYL_ADVECT, "cyl_advect" },
          { ProblemType::CYL_VORTEX, "cyl_vortex" },
          { ProblemType::SHEDDING_FLOW, "shedding_flow" },
          { ProblemType::SOD_SHOCKTUBE, "sod_shocktube" },
          { ProblemType::ROTATED_SOD_SHOCKTUBE, "rotated_sod_shocktube" },
          { ProblemType::SEDOV_BLASTWAVE, "sedov_blastwave" },
          { ProblemType::INTERFACE_ADVECTION, "interface_advection" },
          { ProblemType::GAUSS_HUMP_COMPFLOW, "gauss_hump_compflow" },
          { ProblemType::WATERAIR_SHOCKTUBE, "waterair_shocktube" },
          { ProblemType::SHOCK_HEBUBBLE, "shock_hebubble" },
          { ProblemType::UNDERWATER_EX, "underwater_ex" },
          { ProblemType::SHOCKDENSITY_WAVE, "shockdensity_wave" },
          { ProblemType::EQUILINTERFACE_ADVECT, "equilinterface_advect" },
          { ProblemType::RICHTMYER_MESHKOV, "richtmyer_meshkov" },
          { ProblemType::SINEWAVE_PACKET, "sinewave_packet" }
        },
        //! keywords -> Enums
        { { "user_defined", ProblemType::USER_DEFINED },
          { "shear_diff", ProblemType::SHEAR_DIFF },
          { "vortical_flow", ProblemType::VORTICAL_FLOW },
          { "nl_energy_growth", ProblemType::NL_ENERGY_GROWTH },
          { "rayleigh_taylor", ProblemType::RAYLEIGH_TAYLOR },
          { "taylor_green", ProblemType::TAYLOR_GREEN },
          { "slot_cyl", ProblemType::SLOT_CYL },
          { "gauss_hump", ProblemType::GAUSS_HUMP },
          { "cyl_advect", ProblemType::CYL_ADVECT },
          { "cyl_vortex", ProblemType::CYL_VORTEX },
          { "shedding_flow", ProblemType::SHEDDING_FLOW },
          { "sod_shocktube", ProblemType::SOD_SHOCKTUBE },
          { "rotated_sod_shocktube", ProblemType::ROTATED_SOD_SHOCKTUBE },
          { "sod_shocktube", ProblemType::SOD_SHOCKTUBE },
          { "sedov_blastwave", ProblemType::SEDOV_BLASTWAVE },
          { "interface_advection", ProblemType::INTERFACE_ADVECTION },
          { "gauss_hump_compflow", ProblemType::GAUSS_HUMP_COMPFLOW },
          { "waterair_shocktube", ProblemType::WATERAIR_SHOCKTUBE },
          { "shock_hebubble", ProblemType::SHOCK_HEBUBBLE },
          { "underwater_ex", ProblemType::UNDERWATER_EX },
          { "shockdensity_wave", ProblemType::SHOCKDENSITY_WAVE },
          { "equilinterface_advect", ProblemType::EQUILINTERFACE_ADVECT },
          { "richtmyer_meshkov", ProblemType::RICHTMYER_MESHKOV },
          { "sinewave_packet", ProblemType::SINEWAVE_PACKET }
        } )
    {}
};

} // ctr::
} // inciter::

#endif // ProblemOptions_h
