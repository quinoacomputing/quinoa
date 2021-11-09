// *****************************************************************************
/*!
  \file      src/Control/Walker/Options/CoeffPolicy.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Differential equation coefficients policy options
  \details   Differential equation coefficients policy options
*/
// *****************************************************************************
#ifndef CoeffPolicyOptions_h
#define CoeffPolicyOptions_h

#include <brigand/sequences/list.hpp>
#include <brigand/algorithms/for_each.hpp>

#include "Toggle.hpp"
#include "Keywords.hpp"
#include "PUPUtil.hpp"

namespace walker {
namespace ctr {

//! Differential equation coefficients policy types
enum class CoeffPolicyType : uint8_t { CONST_COEFF=0
                                     , DECAY
                                     , HOMOGENEOUS
                                     , HOMOGENEOUS_DECAY
                                     , MONTE_CARLO_HOMOGENEOUS_DECAY
                                     , HYDROTIMESCALE
                                     , CONST_SHEAR
                                     , STATIONARY
                                     , INSTANTANEOUS_VELOCITY
                                     };

//! Pack/Unpack CoeffPolicyType: forward overload to generic enum class packer
inline void operator|( PUP::er& p, CoeffPolicyType& e ) { PUP::pup( p, e ); }

//! CoeffPolicy options: outsource searches to base templated on enum type
class CoeffPolicy : public tk::Toggle< CoeffPolicyType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    using keywords = brigand::list< kw::constcoeff
                                  , kw::decay
                                  , kw::homogeneous
                                  , kw::homdecay
                                  , kw::montecarlo_homdecay
                                  , kw::hydrotimescale
                                  , kw::const_shear
                                  , kw::stationary
                                  , kw::inst_velocity
                                  >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit CoeffPolicy() :
      tk::Toggle< CoeffPolicyType >(
        //! Group, i.e., options, name
        "Coefficients Policy",
        //! Enums -> names
        { { CoeffPolicyType::CONST_COEFF, kw::constcoeff::name() },
          { CoeffPolicyType::DECAY, kw::decay::name() },
          { CoeffPolicyType::HOMOGENEOUS, kw::homogeneous::name() },
          { CoeffPolicyType::HOMOGENEOUS_DECAY, kw::homdecay::name() },
          { CoeffPolicyType::MONTE_CARLO_HOMOGENEOUS_DECAY,
            kw::montecarlo_homdecay::name() },
          { CoeffPolicyType::HYDROTIMESCALE, kw::hydrotimescale::name() },
          { CoeffPolicyType::CONST_SHEAR, kw::const_shear::name() },
          { CoeffPolicyType::STATIONARY, kw::stationary::name() },
          { CoeffPolicyType::INSTANTANEOUS_VELOCITY,
            kw::inst_velocity::name() } },
        //! keywords -> Enums
        {  { kw::constcoeff::string(), CoeffPolicyType::CONST_COEFF },
           { kw::decay::string(), CoeffPolicyType::DECAY },
           { kw::homogeneous::string(), CoeffPolicyType::HOMOGENEOUS },
           { kw::homdecay::string(), CoeffPolicyType::HOMOGENEOUS_DECAY },
           { kw::montecarlo_homdecay::string(),
             CoeffPolicyType::MONTE_CARLO_HOMOGENEOUS_DECAY },
           { kw::hydrotimescale::string(), CoeffPolicyType::HYDROTIMESCALE },
           { kw::const_shear::string(), CoeffPolicyType::CONST_SHEAR },
           { kw::stationary::string(), CoeffPolicyType::STATIONARY },
           { kw::inst_velocity::string(),
             CoeffPolicyType::INSTANTANEOUS_VELOCITY } } )
    {}
};

} // ctr::
} // walker::

#endif // CoeffPolicyOptions_h
