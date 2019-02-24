// *****************************************************************************
/*!
  \file      src/Control/Walker/Options/CoeffPolicy.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Differential equation coefficients policy options
  \details   Differential equation coefficients policy options
*/
// *****************************************************************************
#ifndef CoeffPolicyOptions_h
#define CoeffPolicyOptions_h

#include <brigand/sequences/list.hpp>
#include <brigand/algorithms/for_each.hpp>

#include "Toggle.h"
#include "Keywords.h"
#include "PUPUtil.h"

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
                                  , kw::instantaneous_velocity
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
            kw::instantaneous_velocity::name() } },
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
           { kw::instantaneous_velocity::string(),
             CoeffPolicyType::INSTANTANEOUS_VELOCITY } } )
    {
       brigand::for_each< keywords >( assertPolicyCodes() );
    }

    //! \brief Return policy code based on Enum
    //! \param[in] p Enum value of the option requested
    //! \return Policy code of the option
    const std::string& code( CoeffPolicyType p ) const {
      using tk::operator<<;
      auto it = policy.find( p );
      Assert( it != end(policy),
              std::string("Cannot find policy code for physics \"") << p <<
                "\"" );
      return it->second;
    }

  private:
    //! Function object for ensuring the existence of policy codes
    struct assertPolicyCodes {
      //! \brief Function call operator templated on the type to assert the
      //!   existence of a policy code
      template< typename U > void operator()( brigand::type_<U> ) {
        static_assert( tk::HasTypedefCode_v< typename U::info >,
                       "Policy code undefined for keyword" );
      }
    };

    //! Enums -> policy code
    std::map< CoeffPolicyType, std::string > policy {
        { CoeffPolicyType::CONST_COEFF, *kw::constcoeff::code() }
      , { CoeffPolicyType::DECAY, *kw::decay::code() }
      , { CoeffPolicyType::HOMOGENEOUS, *kw::homogeneous::code() }
      , { CoeffPolicyType::HOMOGENEOUS_DECAY, *kw::homdecay::code() }
      , { CoeffPolicyType::MONTE_CARLO_HOMOGENEOUS_DECAY,
          *kw::montecarlo_homdecay::code() }
      , { CoeffPolicyType::HYDROTIMESCALE, *kw::hydrotimescale::code() }
      , { CoeffPolicyType::CONST_SHEAR, *kw::const_shear::code() }
      , { CoeffPolicyType::STATIONARY, *kw::stationary::code() }
      , { CoeffPolicyType::INSTANTANEOUS_VELOCITY,
          *kw::instantaneous_velocity::code() }
    };

};

} // ctr::
} // walker::

#endif // CoeffPolicyOptions_h
