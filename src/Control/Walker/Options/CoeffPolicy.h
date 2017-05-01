// *****************************************************************************
/*!
  \file      src/Control/Walker/Options/CoeffPolicy.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Differential equation coefficients policy options
  \details   Differential equation coefficients policy options
*/
// *****************************************************************************
#ifndef CoeffPolicyOptions_h
#define CoeffPolicyOptions_h

#include <boost/mpl/vector.hpp>
#include "NoWarning/for_each.h"

#include "Toggle.h"
#include "Keywords.h"
#include "PUPUtil.h"

namespace walker {
namespace ctr {

//! Differential equation coefficients policy types
//! \author J. Bakosi
enum class CoeffPolicyType : uint8_t { CONSTANT=0
                                     , DECAY
                                     , HOMOGENEOUS_DECAY
                                     , MONTE_CARLO_HOMOGENEOUS_DECAY
                                     , HYDROTIMESCALE_HOMOGENEOUS_DECAY
                                     };

//! Pack/Unpack CoeffPolicyType: forward overload to generic enum class packer
//! \author J. Bakosi
inline void operator|( PUP::er& p, CoeffPolicyType& e ) { PUP::pup( p, e ); }

//! CoeffPolicy options: outsource searches to base templated on enum type
//! \author J. Bakosi
class CoeffPolicy : public tk::Toggle< CoeffPolicyType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    //! \author J. Bakosi
    using keywords = boost::mpl::vector< kw::constant
                                       , kw::decay
                                       , kw::homdecay
                                       , kw::montecarlo_homdecay
                                       , kw::hydrotimescale
                                       >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    //! \author J. Bakosi
    explicit CoeffPolicy() :
      tk::Toggle< CoeffPolicyType >(
        //! Group, i.e., options, name
        "Coefficients Policy",
        //! Enums -> names
        { { CoeffPolicyType::CONSTANT, kw::constant::name() },
          { CoeffPolicyType::DECAY, kw::decay::name() },
          { CoeffPolicyType::HOMOGENEOUS_DECAY, kw::homdecay::name() },
          { CoeffPolicyType::MONTE_CARLO_HOMOGENEOUS_DECAY,
            kw::montecarlo_homdecay::name() },
          { CoeffPolicyType::HYDROTIMESCALE_HOMOGENEOUS_DECAY,
            kw::hydrotimescale::name() } },
        //! keywords -> Enums
        {  { kw::constant::string(), CoeffPolicyType::CONSTANT },
           { kw::decay::string(), CoeffPolicyType::DECAY },
           { kw::homdecay::string(), CoeffPolicyType::HOMOGENEOUS_DECAY },
           { kw::montecarlo_homdecay::string(),
             CoeffPolicyType::MONTE_CARLO_HOMOGENEOUS_DECAY },
           { kw::hydrotimescale::string(),
             CoeffPolicyType::HYDROTIMESCALE_HOMOGENEOUS_DECAY } } )
    {
       boost::mpl::for_each< keywords >( assertPolicyCodes() );
    }

    //! \brief Return policy code based on Enum
    //! \param[in] p Enum value of the option requested
    //! \return Policy code of the option
    //! \author J. Bakosi
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
    //! \author J. Bakosi
    struct assertPolicyCodes {
      //! \brief Function call operator templated on the type to assert the
      //!   existence of a policy code
      template< typename U > void operator()( U ) {
        static_assert( tk::HasTypedefCode< typename U::info >::value,
                       "Policy code undefined for keyword" );
      }
    };

    //! Enums -> policy code
    std::map< CoeffPolicyType, std::string > policy {
        { CoeffPolicyType::CONSTANT, *kw::constant::code() }
      , { CoeffPolicyType::DECAY, *kw::decay::code() }
      , { CoeffPolicyType::HOMOGENEOUS_DECAY, *kw::homdecay::code() }
      , { CoeffPolicyType::MONTE_CARLO_HOMOGENEOUS_DECAY,
          *kw::montecarlo_homdecay::code() }
      , { CoeffPolicyType::HYDROTIMESCALE_HOMOGENEOUS_DECAY,
          *kw::hydrotimescale::code() }
    };

};

} // ctr::
} // walker::

#endif // CoeffPolicyOptions_h
