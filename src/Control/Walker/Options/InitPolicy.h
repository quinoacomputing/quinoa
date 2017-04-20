// *****************************************************************************
/*!
  \file      src/Control/Walker/Options/InitPolicy.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Differential equation initialization policy options for walker
  \details   Differential equation initialization policy options for walker
*/
// *****************************************************************************
#ifndef InitPolicyOptions_h
#define InitPolicyOptions_h

#include <boost/mpl/vector.hpp>
#include "NoWarning/for_each.h"

#include "Toggle.h"
#include "Keywords.h"
#include "PUPUtil.h"

namespace walker {
namespace ctr {

//! Differential equation initializion policy types
//! \author J. Bakosi
enum class InitPolicyType : uint8_t { RAW=0,
                                      ZERO,
                                      JOINTDELTA,
                                      JOINTBETA };

//! Pack/Unpack InitPolicyType: forward overload to generic enum class packer
//! \author J. Bakosi
inline void operator|( PUP::er& p, InitPolicyType& e ) { PUP::pup( p, e ); }

//! InitPolicy options: outsource searches to base templated on enum type
//! \author J. Bakosi
class InitPolicy : public tk::Toggle< InitPolicyType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    //! \author J. Bakosi
    using keywords = boost::mpl::vector< kw::raw
                                       , kw::zero
                                       , kw::jointdelta
                                       , kw::jointbeta
                                       >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    //! \author J. Bakosi
    explicit InitPolicy() :
      tk::Toggle< InitPolicyType >(
        //! Group, i.e., options, name
        "Initialization Policy",
        //! Enums -> names
        { { InitPolicyType::RAW, kw::raw::name() },
          { InitPolicyType::ZERO, kw::zero::name() },
          { InitPolicyType::JOINTDELTA, kw::jointdelta::name() },
          { InitPolicyType::JOINTBETA, kw::jointbeta::name() } },
        //! keywords -> Enums
        { { kw::raw::string(), InitPolicyType::RAW },
          { kw::zero::string(), InitPolicyType::ZERO },
          { kw::jointdelta::string(), InitPolicyType::JOINTDELTA },
          { kw::jointbeta::string(), InitPolicyType::JOINTBETA } } )
    {
       boost::mpl::for_each< keywords >( assertPolicyCodes() );
    }

    //! \brief Return policy code based on Enum
    //! \param[in] p Enum value of the option requested
    //! \return Policy code of the option
    //! \author J. Bakosi
    const std::string& code( InitPolicyType p ) const {
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
    std::map< InitPolicyType, std::string > policy {
        { InitPolicyType::RAW, *kw::raw::code() }
      , { InitPolicyType::ZERO, *kw::zero::code() }
      , { InitPolicyType::JOINTDELTA, *kw::jointdelta::code() }
      , { InitPolicyType::JOINTBETA, *kw::jointbeta::code() }
    };
};

} // ctr::
} // walker::

#endif // InitPolicyOptions_h
