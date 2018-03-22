// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/AMRInitial.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Initial (before t=0) adaptive mesh refinement (AMR) options
  \details   Initial (before t=0) adaptive mesh refinement (AMR) options
*/
// *****************************************************************************
#ifndef InciterAMRInitialOptions_h
#define InciterAMRInitialOptions_h

#include "NoWarning/vector.h"

#include "Toggle.h"
#include "Keywords.h"
#include "PUPUtil.h"

namespace inciter {
namespace ctr {

//! Initial AMR types
enum class AMRInitialType : uint8_t { UNIFORM
                                    , INITIAL_CONDITIONS };

//! Pack/Unpack AMRInitialType: forward overload to generic enum class packer
inline void operator|( PUP::er& p, AMRInitialType& e )
{ PUP::pup( p, e ); }

//! AMRInitial options: outsource searches to base templated on enum type
class AMRInitial : public tk::Toggle< AMRInitialType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    using keywords = boost::mpl::vector< kw::amr_uniform
                                       , kw::amr_initial_conditions >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit AMRInitial() :
      tk::Toggle< AMRInitialType >(
        //! Group, i.e., options, name
        "Initial refinement typelist",
        //! Enums -> names
        { { AMRInitialType::UNIFORM, kw::amr_uniform::name() },
          { AMRInitialType::INITIAL_CONDITIONS,
            kw::amr_initial_conditions::name() } },
        //! keywords -> Enums
        { { kw::amr_uniform::string(), AMRInitialType::UNIFORM },
          { kw::amr_initial_conditions::string(),
            AMRInitialType::INITIAL_CONDITIONS } } )
    {
       boost::mpl::for_each< keywords >( assertPolicyCodes() );
    }

    //! \brief Return policy code based on Enum
    //! \param[in] p Enum value of the option requested
    //! \return Policy code of the option
    const std::string& code( AMRInitialType p ) const {
      using tk::operator<<;
      auto it = policy.find( p );
      Assert( it != end(policy),
              std::string("Cannot find policy code for initial AMR type \"")
              << p << "\"" );
      return it->second;
    }

  private:
    //! Function object for ensuring the existence of policy codes
    struct assertPolicyCodes {
      //! \brief Function call operator templated on the type to assert the
      //!   existence of a policy code
      template< typename U > void operator()( U ) {
        static_assert( tk::HasTypedefCode< typename U::info >::value,
                       "Policy code undefined for keyword" );
      }
    };

    //! Enums -> policy code
    std::map< AMRInitialType, std::string > policy {
        { AMRInitialType::UNIFORM, *kw::amr_uniform::code() }
      , { AMRInitialType::INITIAL_CONDITIONS,
          *kw::amr_initial_conditions::code() }
    };
};

} // ctr::
} // inciter::

#endif // InciterAMRInitialOptions_h
