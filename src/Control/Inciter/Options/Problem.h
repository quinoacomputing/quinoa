// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/Problem.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Problem options for inciter
  \details   Problem options for inciter
*/
// *****************************************************************************
#ifndef ProblemOptions_h
#define ProblemOptions_h

#include <boost/mpl/vector.hpp>
#include "NoWarning/for_each.h"

#include "Toggle.h"
#include "Keywords.h"
#include "PUPUtil.h"

namespace inciter {
namespace ctr {

//! Problem types
//! \author J. Bakosi
enum class ProblemType : uint8_t { USER_DEFINED=0,
                                   SHEAR_DIFF,
                                   DIR_NEU,
                                   VORTICAL_FLOW,
                                   SLOT_CYL };

//! Pack/Unpack ProblemType: forward overload to generic enum class packer
//! \author J. Bakosi
inline void operator|( PUP::er& p, ProblemType& e ) { PUP::pup( p, e ); }

//! \brief Problem options: outsource to base templated on enum type
//! \author J. Bakosi
class Problem : public tk::Toggle< ProblemType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    //! \author J. Bakosi
    using keywords = boost::mpl::vector< kw::user_defined
                                       , kw::shear_diff
                                       , kw::dir_neu
                                       , kw::vortical_flow
                                       , kw::slot_cyl
                                       >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    //! \author J. Bakosi
    explicit Problem() :
      tk::Toggle< ProblemType >(
        //! Group, i.e., options, name
        "Test problem",
        //! Enums -> names
        { { ProblemType::USER_DEFINED, kw::user_defined::name() },
          { ProblemType::SHEAR_DIFF, kw::shear_diff::name() },
          { ProblemType::DIR_NEU, kw::dir_neu::name() },
          { ProblemType::VORTICAL_FLOW, kw::vortical_flow::name() },
          { ProblemType::SLOT_CYL, kw::slot_cyl::name() } },
        //! keywords -> Enums
        { { kw::user_defined::string(), ProblemType::USER_DEFINED },
          { kw::shear_diff::string(), ProblemType::SHEAR_DIFF },
          { kw::dir_neu::string(), ProblemType::DIR_NEU },
          { kw::vortical_flow::string(), ProblemType::VORTICAL_FLOW },
          { kw::slot_cyl::string(), ProblemType::SLOT_CYL } } )
    {
       boost::mpl::for_each< keywords >( assertPolicyCodes() );
    }

    //! \brief Return policy code based on Enum
    //! \param[in] p Enum value of the problem option requested
    //! \return Policy code of the option
    //! \author J. Bakosi
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
    std::map< ProblemType, std::string > policy {
        { ProblemType::USER_DEFINED, *kw::user_defined::code() }
      , { ProblemType::SHEAR_DIFF, *kw::shear_diff::code() }
      , { ProblemType::DIR_NEU, *kw::dir_neu::code() }
      , { ProblemType::VORTICAL_FLOW, *kw::vortical_flow::code() }
      , { ProblemType::SLOT_CYL, *kw::slot_cyl::code() }
    };
};

} // ctr::
} // inciter::

#endif // ProblemOptions_h
