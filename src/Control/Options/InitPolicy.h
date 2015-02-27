//******************************************************************************
/*!
  \file      src/Control/Options/InitPolicy.h
  \author    J. Bakosi
  \date      Wed 25 Feb 2015 08:26:51 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Differential equation initialization policy options
  \details   Differential equation initialization policy options
*/
//******************************************************************************
#ifndef InitPolicyOptions_h
#define InitPolicyOptions_h

#include <boost/mpl/vector.hpp>

#include <Toggle.h>
#include <Keywords.h>
#include <PUPUtil.h>

namespace tk {
namespace ctr {

//! Differential equation initializion policy types
//! \author J. Bakosi
enum class InitPolicyType : uint8_t { RAW=0,
                                      ZERO,
                                      DELTA };

//! Pack/Unpack InitPolicyType: forward overload to generic enum class packer
//! \author J. Bakosi
inline void operator|( PUP::er& p, InitPolicyType& e ) { PUP::pup( p, e ); }

//! \brief InitPolicy options: outsourcs earches to base templated on enum type
//! \author J. Bakosi
class InitPolicy : public tk::Toggle< InitPolicyType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    //! \author J. Bakosi
    using keywords = boost::mpl::vector< kw::raw
                                       , kw::zero
                                       , kw::delta
                                       >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    //! \author J. Bakosi
    explicit InitPolicy() :
      Toggle< InitPolicyType >(
        //! Group, i.e., options, name
        "Initialization Policy",
        //! Enums -> names
        { { InitPolicyType::RAW, kw::raw::name() },
          { InitPolicyType::ZERO, kw::zero::name() },
          { InitPolicyType::DELTA, kw::delta::name() } },
        //! keywords -> Enums
        { { kw::raw::string(), InitPolicyType::RAW },
          { kw::zero::string(), InitPolicyType::ZERO },
          { kw::delta::string(), InitPolicyType::DELTA } } ) {}
};

} // ctr::
} // tk::

#endif // InitPolicyOptions_h
