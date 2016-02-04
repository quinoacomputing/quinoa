//******************************************************************************
/*!
  \file      src/Control/Walker/Options/InitPolicy.h
  \author    J. Bakosi
  \date      Thu 04 Feb 2016 06:14:29 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Differential equation initialization policy options for walker
  \details   Differential equation initialization policy options for walker
*/
//******************************************************************************
#ifndef InitPolicyOptions_h
#define InitPolicyOptions_h

#include <boost/mpl/vector.hpp>

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
      Toggle< InitPolicyType >(
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
          { kw::jointbeta::string(), InitPolicyType::JOINTBETA } } ) {}
};

} // ctr::
} // walker::

#endif // InitPolicyOptions_h
