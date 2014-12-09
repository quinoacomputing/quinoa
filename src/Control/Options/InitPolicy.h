//******************************************************************************
/*!
  \file      src/Control/Options/InitPolicy.h
  \author    J. Bakosi
  \date      Mon 08 Dec 2014 03:08:12 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Differential equation initialization policy options
  \details   Differential equation initialization policy options
*/
//******************************************************************************
#ifndef InitPolicyOptions_h
#define InitPolicyOptions_h

#include <map>

#include <Toggle.h>
#include <Keywords.h>

namespace tk {
namespace ctr {

//! Differential equation initializion policies
enum class InitPolicyType : uint8_t { RAW=0,
                                      ZERO };

//! Pack/Unpack: forward overload to generic enum class packer
inline void operator|( PUP::er& p, InitPolicyType& e ) { PUP::pup( p, e ); }

//! Class with base templated on the above enum class with associations
class InitPolicy : public tk::Toggle< InitPolicyType > {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit InitPolicy() :
      Toggle< InitPolicyType >( "Initialization Policy",
        //! Enums -> names
        { { InitPolicyType::RAW, kw::raw().name() },
          { InitPolicyType::ZERO, kw::zero().name() } },
        //! keywords -> Enums
        { { kw::raw().string(), InitPolicyType::RAW },
          { kw::zero().string(), InitPolicyType::ZERO } } ) {}
};

} // ctr::
} // tk::

#endif // InitPolicyOptions_h
