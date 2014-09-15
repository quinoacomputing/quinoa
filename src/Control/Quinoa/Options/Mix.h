//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Mix.h
  \author    J. Bakosi
  \date      Mon 15 Sep 2014 08:14:07 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Mix model options and associations
  \details   Mix model options and associations
*/
//******************************************************************************
#ifndef QuinoaMixOptions_h
#define QuinoaMixOptions_h

#include <map>

#include <Toggle.h>
#include <Quinoa/InputDeck/Keywords.h>

namespace quinoa {
namespace ctr {

//! Mix model types
enum class MixType : uint8_t { NO_MIX=0,
                               IEM,
                               IECM,
                               DIRICHLET,
                               GENDIR };

//! Pack/Unpack BatteryType: forward overload to generic enum class packer
inline void operator|( PUP::er& p, MixType& e ) { tk::pup( p, e ); }

//! Class with base templated on the above enum class with associations
class Mix : public tk::Toggle< MixType > {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit Mix() :
      Toggle< MixType >( "Material mix",
        //! Enums -> names
        { { MixType::NO_MIX, "n/a" },
          { MixType::IEM, kw::mix_iem().name() },
          { MixType::IECM, kw::mix_iecm().name() },
          { MixType::DIRICHLET, kw::mix_dir().name() },
          { MixType::GENDIR, kw::mix_gendir().name() } },
       //! keywords -> Enums
       { { "no_mix", MixType::NO_MIX },
         { kw::mix_iem().string(), MixType::IEM },
         { kw::mix_iecm().string(), MixType::IECM },
         { kw::mix_dir().string(), MixType::DIRICHLET },
         { kw::mix_gendir().string(), MixType::GENDIR } } ) {}
};

} // ctr::
} // quinoa::

#endif // QuinoaMixOptions_h
