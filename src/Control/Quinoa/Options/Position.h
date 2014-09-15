//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Position.h
  \author    J. Bakosi
  \date      Mon 15 Sep 2014 12:37:28 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Position model options and associations
  \details   Position model options and associations
*/
//******************************************************************************
#ifndef QuinoaPositionOptions_h
#define QuinoaPositionOptions_h

#include <map>

#include <Toggle.h>
#include <Quinoa/InputDeck/Keywords.h>

namespace quinoa {
namespace ctr {

//! Position model types
enum class PositionType : uint8_t { NO_POSITION=0,
                                    INVISCID,
                                    VISCOUS };

//! Pack/Unpack: forward overload to generic enum class packer
inline void operator|( PUP::er& p, PositionType& e ) { tk::pup( p, e ); }

//! Class with base templated on the above enum class with associations
class Position : public tk::Toggle< PositionType > {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit Position() :
      Toggle< PositionType >( "Position",
      //! Enums -> names
      { { PositionType::NO_POSITION, "n/a" },
        { PositionType::INVISCID, kw::pos_inviscid().name() },
        { PositionType::VISCOUS, kw::pos_viscous().name() } },
      //! keywords -> Enums
      { { "no_position", PositionType::NO_POSITION },
        { kw::pos_inviscid().string(), PositionType::INVISCID },
        { kw::pos_viscous().string(), PositionType::INVISCID } } ) {}
};

} // ctr::
} // quinoa:::

#endif // QuinoaPositionOptions_h
