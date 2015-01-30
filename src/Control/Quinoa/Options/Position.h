//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Position.h
  \author    J. Bakosi
  \date      Wed 21 Jan 2015 08:00:34 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Position model options
  \details   Position model options
*/
//******************************************************************************
#ifndef QuinoaPositionOptions_h
#define QuinoaPositionOptions_h

#include <boost/mpl/vector.hpp>

#include <Toggle.h>
#include <Keywords.h>
#include <PUPUtil.h>

namespace quinoa {
namespace ctr {

//! Position model types
//! \author J. Bakosi
enum class PositionType : uint8_t { NO_POSITION=0,
                                    INVISCID,
                                    VISCOUS };

//! \brief Pack/Unpack PositionType: forward overload to generic enum class
//!    packer
//! \author J. Bakosi
inline void operator|( PUP::er& p, PositionType& e ) { PUP::pup( p, e ); }

//! \brief Position model options: outsource searches to base templated on enum
//!   type
//! \author J. Bakosi
class Position : public tk::Toggle< PositionType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    //! \author J. Bakosi
    using keywords = boost::mpl::vector< kw::pos_inviscid
                                       , kw::pos_viscous
                                       >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    //! \author J. Bakosi
    explicit Position() :
      Toggle< PositionType >(
        //! Group, i.e., options, name
        "Position",
        //! Enums -> names
        { { PositionType::NO_POSITION, "n/a" },
          { PositionType::INVISCID, kw::pos_inviscid::name() },
          { PositionType::VISCOUS, kw::pos_viscous::name() } },
        //! keywords -> Enums
        { { "no_position", PositionType::NO_POSITION },
          { kw::pos_inviscid::string(), PositionType::INVISCID },
          { kw::pos_viscous::string(), PositionType::INVISCID } } ) {}
};

} // ctr::
} // quinoa:::

#endif // QuinoaPositionOptions_h
