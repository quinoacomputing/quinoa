//******************************************************************************
/*!
  \file      src/Control/PositionOptions.h
  \author    J. Bakosi
  \date      Mon 27 May 2013 07:12:28 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Position model options and associations
  \details   Position model options and associations
*/
//******************************************************************************
#ifndef PositionOptions_h
#define PositionOptions_h

#include <map>

#include <Exception.h>
#include <Toggle.h>

namespace Quinoa {

namespace select {

//! Position model types
enum class PositionTypes { NO_POSITION=0,
                           INVISCID,
                           VISCOUS };

//! Class with base templated on the above enum class with associations
class Position : public Toggle<PositionTypes> {

  public:
    //! Constructor initializing associations
    // ICC: use initializer lists
    Position() : Toggle<PositionTypes>(names, values) {
      //! Enums -> names
      names[PositionTypes::NO_POSITION] = "No position";
      names[PositionTypes::INVISCID] = "Inviscid";
      names[PositionTypes::VISCOUS] = "Viscous";
      //! keywords -> Enums
      values["no_position"] = PositionTypes::NO_POSITION;
      values["invpos"] = PositionTypes::INVISCID;
      values["vispos"] = PositionTypes::INVISCID;
    }

  private:
    //! Don't permit copy constructor
    Position(const Position&) = delete;
    //! Don't permit copy assigment
    Position& operator=(const Position&) = delete;
    //! Don't permit move constructor
    Position(Position&&) = delete;
    //! Don't permit move assigment
    Position& operator=(Position&&) = delete;

    std::map<PositionTypes, std::string> names;
    std::map<std::string, PositionTypes> values;
};

} // namespace select

} // namespace Quinoa

#endif // PositionOptions_h
