//******************************************************************************
/*!
  \file      src/Control/PositionOptions.h
  \author    J. Bakosi
  \date      Wed Sep 11 16:46:44 2013
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

namespace quinoa {

namespace select {

//! Position model types
enum class PositionType : uint8_t { NO_POSITION=0,
                                    INVISCID,
                                    VISCOUS };

//! Class with base templated on the above enum class with associations
class Position : public Toggle<PositionType> {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit Position() : Toggle<PositionType>(names, values) {}

  private:
    //! Don't permit copy constructor
    Position(const Position&) = delete;
    //! Don't permit copy assigment
    Position& operator=(const Position&) = delete;
    //! Don't permit move constructor
    Position(Position&&) = delete;
    //! Don't permit move assigment
    Position& operator=(Position&&) = delete;

    //! Enums -> names
    const std::map<PositionType, std::string> names {
      { PositionType::NO_POSITION, "" },
      { PositionType::INVISCID, "Inviscid" },
      { PositionType::VISCOUS, "Viscous" }
    };

    //! keywords -> Enums
    const std::map<std::string, PositionType> values {
      { "no_position", PositionType::NO_POSITION },
      { "pos_inviscid", PositionType::INVISCID },
      { "pos_viscous", PositionType::INVISCID }
    };
};

} // namespace select

} // namespace quinoa

#endif // PositionOptions_h
