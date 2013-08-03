//******************************************************************************
/*!
  \file      src/Control/PositionOptions.h
  \author    J. Bakosi
  \date      Fri Aug  2 15:42:47 2013
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
enum class PositionType : uint8_t { NO_POSITION=0,
                                    INVISCID,
                                    VISCOUS };

//! Class with base templated on the above enum class with associations
class Position : public Toggle<PositionType> {

  public:
    //! Constructor initializing associations
    // ICC: use initializer lists
    Position() : Toggle<PositionType>(names, values) {
      //! Enums -> names
      names[PositionType::NO_POSITION] = "No position";
      names[PositionType::INVISCID] = "Inviscid";
      names[PositionType::VISCOUS] = "Viscous";
      //! keywords -> Enums
      values["no_position"] = PositionType::NO_POSITION;
      values["pos_inviscid"] = PositionType::INVISCID;
      values["pos_viscous"] = PositionType::INVISCID;
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

    std::map<PositionType, std::string> names;
    std::map<std::string, PositionType> values;
};

} // namespace select

} // namespace Quinoa

#endif // PositionOptions_h
