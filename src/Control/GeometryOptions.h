//******************************************************************************
/*!
  \file      src/Control/GeometryOptions.h
  \author    J. Bakosi
  \date      Fri Aug  2 15:42:06 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Geometry options and associations
  \details   Geometry options and associations
*/
//******************************************************************************
#ifndef GeometryOptions_h
#define GeometryOptions_h

#include <map>

#include <Exception.h>
#include <Toggle.h>

namespace Quinoa {

namespace select {

//! Geometry definition types
enum class GeometryType : uint8_t { NO_GEOMETRY=0,
                                    ANALYTIC,
                                    DISCRETE };

//! Class with base templated on the above enum class with associations
class Geometry : public Toggle<GeometryType> {

  public:
    //! Constructor initializing associations
    // ICC: use initializer lists
    explicit Geometry() : Toggle<GeometryType>(names, values) {
      //! Enums -> names
      names[GeometryType::NO_GEOMETRY] = "No geometry";
      names[GeometryType::ANALYTIC] = "Analytic";
      names[GeometryType::DISCRETE] = "Discrete";
      //! keywords -> Enums
      values["no_geometry"] = GeometryType::NO_GEOMETRY;
      values["analytic_geometry"] = GeometryType::ANALYTIC;
      values["discrete_geometry"] = GeometryType::DISCRETE;
    }

  private:
    //! Don't permit copy constructor
    Geometry(const Geometry&) = delete;
    //! Don't permit copy assigment
    Geometry& operator=(const Geometry&) = delete;
    //! Don't permit move constructor
    Geometry(Geometry&&) = delete;
    //! Don't permit move assigment
    Geometry& operator=(Geometry&&) = delete;

    std::map<GeometryType, std::string> names;
    std::map<std::string, GeometryType> values;
};

} // namespace select

} // namespace Quinoa

#endif // GeometryOptions_h
