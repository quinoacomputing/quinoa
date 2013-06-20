//******************************************************************************
/*!
  \file      src/Control/GeometryOptions.h
  \author    J. Bakosi
  \date      Wed 19 Jun 2013 08:00:49 PM MDT
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
enum class GeometryTypes : uint8_t { NO_GEOMETRY=0,
                                     ANALYTIC,
                                     DISCRETE };

//! Class with base templated on the above enum class with associations
class Geometry : public Toggle<GeometryTypes> {

  public:
    //! Constructor initializing associations
    // ICC: use initializer lists
    explicit Geometry() : Toggle<GeometryTypes>(names, values) {
      //! Enums -> names
      names[GeometryTypes::NO_GEOMETRY] = "No geometry";
      names[GeometryTypes::ANALYTIC] = "Analytic";
      names[GeometryTypes::DISCRETE] = "Discrete";
      //! keywords -> Enums
      values["no_geometry"] = GeometryTypes::NO_GEOMETRY;
      values["analytic_geometry"] = GeometryTypes::ANALYTIC;
      values["discrete_geometry"] = GeometryTypes::DISCRETE;
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

    std::map<GeometryTypes, std::string> names;
    std::map<std::string, GeometryTypes> values;
};

} // namespace select

} // namespace Quinoa

#endif // GeometryOptions_h
