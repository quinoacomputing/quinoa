//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Geometry.h
  \author    J. Bakosi
  \date      Wed 06 Aug 2014 03:27:20 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Geometry options and associations
  \details   Geometry options and associations
*/
//******************************************************************************
#ifndef QuinoaGeometryOptions_h
#define QuinoaGeometryOptions_h

#include <map>

#include <Toggle.h>
#include <Quinoa/InputDeck/Keywords.h>

namespace quinoa {
namespace ctr {

//! Geometry definition types
enum class GeometryType : uint8_t { NO_GEOMETRY=0,
                                    ANALYTIC,
                                    DISCRETE };

//! Class with base templated on the above enum class with associations
class Geometry : public tk::Toggle< GeometryType > {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit Geometry() :
      Toggle< GeometryType >( "Geometry",
        //! Enums -> names
        { { GeometryType::NO_GEOMETRY, "n/a" },
          { GeometryType::ANALYTIC, kw::analytic_geometry().name() },
          { GeometryType::DISCRETE, kw::discrete_geometry().name() } },
        //! keywords -> Enums
        { { "no_geometry", GeometryType::NO_GEOMETRY },
          { kw::analytic_geometry().string(), GeometryType::ANALYTIC },
          { kw::discrete_geometry().string(), GeometryType::DISCRETE } } ) {}
};

} // ctr::
} // quinoa::

#endif // QuinoaGeometryOptions_h
