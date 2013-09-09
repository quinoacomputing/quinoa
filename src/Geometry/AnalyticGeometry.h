//******************************************************************************
/*!
  \file      src/Geometry/AnalyticGeometry.h
  \author    J. Bakosi
  \date      Mon Sep  9 08:20:33 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Analytic geometry definition
  \details   Analytic geometry definition
*/
//******************************************************************************
#ifndef AnalyticGeometry_h
#define AnalyticGeometry_h

#include <vector>

#include <Exception.h>
#include <Geometry.h>

namespace quinoa {

class Primitive;

//! Analytic geometry definition
class AnalyticGeometry : public Geometry {

  public:
    //! Constructor
    explicit AnalyticGeometry( Memory* const memory,
                               Paradigm* const paradigm,
                               const QuinoaControl& control,
                               Timer* const timer);

    //! Destructor
    ~AnalyticGeometry() noexcept override = default;

    //! Initialize analytic geometry
    void init() override;

    //! Space-fill analytic geometry
    void fill() override;

  private:
    //! Don't permit copy constructor
    AnalyticGeometry(const AnalyticGeometry&) = delete;
    //! Don't permit copy assigment
    AnalyticGeometry& operator=(const AnalyticGeometry&) = delete;
    //! Don't permit move constructor
    AnalyticGeometry(AnalyticGeometry&&) = delete;
    //! Don't permit move assigment
    AnalyticGeometry& operator=(AnalyticGeometry&&) = delete;

    //! Echo information on analytic geometry
    void echo();

    std::vector<Primitive*> m_primitive;     //!< Vector of geometric primitives
};

} // namespace quinoa

#endif // AnalyticGeometry_h
