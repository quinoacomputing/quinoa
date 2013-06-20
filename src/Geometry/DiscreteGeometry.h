//******************************************************************************
/*!
  \file      src/Geometry/DiscreteGeometry.h
  \author    J. Bakosi
  \date      Thu Jun 20 06:49:03 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Discrete geometry definition
  \details   Discrete geometry definition
*/
//******************************************************************************
#ifndef DiscreteGeometry_h
#define DiscreteGeometry_h

#include <Exception.h>
#include <Geometry.h>

namespace Quinoa {

//! Discrete geometry definition
class DiscreteGeometry : public Geometry {

  public:
    //! Constructor
    explicit DiscreteGeometry( Memory* const memory,
                               Paradigm* const paradigm,
                               Control* const control,
                               Timer* const timer) :
      Geometry(memory, paradigm, control, timer) {}

    //! Destructor
    ~DiscreteGeometry() noexcept = default;

    //! Initialize analytic geometry
    virtual void init() {}

    //! Space-fill discrete geometry
    virtual void fill() {}

  private:
    //! Don't permit copy constructor
    DiscreteGeometry(const DiscreteGeometry&) = delete;
    //! Don't permit copy assigment
    DiscreteGeometry& operator=(const DiscreteGeometry&) = delete;
    //! Don't permit move constructor
    DiscreteGeometry(DiscreteGeometry&&) = delete;
    //! Don't permit move assigment
    DiscreteGeometry& operator=(DiscreteGeometry&&) = delete;
};

} // namespace Quinoa

#endif // DiscreteGeometry_h
