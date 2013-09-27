//******************************************************************************
/*!
  \file      src/Geometry/Geometry.h
  \author    J. Bakosi
  \date      Fri Sep 27 10:58:52 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Geometry base
  \details   Geometry base
*/
//******************************************************************************
#ifndef Geometry_h
#define Geometry_h

#include <Base.h>

namespace quinoa {

//! Geometry base
class Geometry {

  public:
    //! Constructor
    explicit Geometry(const Base& base);

    //! Destructor
    virtual ~Geometry() noexcept = default;

    //! Initialize geometry
    virtual void init() = 0;

    //! Space-fill geometry
    virtual void fill() = 0;

  protected:
    //! Echo information on geometry
    void echo();

    const Base& m_base;                      //!< Essentials

  private:
    //! Don't permit copy constructor
    Geometry(const Geometry&) = delete;
    //! Don't permit copy assigment
    Geometry& operator=(const Geometry&) = delete;
    //! Don't permit move constructor
    Geometry(Geometry&&) = delete;
    //! Don't permit move assigment
    Geometry& operator=(Geometry&&) = delete;
};

} // namespace quinoa

#endif // Geometry_h
