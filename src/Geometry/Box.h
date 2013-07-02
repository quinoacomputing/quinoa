//******************************************************************************
/*!
  \file      src/Geometry/Box.h
  \author    J. Bakosi
  \date      Tue Jul  2 14:28:31 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Box primitive
  \details   Box primitive
*/
//******************************************************************************
#ifndef Box_h
#define Box_h

#include <QuinoaTypes.h>
#include <Primitive.h>
#include <Point.h>

namespace Quinoa {

//! Box primitive
class Box : public Primitive {

  public:
    //! Zero constructor
    explicit Box() : m_p1(0.0,0.0,0.0), m_p2(0.0,0.0,0.0) {}
    //! Fill constructor
    explicit Box(const Point& p1, const Point& p2) : m_p1(p1), m_p2(p2) {}

    //! Destructor
    virtual ~Box() noexcept {}

  private:
    //! Don't permit copy constructor
    Box(const Box&) = delete;
    //! Don't permit copy assigment
    Box& operator=(const Box&) = delete;
    //! Don't permit move constructor
    Box(Box&&) = delete;
    //! Don't permit move assigment
    Box& operator=(Box&&) = delete;

    Point m_p1;                 //!< Points determining the box in space
    Point m_p2;
};

} // namespace Quinoa

#endif // Box_h
