//******************************************************************************
/*!
  \file      src/Geometry/Primitive.h
  \author    J. Bakosi
  \date      Sat 24 Aug 2013 08:43:24 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Base for geometry primitives
  \details   Base for geometry primitives
*/
//******************************************************************************
#ifndef Primitive_h
#define Primitive_h

namespace Quinoa {

//! Primitive base
class Primitive {

  public:
    //! Constructor
    explicit Primitive() {}

    //! Destructor
    virtual ~Primitive() noexcept = default;

  private:
    //! Don't permit copy constructor
    Primitive(const Primitive&) = delete;
    //! Don't permit copy assigment
    Primitive& operator=(const Primitive&) = delete;
    //! Don't permit move constructor
    Primitive(Primitive&&) = delete;
    //! Don't permit move assigment
    Primitive& operator=(Primitive&&) = delete;
};

} // namespace Quinoa

#endif // Primitive_h
