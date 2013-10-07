//******************************************************************************
/*!
  \file      src/Geometry/DiscreteGeometry.h
  \author    J. Bakosi
  \date      Mon Oct  7 08:22:22 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Discrete geometry definition
  \details   Discrete geometry definition
*/
//******************************************************************************
#ifndef DiscreteGeometry_h
#define DiscreteGeometry_h

#include <Geometry.h>
#include <STLMesh.h>

namespace quinoa {

//! Discrete geometry definition
class DiscreteGeometry : public Geometry {

  public:
    //! Constructor
    explicit DiscreteGeometry(const Base& base);

    //! Destructor
    ~DiscreteGeometry() noexcept override {}

    //! Initialize discrete geometry
    void init() override {}

    //! Space-fill discrete geometry
    void fill() override {}

  private:
    //! Don't permit copy constructor
    DiscreteGeometry(const DiscreteGeometry&) = delete;
    //! Don't permit copy assigment
    DiscreteGeometry& operator=(const DiscreteGeometry&) = delete;
    //! Don't permit move constructor
    DiscreteGeometry(DiscreteGeometry&&) = delete;
    //! Don't permit move assigment
    DiscreteGeometry& operator=(DiscreteGeometry&&) = delete;

    STLMesh m_mesh;             //!< Mesh object
};

} // quinoa::

#endif // DiscreteGeometry_h
