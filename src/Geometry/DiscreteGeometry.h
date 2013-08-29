//******************************************************************************
/*!
  \file      src/Geometry/DiscreteGeometry.h
  \author    J. Bakosi
  \date      Thu Aug 29 15:05:12 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Discrete geometry definition
  \details   Discrete geometry definition
*/
//******************************************************************************
#ifndef DiscreteGeometry_h
#define DiscreteGeometry_h

#include <Geometry.h>

namespace quinoa {

class STLMesh;
class QuinoaControl;

//! Discrete geometry definition
class DiscreteGeometry : public Geometry {

  public:
    //! Constructor
    explicit DiscreteGeometry(Memory* const memory,
                              Paradigm* const paradigm,
                              QuinoaControl* const control,
                              Timer* const timer);

    //! Destructor
    ~DiscreteGeometry() noexcept;

    //! Finalize, single exit point, called implicitly from destructor or
    //! explicitly from anywhere else
    virtual void finalize() noexcept;

    //! Initialize discrete geometry
    virtual void init();

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

    STLMesh* m_mesh;                        //!< Mesh object
};

} // namespace quinoa

#endif // DiscreteGeometry_h
