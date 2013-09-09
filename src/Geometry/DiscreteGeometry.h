//******************************************************************************
/*!
  \file      src/Geometry/DiscreteGeometry.h
  \author    J. Bakosi
  \date      Mon Sep  9 08:19:05 2013
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
                              const QuinoaControl& control,
                              Timer* const timer);

    //! Destructor
    ~DiscreteGeometry() noexcept;

    //! Finalize, single exit point, called implicitly from destructor or
    //! explicitly from anywhere else
    void finalize() noexcept;

    //! Initialize discrete geometry
    void init() override;

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

    //! Echo information on discrete geometry
    void echo();

    STLMesh* m_mesh;                        //!< Mesh object
};

} // namespace quinoa

#endif // DiscreteGeometry_h
