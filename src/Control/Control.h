//******************************************************************************
/*!
  \file      src/Control/Control.h
  \author    J. Bakosi
  \date      Sat 26 Jan 2013 09:28:47 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Main control category
  \details   Main control catgeory
*/
//******************************************************************************
#ifndef Control_h
#define Control_h

namespace Quinoa {

//! Control base
class Control {

  public:
    //! Constructor
    Control() = default;

    //! Destructor
    ~Control() = default;

    //! Set title
    void setTitle(const string& title) { m_title = title; }

  private:
    //! Don't permit copy constructor
    Control(const Control&) = delete;
    //! Don't permit copy assigment
    Control& operator=(const Control&) = delete;
    //! Don't permit move constructor
    Control(Control&&) = delete;
    //! Don't permit move assigment
    Control& operator=(Control&&) = delete;

    string m_title;               //!< Title
};


//! Available physics (methods) options
enum class PhysicsType {
  HOMOGENEOUS_DIRICHLET,     //!< Homogeneous Dirichlet
  HOMOGENEOUS_GENDIRICHLET,  //!< Homogeneous Generalized Dirichlet
  SPINSFLOW           //!< Standalone-Particle Incompressible Navier-Stokes Flow
};

//! Available hydrodynamics models
enum class HydroType {
  SLM,     //!< Simplified Langevin model
  GLM      //!< Generalized Langevin model
};

} // namespace Quinoa

#endif // Control_h
