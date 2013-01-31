//******************************************************************************
/*!
  \file      src/Control/Control.h
  \author    J. Bakosi
  \date      Thu 31 Jan 2013 06:57:00 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Main control category
  \details   Main control catgeory
*/
//******************************************************************************
#ifndef Control_h
#define Control_h

#include<string>

#include<QuinoaTypes.h>

using namespace std;

namespace Quinoa {

//! Physics (methods: collection of models) types
enum class PhysicsType { NO_PHYSICS,
                         HOMOGENEOUS_DIRICHLET,
                         HOMOGENEOUS_GENERALIZED_DIRICHLET,
                         SPINSFLOW
};

//! Hydrodynamics model types
enum class HydroType { NO_HYDRO,
                       SLM,
                       GLM
};

//! Material mix model types
enum class MixType { NO_MIX,
                     IEM,
                     IECM,
                     DIRICHLET,
                     GENERALIZED_DIRICHLET
};

//! Control base
class Control {

  public:
    //! Constructor
    Control();

    //! Destructor
    ~Control() = default;

    //! Set problem title
    void setTitle(const string& title) { m_title = title; }
    //! Set physics
    void setPhysics(const PhysicsType& physics) { m_physics = physics; }
    //! Set hydrodynamics model
    void setHydro(const HydroType& hydro) { m_hydro = hydro; }
    //! Set material mix model
    void setMix(const MixType& mix);
    //! Set number of time steps to take
    void setNstep(const int& nstep);
    //! Set value at which to stop simulation
    void setTerm(const real& term);
    //! Set size of time step
    void setDt(const real& dt);
    //! Set number of mixing scalars
    void setNscalar(const int& nscalar);
    //! Set total number of particles
    void setNpar(const int& npar);
    //! Set echo interval
    void setEcho(const int& echo);

    //! Get problem title
    const string& title() { return m_title; }
    //! Get physics
    const PhysicsType& physics() { return m_physics; }
    //! Get hydrodynamics model
    const HydroType& hydro() { return m_hydro; }
    //! Get material mix model
    const MixType& mix() { return m_mix; }
    //! Get number of time steps to take
    const int& nstep() { return m_nstep; }
    //! Get value at which to stop simulation
    const real& term() { return m_term; }
    //! Get size of time step
    const real& dt() { return m_dt; }
    //! Get number of mixing scalars
    const int& nscalar() { return m_nscalar; }
    //! Get total number of particles
    const int& npar() { return m_npar; }
    //! Get echo interval
    const int& echo() { return m_echo; }

  private:
    //! Don't permit copy constructor
    Control(const Control&) = delete;
    //! Don't permit copy assigment
    Control& operator=(const Control&) = delete;
    //! Don't permit move constructor
    Control(Control&&) = delete;
    //! Don't permit move assigment
    Control& operator=(Control&&) = delete;

    //! The parser stores everything the user selected in these variables
    string m_title;               //!< Title
    PhysicsType m_physics;        //!< Physics
    HydroType m_hydro;            //!< Hydrodynamics model
    MixType m_mix;                //!< Material mix model
    int m_nstep;                  //!< Number of time steps to take
    real m_term;                  //!< Terminate time stepping at this time
    real m_dt;                    //!< Size of time step
    int m_nscalar;                //!< Number of mixing scalars
    int m_npar;                   //!< Total number of particles
    int m_echo;                   //!< One-line info every few time steps
};

} // namespace Quinoa

#endif // Control_h
