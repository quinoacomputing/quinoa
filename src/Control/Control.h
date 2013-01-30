//******************************************************************************
/*!
  \file      src/Control/Control.h
  \author    J. Bakosi
  \date      Tue 29 Jan 2013 09:23:59 PM MST
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
                     DIRICLET,
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
    //! Get problem title
    const string& title() { return m_title; }

    //! Set physics
    void setPhysics(const string& physics);
    //! Get physics
    const PhysicsType& physics() { return m_physics; }

    //! Set hydrodynamics model
    void setHydro(const string& hydro);
    //! Get hydrodynamics model
    const HydroType& hydro() { return m_hydro; }

    //! Set material mix model
    void setMix(const string& mix);
    //! Get material mix model
    const MixType& mix() { return m_mix; }

    //! Set number of time steps to take
    void setNstep(const string& nstep);
    //! Get number of time steps to take
    const int& nstep() { return m_nstep; }

    //! Set value at which to stop simulation
    void setTerm(const string& term);
    //! Get value at which to stop simulation
    const real& term() { return m_term; }

    //! Set size of time step
    void setDt(const string& dt);
    //! Get size of time step
    const real& dt() { return m_dt; }

    //! Set number of mixing scalars
    void setNscalar(const string& nscalar);
    //! Get number of mixing scalars
    const int& nscalar() { return m_nscalar; }

    //! Set total number of particles
    void setNpar(const string& npar);
    //! Get total number of particles
    const int& npar() { return m_npar; }

    //! Set echo interval
    void setEcho(const string& echo);
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
