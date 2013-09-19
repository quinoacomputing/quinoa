//******************************************************************************
/*!
  \file      src/Physics/Physics.h
  \author    J. Bakosi
  \date      Thu Sep 19 10:04:41 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Physics base
  \details   Physics base
*/
//******************************************************************************
#ifndef Physics_h
#define Physics_h

#include <QuinoaConfig.h>
#include <Base.h>
#include <Mass.h>
#include <Hydro.h>
#include <Mix.h>

namespace quinoa {

class Statistics;
class GlobWriter;
class TxtStatWriter;

//! Physics base
class Physics {

  public:
    //! Destructor
    virtual ~Physics() noexcept;

    //! Initialize physics
    virtual void init() = 0;

    //! Solve physics
    virtual void solve() = 0;

    //! Constant accessor to control object
    //! \return Control object
    const QuinoaControl& control() const noexcept { return m_base.control; }

    //! Constant accessor to timer object pointer
    //! \return Pointer to timer object
    const Timer& timer() const noexcept { return m_base.timer; }

    //! Constant accessor to mass model
    //! \return Pointer to mass model
    Mass* mass() const noexcept { return m_mass.get(); }

    //! Constant accessor to hydro model
    //! \return Pointer to hydro model
    Hydro* hydro() const noexcept { return m_hydro.get(); }

    //! Constant accessor to mix model
    //! \return Pointer to mix model
    Mix* mix() const noexcept { return m_mix.get(); }

    //! Constant accessor to statistics estimator
    //! \return Pointer to statistics estimator
    Statistics* statistics() const noexcept { return m_statistics; }

    //! Constant accessor to glob file writer
    //! \return Pointer to glob file writer
    GlobWriter* globWriter() const noexcept { return m_glob; }

    //! Constant accessor to statistics file writer
    //! \return Pointer to statistics file writer
    TxtStatWriter* statWriter() const noexcept { return m_stat; }

    //! Constant accessor to particle properties pointer
    //! \return Raw pointer to particle properties array
    real* particles() const noexcept { return m_particles.ptr; }

  protected:
    //! Constructor: protected, designed to be base-only
    explicit Physics(const Base& base);

    //! Echo information on physics
    void echo();

    const int m_nposition;                //!< Number of position components
    const int m_ndensity;                 //!< Number of density components
    const int m_nvelocity;                //!< Number of velocity components
    const int m_nscalar;                  //!< Number of scalar components
    const uint64_t m_npar;                //!< Number of particles
    const real m_term;                    //!< Maximum time to simulate

  private:
    //! Don't permit copy constructor
    Physics(const Physics&) = delete;
    //! Don't permit copy assigment
    Physics& operator=(const Physics&) = delete;
    //! Don't permit move constructor
    Physics(Physics&&) = delete;
    //! Don't permit move assigment
    Physics& operator=(Physics&&) = delete;

    //! Finalize, single exit point, called implicitly from destructor or
    //! explicitly from anywhere else
    void finalize() noexcept;

    //! Initialize model factory
    void initFactory();

    const Base& m_base;                   //!< Essentials

    //! Model factories
    std::map<sel::MassType, std::function<Mass*()>> m_massFactory;
    std::map<sel::HydroType, std::function<Hydro*()>> m_hydroFactory;
    std::map<sel::MixType, std::function<Mix*()>> m_mixFactory;

    std::unique_ptr<Mass> m_mass;         //!< Mass model
    std::unique_ptr<Hydro> m_hydro;       //!< Hydro model
    std::unique_ptr<Mix> m_mix;           //!< Mix model

    Statistics* m_statistics;             //!< Statistics estimator object
    GlobWriter* m_glob;                   //!< Glob file writer
    TxtStatWriter* m_stat;                //!< Statistics file writer
    Data<real> m_particles;               //!< Particle properties
};

} // namespace quinoa

#endif // Physics_h
