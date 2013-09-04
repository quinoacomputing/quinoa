//******************************************************************************
/*!
  \file      src/Model/Mix/Mix.h
  \author    J. Bakosi
  \date      Wed Sep  4 12:15:14 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mix model base
  \details   Mix mode lbase
*/
//******************************************************************************
#ifndef Mix_h
#define Mix_h

#include <cstring>

#include <QuinoaTypes.h>
#include <QuinoaControl.h>
#include <Model.h>
#include <MKLRandom.h>
#include <MKLRndStream.h>

namespace quinoa {

//! Mix model base for CRTP
//! See: http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
template< typename MixType >
class Mix : public Model {

  public:
    //! Constructor
    explicit Mix(Memory* const memory,
                 Paradigm* const paradigm,
                 const QuinoaControl& control,
                 real* const particles) :
      Model(memory,
            paradigm,
            control,
            particles,
            control.get<control::component, control::npar>(),
            control.nprop()),
      m_offset(control.scalarOffset()),
      m_nscalar(control.get<control::component, control::nscalar>()) {
      ErrChk(m_nscalar > 0, ExceptType::FATAL, "Wrong number of scalars");
    }

    //! Destructor
    ~Mix() noexcept override = default;

    //! CRTP interface: Initialize particles
    void init(int p, int tid) { static_cast<MixType*>(this)->init(p, tid); }

    //! CRTP interface: Advance particles in mix model
    void advance(int p, int tid, real dt) {
      static_cast<MixType*>(this)->advance(p, tid, dt);
    }

  protected:
    const int m_offset;             //!< Scalar-offset relative to base
    const int m_nscalar;            //!< Number of mixing scalars

    //! Initialize scalars with zero
    void initZero(int p) {
      memset(m_particles + p*m_nprop + m_offset, 0, m_nscalar*sizeof(real));
    }

  private:
    //! Don't permit copy constructor
    Mix(const Mix&) = delete;
    //! Don't permit copy assigment
    Mix& operator=(const Mix&) = delete;
    //! Don't permit move constructor
    Mix(Mix&&) = delete;
    //! Don't permit move assigment
    Mix& operator=(Mix&&) = delete;
};

} // namespace quinoa

#endif // Mix_h
