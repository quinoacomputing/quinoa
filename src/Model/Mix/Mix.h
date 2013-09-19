//******************************************************************************
/*!
  \file      src/Model/Mix/Mix.h
  \author    J. Bakosi
  \date      Thu Sep 19 09:14:08 2013
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

//! Mix model base
class Mix : public Model {

  public:
    //! Constructor
    explicit Mix(const Base& base, real* const particles) :
      Model(base,
            particles,
            base.control.get<ctr::component, ctr::npar>(),
            base.control.nprop()),
      m_offset(base.control.scalarOffset()),
      m_nscalar(base.control.get<ctr::component, ctr::nscalar>()) {
      ErrChk(m_nscalar > 0, ExceptType::FATAL, "Wrong number of scalars");
    }

    //! Destructor
    ~Mix() noexcept override = default;

    //! Initialize particles
    virtual void init(int p, int tid) = 0;

    //! Advance particles in mix model
    virtual void advance(int p, int tid, real dt) = 0;

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
