//******************************************************************************
/*!
  \file      src/Model/Mix/Mix.h
  \author    J. Bakosi
  \date      Sun 12 May 2013 06:56:25 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mix model base
  \details   Mix mode lbase
*/
//******************************************************************************
#ifndef Mix_h
#define Mix_h

#include <QuinoaTypes.h>
#include <Model.h>
#include <Control.h>

namespace Quinoa {

using namespace std;

//! Mix model base for CRTP
//! See: http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
template< typename MixType >
class Mix : public Model {

  public:
    //! Constructor
    explicit Mix(Memory* const memory,
                 Paradigm* const paradigm,
                 Control* const control,
                 const int nscalar,
                 real* const scalars) :
      Model(memory, paradigm, control, control->get<control::NPAR>()),
      m_nscalar(nscalar),
      m_scalars(scalars) {
      ErrChk(m_nscalar > 0, FATAL, "Wrong number of scalars");      
      Assert(m_scalars != nullptr, FATAL, "Scalar pointer null?");
    }

    //! Destructor
    virtual ~Mix() noexcept = default;

    //! Initialize particles
    void init() { static_cast<MixType*>(this)->init(); }

    //! Advance particles in mix model
    void advance(const real& dt) { static_cast<MixType*>(this)->advance(dt); }

    //! Echo information on mix model
    void echo() { static_cast<MixType*>(this)->echo(); }

  protected:
    const int m_nscalar;            //!< Number of mixing scalars
    real* const m_scalars;          //!< Raw pointer to particle scalars

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

} // namespace Quinoa

#endif // Mix_h
