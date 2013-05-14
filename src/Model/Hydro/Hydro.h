//******************************************************************************
/*!
  \file      src/Model/Hydro/Hydro.h
  \author    J. Bakosi
  \date      Mon 13 May 2013 10:40:48 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Hydro base
  \details   Hydro base
*/
//******************************************************************************
#ifndef Hydro_h
#define Hydro_h

#include <QuinoaTypes.h>
#include <Model.h>
#include <Control.h>
#include <MKLRandom.h>
#include <MKLRndStream.h>

namespace Quinoa {

using namespace std;

//! Hydro model base for CRTP
//! See: http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
template< typename HydroType >
class Hydro : public Model {

  public:
    //! Constructor
    explicit Hydro(Memory* const memory,
                   Paradigm* const paradigm,
                   Control* const control,
                   int nvelocity,
                   real* const velocities)
      try :
        Model(memory, paradigm, control, control->get<control::NPAR>()),
        m_nvelocity(nvelocity),
        m_velocities(velocities)
      {

        ErrChk(m_nvelocity > 0, FATAL,
               "Wrong number of particle velocity components");
        Assert(m_velocities != nullptr, FATAL, "Velocity pointer null?");

        // Instantiate random number generator
        m_random = new (nothrow) MKLRandom(memory, paradigm);
        ErrChk(m_random != nullptr, FATAL,
               "Cannot allocate memory for random number generator");

        // Create random number leapfrog stream
        m_rndStr = m_random->addStream(VSL_BRNG_MCG59, 0);
        // Get array of MKL VSL stream state pointers right away
        m_str = m_random->getStr(m_rndStr);

      } // Roll back changes and rethrow on error
        catch (Exception& e) {
          // No need to clean up if exception thrown from base constructor
          if (e.func() == __PRETTY_FUNCTION__) finalize();
          throw;
        }
        catch (exception&) {
          finalize();
          throw;
        }
        catch (...) {
          finalize();
          Throw(UNCAUGHT, "Non-standard exception");
        }

    //! Destructor
    virtual ~Hydro() noexcept { finalize(); }

    //! CRTP interface: Initialize particles
    void init() { static_cast<HydroType*>(this)->init(); }

    //! CRTP interface: Advance particles in hydro model
    void advance(const real& dt) { static_cast<HydroType*>(this)->advance(dt); }

  protected:
    const int m_nvelocity;          //!< Number of velocities
    real* const m_velocities;       //!< Raw pointer to particle velocities
    const VSLStreamStatePtr* m_str; //!< Array of MKL VSL stream state pointers

    //! Constant accessor to random number stream object pointer
    //! \return Pointer to random number stream
    MKLRndStream* rndstr() const noexcept { return m_rndStr; }

  private:
    //! Don't permit copy constructor
    Hydro(const Hydro&) = delete;
    //! Don't permit copy assigment
    Hydro& operator=(const Hydro&) = delete;
    //! Don't permit move constructor
    Hydro(Hydro&&) = delete;
    //! Don't permit move assigment
    Hydro& operator=(Hydro&&) = delete;

    //! Finalize, single exit point, called implicitly from destructor or
    //! explicitly from anywhere else
    void finalize() noexcept {
      if (m_random) { delete m_random; m_random = nullptr; }  
    }

    MKLRandom* m_random;            //!< Random number generator object
    MKLRndStream* m_rndStr;         //!< Random number stream object
};

} // namespace Quinoa

#endif // Hydro_h
