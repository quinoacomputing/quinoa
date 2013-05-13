//******************************************************************************
/*!
  \file      src/Model/Mass/Mass.h
  \author    J. Bakosi
  \date      Sun 12 May 2013 09:08:30 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mass model base
  \details   Mass mode lbase
*/
//******************************************************************************
#ifndef Mass_h
#define Mass_h

#include <cstring>

#include <QuinoaTypes.h>
#include <Model.h>
#include <Control.h>
#include <MKLRandom.h>
#include <MKLRndStream.h>

namespace Quinoa {

using namespace std;

//! Mass model base for CRTP
//! See: http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
template< typename MassType >
class Mass : public Model {

  public:
    //! Constructor
    explicit Mass(Memory* const memory,
                  Paradigm* const paradigm,
                  Control* const control,
                  const int nscalar,
                  real* const densities)
      try :
        Model(memory, paradigm, control, control->get<control::NPAR>()),
        m_densities(densities),
        m_str(nullptr),
        m_random(nullptr),
        m_rndStr(nullptr)
      {

        Assert(m_densities != nullptr, FATAL, "Density pointer null?");

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
    virtual ~Mass() noexcept { finalize(); }

    //! CRTP interface: Initialize particles
    void init() { static_cast<MassType*>(this)->init(); }

    //! CRTP interface: Advance particles in mix model
    void advance(const real& dt) { static_cast<MassType*>(this)->advance(dt); }

  protected:
    real* const m_densities;        //!< Raw pointer to particle densities
    const VSLStreamStatePtr* m_str; //!< Array of MKL VSL stream state pointers

    //! Initialize densities with beta symmetric PDF
    //! \param[in] alpha  First shape parameter
    //! \param[in] beta   Second shape parameter
    //! \param[in] disp   Displacement (i.e., shift) parameter
    //! \param[in] scale  Scale parameter
    void initBeta(const real alpha,
                  const real beta,
                  const real disp,
                  const real scale) {
      for (int p=0; p<m_npar; ++p) {
        m_rndStr->beta(VSL_RNG_METHOD_BETA_CJA,
                       m_str[0], 1, m_densities+p, alpha, beta, disp, scale);
      }
    }

    //! Constant accessor to random number stream object pointer
    //! \return Pointer to random number stream
    MKLRndStream* rndstr() const noexcept { return m_rndStr; }

  private:
    //! Don't permit copy constructor
    Mass(const Mass&) = delete;
    //! Don't permit copy assigment
    Mass& operator=(const Mass&) = delete;
    //! Don't permit move constructor
    Mass(Mass&&) = delete;
    //! Don't permit move assigment
    Mass& operator=(Mass&&) = delete;

    //! Finalize, single exit point, called implicitly from destructor or
    //! explicitly from anywhere else
    void finalize() noexcept {
      if (m_random) { delete m_random; m_random = nullptr; }  
    }

    MKLRandom* m_random;            //!< Random number generator object
    MKLRndStream* m_rndStr;         //!< Random number stream object
};

} // namespace Quinoa

#endif // Mass_h
