//******************************************************************************
/*!
  \file      src/Model/Mix/Mix.h
  \author    J. Bakosi
  \date      Sun 12 May 2013 03:42:30 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mix model base
  \details   Mix mode lbase
*/
//******************************************************************************
#ifndef Mix_h
#define Mix_h

#include <cstring>

#include <QuinoaTypes.h>
#include <Model.h>
#include <Control.h>
#include <MKLRandom.h>
#include <MKLRndStream.h>

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
                 real* const scalars)
      try :
        Model(memory, paradigm, control, control->get<control::NPAR>()),
        m_nscalar(nscalar),
        m_scalars(scalars),
        m_str(nullptr),
        m_random(nullptr),
        m_rndStr(nullptr)
      {

        ErrChk(m_nscalar > 0, FATAL, "Wrong number of scalars");
        Assert(m_scalars != nullptr, FATAL, "Scalar pointer null?");

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
    virtual ~Mix() noexcept { finalize(); }

    //! CRTP interface: Initialize particles
    void init() { static_cast<MixType*>(this)->init(); }

    //! CRTP interface: Advance particles in mix model
    void advance(const real& dt) { static_cast<MixType*>(this)->advance(dt); }

  protected:
    const int m_nscalar;            //!< Number of mixing scalars
    real* const m_scalars;          //!< Raw pointer to particle scalars
    const VSLStreamStatePtr* m_str; //!< Array of MKL VSL stream state pointers

    //! Initialize scalars with uniform PDF with the last one constrained
    void initUniform() {
      real r[m_nscalar];

      // Generate initial values for all scalars for all particles
      for (int p=0; p<m_npar; ++p) {
        bool accept = false;
        while (!accept) {
          // Generate scalars
          m_rndStr->uniform(VSL_RNG_METHOD_UNIFORM_STD,
                            m_str[0], m_nscalar, r, 0.0, 1.0);

          // Compute their sum
          real sum = r[0];
          for (int i=1; i<m_nscalar; ++i) sum += r[i];

          // Accept if sum is less then 1.0
          if (sum < 1.0) {
            memcpy(m_scalars + p*m_nscalar, r, m_nscalar*sizeof(real));
            accept = true;
          }
        }
      }
    }

    //! Initialize scalars with Gaussian PDF
    void initGaussian() {
      real r[m_nscalar];

      // Generate initial values for all scalars for all particles
      for (int p=0; p<m_npar; ++p) {
        m_rndStr->gaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER,
                           m_str[0], m_nscalar, r, 0.0, 1.0);
        memcpy(m_scalars + p*m_nscalar, r, m_nscalar*sizeof(real));
      }
    }

    //! Constant accessor to random number stream object pointer
    //! \return Pointer to random number stream
    MKLRndStream* rndstr() const noexcept { return m_rndStr; }

  private:
    //! Don't permit copy constructor
    Mix(const Mix&) = delete;
    //! Don't permit copy assigment
    Mix& operator=(const Mix&) = delete;
    //! Don't permit move constructor
    Mix(Mix&&) = delete;
    //! Don't permit move assigment
    Mix& operator=(Mix&&) = delete;

    //! Finalize, single exit point, called implicitly from destructor or
    //! explicitly from anywhere else
    void finalize() noexcept {
      if (m_random) { delete m_random; m_random = nullptr; }  
    }

    MKLRandom* m_random;            //!< Random number generator object
    MKLRndStream* m_rndStr;         //!< Random number stream object
};

} // namespace Quinoa

#endif // Mix_h
