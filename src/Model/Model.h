//******************************************************************************
/*!
  \file      src/Model/Model.h
  \author    J. Bakosi
  \date      Sun 15 Sep 2013 05:25:30 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Model base
  \details   Model base
*/
//******************************************************************************
#ifndef Model_h
#define Model_h

#include <QuinoaTypes.h>
#include <Base.h>
#include <Exception.h>
#include <MKLRandom.h>
#include <MKLRndStream.h>

namespace quinoa {

class Memory;
class Paradigm;
class QuinoaControl;

//! Model base
class Model {

  protected:
    //! Constructor: protected, designed to be base-only
    explicit Model(const Base& base,
                   real* const particles,
                   const uint64_t npar,
                   int nprop)
    try :
      m_base(base),
      m_particles(particles),
      m_npar(npar),
      m_nprop(nprop),
      m_str(nullptr),
      m_random(nullptr),
      m_rndStr(nullptr) {

      Assert(m_nprop != 0, ExceptType::FATAL,
             "Number of particle properties zero?");
      Assert(m_particles != nullptr, ExceptType::FATAL,
             "Particles pointer null?");
      ErrChk(m_npar > 0, ExceptType::FATAL,
             "Wrong number of particles");

      // Instantiate random number generator
      m_random = new (std::nothrow) MKLRandom(m_base);
      ErrChk(m_random != nullptr, ExceptType::FATAL,
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
      catch (std::exception&) {
        finalize();
        throw;
      }
      catch (...) {
        finalize();
        Throw(ExceptType::UNCAUGHT, "Non-standard exception");
      }

    //! Destructor: protected, designed to be freed via chlidren-only
    virtual ~Model() noexcept { finalize(); }

    const Base& m_base;             //!< Essentials
    real* const m_particles;        //!< Particles
    const uint64_t m_npar;          //!< Number of particles
    const int m_nprop;              //!< Number of particle properties
    const VSLStreamStatePtr* m_str; //!< Array of MKL VSL stream state pointers

    //! Initialize with uncorrelated joint Gaussian
    void initGaussian(int numvar,
                      MKLRndStream* const rndstr,
                      const VSLStreamStatePtr& str,
                      real mean = 0.0,
                      real rms = 1.0);

    //! Constant accessor to random number stream object pointer
    //! \return Pointer to random number stream
    MKLRndStream* rndstr() const noexcept { return m_rndStr; }

  private:
    //! Don't permit copy constructor
    Model(const Model&) = delete;
    //! Don't permit copy assigment
    Model& operator=(const Model&) = delete;
    //! Don't permit move constructor
    Model(Model&&) = delete;
    //! Don't permit move assigment
    Model& operator=(Model&&) = delete;

    //! Finalize, single exit point, called implicitly from destructor or
    //! explicitly from anywhere else
    void finalize() noexcept {
      if (m_random) { delete m_random; m_random = nullptr; }  
    }

    MKLRandom* m_random;            //!< Random number generator object
    MKLRndStream* m_rndStr;         //!< Random number stream object
};

} // namespace quinoa

#endif // Model_h
