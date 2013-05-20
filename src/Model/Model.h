//******************************************************************************
/*!
  \file      src/Model/Model.h
  \author    J. Bakosi
  \date      Sun 19 May 2013 05:52:21 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Model base
  \details   Model base
*/
//******************************************************************************
#ifndef Model_h
#define Model_h

#include <QuinoaTypes.h>
#include <Exception.h>
#include <MKLRandom.h>
#include <MKLRndStream.h>

using namespace std;

namespace Quinoa {

class Memory;
class Paradigm;
class Control;

//! Model base
class Model {

  public:
    //! Constructor
    explicit Model(Memory* const memory,
                   Paradigm* const paradigm,
                   Control* const control,
                   real* const particles,
                   const int npar,
                   int nprop)
    try :
      m_memory(memory),
      m_paradigm(paradigm),
      m_control(control),
      m_particles(particles),
      m_npar(npar),
      m_nprop(nprop),
      m_str(nullptr),
      m_random(nullptr),
      m_rndStr(nullptr) {

      Assert(m_nprop != 0, FATAL, "Number of particle properties zero?");
      Assert(m_particles != nullptr, FATAL, "Particles pointer null?");
      ErrChk(m_npar > 0, FATAL, "Wrong number of particles");

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
    virtual ~Model() noexcept { finalize(); }

  protected:
    Memory* const m_memory;         //!< Memory object pointer
    Paradigm* const m_paradigm;     //!< Parallel programming object pointer
    Control* const m_control;       //!< Parallel programming object pointer
    real* const m_particles;        //!< Particles
    const int m_npar;               //!< Number of particles
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

} // namespace Quinoa

#endif // Model_h
