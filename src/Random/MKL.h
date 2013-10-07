//******************************************************************************
/*!
  \file      src/Random/MKL.h
  \author    J. Bakosi
  \date      Mon Oct  7 10:23:30 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKL-call wrappers with error handling
  \details   MKL-call wrappers with error handling
*/
//******************************************************************************
#ifndef MKL_h
#define MKL_h

#include <mkl_vsl.h>

#include <Types.h>

namespace tk {

//! MKL-call wrappers with error handling
class MKL {

  public:
    //! Constructor: Default, compiler-generated
    explicit MKL() noexcept = default;

    //! Destructor: Default, compiler-generated
    virtual ~MKL() noexcept = default;

    //! Call MKL's vdRngUniform() and handle error
    void uniform(const int& method,
                 const VSLStreamStatePtr& stream,
                 const int& n,
                 tk::real* r,
                 const tk::real& a,
                 const tk::real& b) const;

    //! Call MKL's vdRngGaussian() and handle error
    void gaussian(const int& method,
                  const VSLStreamStatePtr& stream,
                  const int& n,
                  tk::real* r,
                  const tk::real& a,
                  const tk::real& b) const;

    //! Call MKL's vdRngGamma() and handle error
    void gamma(const int& method,
               const VSLStreamStatePtr& stream,
               const int& n,
               tk::real* r,
               const tk::real& alpha,
               const tk::real& a,
               const tk::real& beta) const;

    //! Call MKL's vdRngBeta() and handle error
    void beta(const int& method,
              const VSLStreamStatePtr& stream,
              const int& n,
              tk::real* r,
              const tk::real& alpha,
              const tk::real& beta,
              const tk::real& disp,
              const tk::real& scale) const;

  protected:
    //! Call MKL's vslNewStream() and handle error
    void newStream(VSLStreamStatePtr* const stream,
                   const int& brng,
                   const unsigned int& seed) const;

    //! Call MKL's vslCopyStream() and handle error
    void copyStream(VSLStreamStatePtr* const newstream,
                    const VSLStreamStatePtr& srcstream) const;

    //! Call MKL's vslSkipaheadStream() and handle error
    void skipAheadStream(VSLStreamStatePtr& stream,
                         const long long int& nskip) const;

    //! Call MKL's vslLeapfrogStream() and handle error
    void leapfrogStream(VSLStreamStatePtr& stream,
                        const int& k,
                        const int& nstreams) const;

  private:
    //! Don't permit copy constructor
    MKL(const MKL&) = delete;
    //! Don't permit copy assigment
    MKL& operator=(const MKL&) = delete;
    //! Don't permit move constructor
    MKL(MKL&&) = delete;
    //! Don't permit move assigment
    MKL& operator=(MKL&&) = delete;

    //! Special error handler for MKL calls
    void MKLErrChk(int vslerr) const;
};

} // tk::

#endif // MKL_h
