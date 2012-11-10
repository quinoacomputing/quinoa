//******************************************************************************
/*!
  \file      src/Random/MKL.h
  \author    J. Bakosi
  \date      Sat 10 Nov 2012 09:19:56 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKL-call wrappers with error handling
  \details   MKL-call wrappers with error handling
*/
//******************************************************************************
#ifndef MKL_h
#define MKL_h

#include <mkl_vsl.h>

#include <QuinoaTypes.h>

namespace Quinoa {

//! MKL-call wrappers with error handling
class MKL {

  public:
    //! Constructor: Default, compiler-generated
    MKL() = default;

    //! Destructor: Default, compiler-generated
    ~MKL() = default;

    //! Call MKL's vdRngUniform() and handle error
    void uniform(const int& method,
                 const VSLStreamStatePtr& stream,
                 const int& n,
                 real* r,
                 const real& a,
                 const real& b);

    //! Call MKL's vdRngGaussian() and handle error
    void gaussian(const int& method,
                  const VSLStreamStatePtr& stream,
                  const int& n,
                  real* r,
                  const real& a,
                  const real& b);

    //! Call MKL's vdRngGamma() and handle error
    void gamma(const int& method,
               const VSLStreamStatePtr& stream,
               const int& n,
               real* r,
               const real& alpha,
               const real& a,
               const real& beta);

  protected:
    //! Call MKL's vslNewStream() and handle error
    void newStream(VSLStreamStatePtr* stream,
                   const int& brng,
                   const unsigned int& seed);

    //! Call MKL's vslCopyStream() and handle error
    void copyStream(VSLStreamStatePtr* newstream,
                    const VSLStreamStatePtr& srcstream);

    //! Call MKL's vslSkipaheadStream() and handle error
    void skipAheadStream(VSLStreamStatePtr& stream,
                         const long long int& nskip);

    //! Call MKL's vslLeapfrogStream() and handle error
    void leapfrogStream(VSLStreamStatePtr& stream,
                        const int& k,
                        const int& nstreams);

  private:
    //! Don't permit copy constructor
    MKL(const MKL&) = delete;
    //! Don't permit copy assigment
    MKL& operator=(const MKL&) = delete;
    //! Don't permit move constructor
    MKL(MKL&&) = delete;
    //! Don't permit move assigment
    MKL& operator=(MKL&&) = delete;
};

} // namespace Quinoa

#endif // MKL_h
