//******************************************************************************
/*!
  \file      src/Random/MKL.C
  \author    J. Bakosi
  \date      Thu Aug 29 15:35:31 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKL-call wrappers with error handling
  \details   MKL-call wrappers with error handling
*/
//******************************************************************************

#include <sstream>
#include <iostream>

#include <MKL.h>
#include <Exception.h>

using namespace quinoa;

void
MKL::uniform(const int& method,
             const VSLStreamStatePtr& stream,
             const int& n,
             real* r,
             const real& a,
             const real& b) const
//******************************************************************************
//  Call MKL's vdRngUniform() and handle error
//! \param[in]  method   Generation method
//! \param[in]  stream   Pointer to the stream state structure
//! \param[in]  n        Number of random values to be generated
//! \param[out] r        Vector of n uniform random numbers
//! \param[in]  a        Left bound
//! \param[in]  b        Right bound
//! \author  J. Bakosi
//******************************************************************************
{
#ifdef MKL_CALLS
#ifdef NDEBUG
  vdRngUniform(method, stream, n, r, a, b);
#else  // NDEBUG
  MKLErrChk(vdRngUniform(method, stream, n, r, a, b));
#endif // NDEBUG
#endif
}

void
MKL::gaussian(const int& method,
              const VSLStreamStatePtr& stream,
              const int& n,
              real* r,
              const real& a,
              const real& b) const
//******************************************************************************
//  Call MKL's vdRngGaussian() and handle error
//! \param[in]  method   Generation method
//! \param[in]  stream   Pointer to the stream state structure
//! \param[in]  n        Number of random values to be generated
//! \param[out] r        Vector of n Gaussian random numbers
//! \param[in]  a        mean
//! \param[in]  b        standard deviation
//! \author  J. Bakosi
//******************************************************************************
{
#ifdef MKL_CALLS
#ifdef NDEBUG
  vdRngGaussian(method, stream, n, r, a, b);
#else  // NDEBUG
  MKLErrChk(vdRngGaussian(method, stream, n, r, a, b));
#endif // NDEBUG
#endif
}

void
MKL::gamma(const int& method,
           const VSLStreamStatePtr& stream,
           const int& n,
           real* r,
           const real& alpha,
           const real& a,
           const real& beta) const
//******************************************************************************
//  Call MKL's vdRngGamma() and handle error
//! \param[in]  method   Generation method
//! \param[in]  stream   Pointer to the stream state structure
//! \param[in]  n        Number of random values to be generated
//! \param[out] r        Vector of n gamma-distributed random numbers
//! \param[in]  alpha    Shape
//! \param[in]  a        Displacement
//! \param[in]  beta     Scale factor
//! \author  J. Bakosi
//******************************************************************************
{
#ifdef MKL_CALLS
#ifdef NDEBUG
  vdRngGamma(method, stream, n, r, alpha, a, beta);
#else  // NDEBUG
  MKLErrChk(vdRngGamma(method, stream, n, r, alpha, a, beta));
#endif // NDEBUG
#endif
}

void
MKL::beta(const int& method,
          const VSLStreamStatePtr& stream,
          const int& n,
          real* r,
          const real& alpha,
          const real& beta,
          const real& disp,
          const real& scale) const
//******************************************************************************
//  Call MKL's vdRngBeta() and handle error
//! \param[in]  method   Generation method
//! \param[in]  stream   Pointer to the stream state structure
//! \param[in]  n        Number of random values to be generated
//! \param[out] r        Vector of n beta-distributed random numbers
//! \param[in] alpha     First shape parameter
//! \param[in] beta      Second shape parameter
//! \param[in] disp      Displacement (i.e., shift) parameter
//! \param[in] scale     Scale parameter
//! \author  J. Bakosi
//******************************************************************************
{
#ifdef MKL_CALLS
#ifdef NDEBUG
  vdRngBeta(method, stream, n, r, alpha, beta, disp, scale);
#else  // NDEBUG
  MKLErrChk(vdRngBeta(method, stream, n, r, alpha, beta, disp, scale));
#endif // NDEBUG
#endif
}

void
MKL::newStream(VSLStreamStatePtr* const stream,
               const int& brng,
               const unsigned int& seed) const
//******************************************************************************
//  Call MKL's vslNewStream() and handle error
//! \param[out]  stream  VSL stream state descriptor
//! \param[in]   brng    Index of the basic generator to initialize the stream
//! \param[in]   seed    Initial condition of the stream
//! \author  J. Bakosi
//******************************************************************************
{
#ifdef MKL_CALLS
#ifdef NDEBUG
  vslNewStream(stream, brng, seed);
#else  // NDEBUG
  MKLErrChk(vslNewStream(stream, brng, seed));
#endif // NDEBUG
#endif
}

void
MKL::copyStream(VSLStreamStatePtr* const newstream,
                const VSLStreamStatePtr& srcstream) const
//******************************************************************************
//  Call MKL's vslCopyStream() and handle error
//! \param[out]  newstream  Copied random stream descriptor
//! \param[in]   srcstream  Pointer to the stream state structure to be copied
//! \author  J. Bakosi
//******************************************************************************
{
#ifdef MKL_CALLS
#ifdef NDEBUG
  vslCopyStream(newstream, srcstream);
#else  // NDEBUG
  MKLErrChk(vslCopyStream(newstream, srcstream));
#endif // NDEBUG
#endif
}

void
MKL::skipAheadStream(VSLStreamStatePtr& stream,
                     const long long int& nskip) const
//******************************************************************************
//  Call MKL's vslCopyStream() and handle error
//! \param[in]  stream  Pointer to the stream state structure to which the
//!                     block-splitting method is applied
//! \param[in]  nskip   Number of skipped elements
//! \author  J. Bakosi
//******************************************************************************
{
#ifdef MKL_CALLS
#ifdef NDEBUG
  vslSkipAheadStream(stream, nskip);
#else  // NDEBUG
  MKLErrChk(vslSkipAheadStream(stream, nskip));
#endif // NDEBUG
#endif
}

void
MKL::leapfrogStream(VSLStreamStatePtr& stream,
                    const int& k,
                    const int& nstreams) const
//******************************************************************************
//  Call MKL's vslCopyStream() and handle error
//! \param[in]  stream   Pointer to the stream state structure to which leapfrog
//                       method is applied
//! \param[in]  k        Index of the computational node, or stream number
//! \param[in]  nstreams Largest number of computational nodes, or stride
//! \author  J. Bakosi
//******************************************************************************
{
#ifdef MKL_CALLS
#ifdef NDEBUG
  vslLeapfrogStream(stream, k, nstreams);
#else  // NDEBUG
  MKLErrChk(vslLeapfrogStream(stream, k, nstreams));
#endif // NDEBUG
#endif
}

void
MKL::MKLErrChk(int vslerr) const
//******************************************************************************
//  Special error handler for MKL
//! \param[in]  vslerr     Error code
//! \author  J. Bakosi
//******************************************************************************
{
  if (vslerr != VSL_STATUS_OK)
    try {

      std::stringstream s;
      s << "MKL VSL Error: code " << vslerr;
      Throw(ExceptType::FATAL, s.str());

    } catch (Exception&) {
        throw;
      }
      catch (std::exception& e) {
        Throw(ExceptType::FATAL, e.what());
      }
      catch (...) {
        Throw(ExceptType::UNCAUGHT, "non-standard exception");
      }
}
