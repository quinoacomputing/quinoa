//******************************************************************************
/*!
  \file      src/Random/MKLTest.C
  \author    J. Bakosi
  \date      Sat 24 Aug 2013 07:29:51 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKL random number generator tests
  \details   MKL random number generator tests
*/
//******************************************************************************

#include <mkl_vsl.h>

#include <MKLTest.h>
#include <Exception.h>
#include <Control.h>

using namespace Quinoa;

MKLTest::MKLTest(Control* const control) :
  m_rng(),
  m_testrng(control->get<control::RNGS>()),
  m_stream()
//******************************************************************************
//  Constructor
//! \author  J. Bakosi
//******************************************************************************
{
//   for (auto& r : m_testrng) {
//     if (lib(r) == RNGLibType::MKL) {
//   #ifdef MKL_CALLS
//   #ifdef NDEBUG
//     vslNewStream(&stream, brng, seed);
//   #else  // NDEBUG
//     errchk(vslNewStream(&stream, brng, seed));
//   #endif // NDEBUG
//   #endif
//     }
//   }
}

MKLTest::~MKLTest()
//******************************************************************************
//  Destructor
//! \author  J. Bakosi
//******************************************************************************
{

}

void
MKLTest::MKLErrChk(int vslerr) const
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
