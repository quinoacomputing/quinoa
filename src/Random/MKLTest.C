//******************************************************************************
/*!
  \file      src/Random/MKLTest.C
  \author    J. Bakosi
  \date      Mon Oct  7 14:33:31 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKL random number generator tests
  \details   MKL random number generator tests
*/
//******************************************************************************

#include <mkl_vsl.h>

#include <MKLTest.h>
#include <Exception.h>
#include <RNGTest/InputDeck/InputDeck.h>

using namespace rngtest;

MKLTest::MKLTest(const ctr::InputDeck& control) :
  m_rng(),
  m_testrng(control.get<ctr::generator>()),
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
  using tk::ExceptType;

  if (vslerr != VSL_STATUS_OK)
    try {

      std::stringstream s;
      s << "MKL VSL Error: code " << vslerr;
      Throw(ExceptType::FATAL, s.str());

    } catch (tk::Exception&) {
        throw;
      }
      catch (std::exception& e) {
        Throw(ExceptType::FATAL, e.what());
      }
      catch (...) {
        Throw(ExceptType::UNCAUGHT, "non-standard exception");
      }
}
