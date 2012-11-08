//******************************************************************************
/*!
  \file      src/Random/VSLException.C
  \author    J. Bakosi
  \date      Wed Nov  7 17:37:01 2012
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     VSLException class definition
  \details   VSLException class definition
*/
//******************************************************************************

#include <mkl_vsl.h>

#include <VSLException.h>

using namespace Quinoa;

VSLException::VSLException(ExceptType except, int vslerr) :
  MKLException(except, MKL_VSL_ERROR)
//******************************************************************************
//  Constructor: zero memory entry pointers held
//! \author J. Bakosi
//******************************************************************************
{
  // Fill VSLError -> VSLExceptType map
  // ICC: once initializer lists are supported this should be a constant map
  m_VSLErrMap.insert(
    make_pair<int,VSLExceptType>(VSL_ERROR_FEATURE_NOT_IMPLEMENTED,
                                 VSL_UNIMPLEMENTED));
  m_VSLErrMap.insert(
    make_pair<int,VSLExceptType>(VSL_ERROR_UNKNOWN,
                                 VSL_UNKNOWN));
  m_VSLErrMap.insert(
    make_pair<int,VSLExceptType>(VSL_ERROR_BADARGS,
                                 VSL_BADARGS));
  m_VSLErrMap.insert(
    make_pair<int,VSLExceptType>(VSL_ERROR_MEM_FAILURE,
                                 VSL_MEM_FAILURE));
  m_VSLErrMap.insert(
    make_pair<int,VSLExceptType>(VSL_ERROR_NULL_PTR,
                                 VSL_NULL_PTR));
  m_VSLErrMap.insert(
    make_pair<int,VSLExceptType>(VSL_RNG_ERROR_INVALID_BRNG_INDEX,
                                 VSL_INVALID_BRNG_INDEX));
  m_VSLErrMap.insert(
    make_pair<int,VSLExceptType>(VSL_RNG_ERROR_LEAPFROG_UNSUPPORTED,
                                 VSL_LEAPFROG_UNSUPPORTED));
  m_VSLErrMap.insert(
    make_pair<int,VSLExceptType>(VSL_RNG_ERROR_SKIPAHEAD_UNSUPPORTED,
                                 VSL_SKIPAHEAD_UNSUPPORTED));
  m_VSLErrMap.insert(
    make_pair<int,VSLExceptType>(VSL_RNG_ERROR_BRNGS_INCOMPATIBLE,
                                 VSL_BRNGS_INCOMPATIBLE));
  m_VSLErrMap.insert(
    make_pair<int,VSLExceptType>(VSL_RNG_ERROR_BAD_STREAM,
                                 VSL_BAD_STREAM));
  m_VSLErrMap.insert(
    make_pair<int,VSLExceptType>(VSL_RNG_ERROR_BRNG_TABLE_FULL,
                                 VSL_BRNG_TABLE_FULL));
  m_VSLErrMap.insert(
    make_pair<int,VSLExceptType>(VSL_RNG_ERROR_BAD_STREAM_STATE_SIZE,
                                 VSL_BAD_STREAM_STATE_SIZE));
  m_VSLErrMap.insert(
    make_pair<int,VSLExceptType>(VSL_RNG_ERROR_BAD_WORD_SIZE,
                                 VSL_BAD_WORD_SIZE));
  m_VSLErrMap.insert(
    make_pair<int,VSLExceptType>(VSL_RNG_ERROR_BAD_NSEEDS,
                                 VSL_BAD_NSEEDS));
  m_VSLErrMap.insert(
    make_pair<int,VSLExceptType>(VSL_RNG_ERROR_BAD_NBITS,
                                 VSL_BAD_NBITS));
  m_VSLErrMap.insert(
    make_pair<int,VSLExceptType>(VSL_RNG_ERROR_BAD_UPDATE,
                                 VSL_BAD_UPDATE));
  m_VSLErrMap.insert(
    make_pair<int,VSLExceptType>(VSL_RNG_ERROR_NO_NUMBERS,
                                 VSL_NO_NUMBERS));
  m_VSLErrMap.insert(
    make_pair<int,VSLExceptType>(VSL_RNG_ERROR_INVALID_ABSTRACT_STREAM,
                                 VSL_INVALID_ABSTRACT_STREAM));
  m_VSLErrMap.insert(
    make_pair<int,VSLExceptType>(VSL_RNG_ERROR_FILE_CLOSE,
                                 VSL_FILE_CLOSE));
  m_VSLErrMap.insert(
    make_pair<int,VSLExceptType>(VSL_RNG_ERROR_FILE_OPEN,
                                 VSL_FILE_OPEN));
  m_VSLErrMap.insert(
    make_pair<int,VSLExceptType>(VSL_RNG_ERROR_FILE_WRITE,
                                 VSL_FILE_WRITE));
  m_VSLErrMap.insert(
    make_pair<int,VSLExceptType>(VSL_RNG_ERROR_FILE_READ,
                                 VSL_FILE_READ));
  m_VSLErrMap.insert(
    make_pair<int,VSLExceptType>(VSL_RNG_ERROR_BAD_FILE_FORMAT,
                                 VSL_BAD_FILE_FORMAT));
  m_VSLErrMap.insert(
    make_pair<int,VSLExceptType>(VSL_RNG_ERROR_UNSUPPORTED_FILE_VER,
                                 VSL_UNSUPPORTED_FILE_VER));

  // Set VSL exception type based on VSL error
  m_except = getException(vslerr);
}

VSLExceptType
VSLException::getException(int vslerr)
//******************************************************************************
//  Get VSLException based on VSLError
//! \author J. Bakosi
//******************************************************************************
{
  auto it = m_VSLErrMap.find(vslerr);
  if (it == m_VSLErrMap.end())
    return VSL_UNIMPLEMENTED;
  else
    return it->second;
}

ErrCode
VSLException::handleException(Driver* driver)
//******************************************************************************
//  Handle VSLException
//! \author J. Bakosi
//******************************************************************************
{
  // Start error message
  m_message = VSLMsg[static_cast<int>(m_except)];

  // Handle Exception (criticality)
  return MKLException::handleException(driver);
}
