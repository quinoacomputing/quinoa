//******************************************************************************
/*!
  \file      src/Random/MKLException.C
  \author    J. Bakosi
  \date      Mon 15 Oct 2012 09:23:09 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKLException class definition
  \details   MKLException class definition
*/
//******************************************************************************

#include <mkl_vsl.h>

#include <MKLException.h>

using namespace Quinoa;

MKLException::MKLException(ExceptType except, Int vslerr) :
  RandomException(except, RND_MKL)
//******************************************************************************
//  Constructor: zero memory entry pointers held
//! \author J. Bakosi
//******************************************************************************
{
  // Fill VSLError -> MKLExceptType map
  // ICC: once initializer lists are supported this should be a constant map
  m_VSLErrMap.insert(
    make_pair<Int,MKLExceptType>(VSL_ERROR_FEATURE_NOT_IMPLEMENTED,
                                 MKL_UNIMPLEMENTED));
  m_VSLErrMap.insert(
    make_pair<Int,MKLExceptType>(VSL_ERROR_UNKNOWN,
                                 MKL_UNKNOWN));
  m_VSLErrMap.insert(
    make_pair<Int,MKLExceptType>(VSL_ERROR_BADARGS,
                                 MKL_BADARGS));
  m_VSLErrMap.insert(
    make_pair<Int,MKLExceptType>(VSL_ERROR_MEM_FAILURE,
                                 MKL_MEM_FAILURE));
  m_VSLErrMap.insert(
    make_pair<Int,MKLExceptType>(VSL_ERROR_NULL_PTR,
                                 MKL_NULL_PTR));
  m_VSLErrMap.insert(
    make_pair<Int,MKLExceptType>(VSL_ERROR_INVALID_BRNG_INDEX,
                                 MKL_INVALID_BRNG_INDEX));
  m_VSLErrMap.insert(
    make_pair<Int,MKLExceptType>(VSL_ERROR_LEAPFROG_UNSUPPORTED,
                                 MKL_LEAPFROG_UNSUPPORTED));
  m_VSLErrMap.insert(
    make_pair<Int,MKLExceptType>(VSL_ERROR_SKIPAHEAD_UNSUPPORTED,
                                 MKL_SKIPAHEAD_UNSUPPORTED));
  m_VSLErrMap.insert(
    make_pair<Int,MKLExceptType>(VSL_ERROR_BRNGS_INCOMPATIBLE,
                                 MKL_BRNGS_INCOMPATIBLE));
  m_VSLErrMap.insert(
    make_pair<Int,MKLExceptType>(VSL_ERROR_BAD_STREAM,
                                 MKL_BAD_STREAM));
  m_VSLErrMap.insert(
    make_pair<Int,MKLExceptType>(VSL_ERROR_BRNG_TABLE_FULL,
                                 MKL_BRNG_TABLE_FULL));
  m_VSLErrMap.insert(
    make_pair<Int,MKLExceptType>(VSL_ERROR_BAD_STREAM_STATE_SIZE,
                                 MKL_BAD_STREAM_STATE_SIZE));
  m_VSLErrMap.insert(
    make_pair<Int,MKLExceptType>(VSL_ERROR_BAD_WORD_SIZE,
                                 MKL_BAD_WORD_SIZE));
  m_VSLErrMap.insert(
    make_pair<Int,MKLExceptType>(VSL_ERROR_BAD_NSEEDS,
                                 MKL_BAD_NSEEDS));
  m_VSLErrMap.insert(
    make_pair<Int,MKLExceptType>(VSL_ERROR_BAD_NBITS,
                                 MKL_BAD_NBITS));
  m_VSLErrMap.insert(
    make_pair<Int,MKLExceptType>(VSL_ERROR_BAD_UPDATE,
                                 MKL_BAD_UPDATE));
  m_VSLErrMap.insert(
    make_pair<Int,MKLExceptType>(VSL_ERROR_NO_NUMBERS,
                                 MKL_NO_NUMBERS));
  m_VSLErrMap.insert(
    make_pair<Int,MKLExceptType>(VSL_ERROR_INVALID_ABSTRACT_STREAM,
                                 MKL_INVALID_ABSTRACT_STREAM));
  m_VSLErrMap.insert(
    make_pair<Int,MKLExceptType>(VSL_ERROR_FILE_CLOSE,
                                 MKL_ERROR_FILE_CLOSE));
  m_VSLErrMap.insert(
    make_pair<Int,MKLExceptType>(VSL_ERROR_FILE_OPEN,
                                 MKL_ERROR_FILE_OPEN));
  m_VSLErrMap.insert(
    make_pair<Int,MKLExceptType>(VSL_ERROR_FILE_WRITE,
                                 MKL_ERROR_FILE_WRITE));
  m_VSLErrMap.insert(
    make_pair<Int,MKLExceptType>(VSL_ERROR_FILE_READ,
                                 MKL_ERROR_FILE_READ));
  m_VSLErrMap.insert(
    make_pair<Int,MKLExceptType>(VSL_ERROR_BAD_FILE_FORMAT,
                                 MKL_BAD_FILE_FORMAT));
  m_VSLErrMap.insert(
    make_pair<Int,MKLExceptType>(VSL_ERROR_UNSUPPORTED_FILE_VER,
                                 MKL_UNSUPPORTED_FILE_VER));

  // Set MKL exception type based on VSL error
  m_except = getException(vslerr);
}

MKLExceptType
MKLException::getException(Int vslerr)
//******************************************************************************
//  Get MKLException based on VSLError
//! \author J. Bakosi
//******************************************************************************
{
  auto it = m_VSLErrMap.find(vslerr);
  if (it == m_VSLErrMap.end())
    return MKLEXCEPT_UNIMPLEMENTED;
  else
    return it->second;
}

ErrCode
MKLException::handleException(Driver* driver)
//******************************************************************************
//  Handle MKLException
//! \author J. Bakosi
//******************************************************************************
{
  // Start error message
  message = MKLMsg[static_cast<Int>(m_except)];

  // Handle Exception (criticality)
  return RandomException::handleException(driver);
}
