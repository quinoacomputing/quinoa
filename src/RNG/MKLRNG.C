//******************************************************************************
/*!
  \file      src/RNG/MKLRNG.C
  \author    J. Bakosi
  \date      Tue Oct 22 15:43:10 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKL-based random number generator
  \details   MKL-based random number generator
*/
//******************************************************************************

#include <iostream>

#include <MKLRNG.h>
#include <Paradigm.h>
#include <Exception.h>

using namespace quinoa;

MKLRNG::MKLRNG(const Base& base) noexcept :
  m_nOMPthreads(base.paradigm.getOpenMP()->nthread()),
  m_table(),
  m_stream()
//******************************************************************************
//  Constructor
//! \param[in] base     Essentials
//! \details   No-throw guarantee: this member function never throws exceptions.
//! \author  J. Bakosi
//******************************************************************************
{
}

MKLRNG::~MKLRNG() noexcept
//******************************************************************************
//  Destructor: Free all random number tables and streams
//! \details No-throw guarantee: this member function never throws exceptions.
//! \author  J. Bakosi
//******************************************************************************
{
  for (tk::MKLRndTable* t : m_table) { delete t; }
  m_table.clear();
  for (tk::MKLRndStream* s : m_stream) { delete s; }
  m_stream.clear();
}

tk::MKLRndTable*
MKLRNG::addTable(const int brng,
                 const tk::RndDist dist,
                 const int method,
                 const unsigned int seed,
                 const long long int number,
                 const std::string name)
//******************************************************************************
//  Add a random number table
//! \param[in]  brng     Basic VSL generator type
//! \param[in]  dist     Type of distribution
//! \param[in]  method   Generation method
//! \param[in]  seed     Seed
//! \param[in]  number   Total number of random numbers in table
//! \param[in]  name     MemoryEntry name of the random number table
//! \author  J. Bakosi
//******************************************************************************
{
  // Create new table
  tk::MKLRndTable* table = new (std::nothrow)
    tk::MKLRndTable(m_nOMPthreads, brng, dist, method, seed, number, name);
  ErrChk(table != nullptr, tk::ExceptType::FATAL, "Cannot allocate memory");

  // Store new table
  std::pair<Tables::iterator,bool> e = m_table.insert(table);
  if (!e.second) {
    if (table) delete table;
    Throw(tk::ExceptType::FATAL, "Cannot store random number table");
  }

  // Return key to caller
  return table;
}

void
MKLRNG::eraseTable(tk::MKLRndTable* table) noexcept
//******************************************************************************
//  Erase a random number table
//! \param[in]  table    Pointer to table to erase
//! \details    Exception safety: no-throw guarantee: never throws exceptions.
//! \author  J. Bakosi
//******************************************************************************
{
  try {

    auto it = m_table.find(table);
    Assert(it != m_table.end(), tk::ExceptType::FATAL,
           "Cannot find random number table in MKLRNG::eraseTable()");

    delete table;
    m_table.erase(it);

  } // emit only a warning on error
    catch (tk::Exception& e) {
      e.echo("WARNING");
    }
    catch (std::exception& e) {
      std::cout << ">>> std::exception in MKLRNG::eraseTable(): "
                << e.what() << std::endl;
    }
    catch (...) {
      std::cout << ">>> UNKNOWN EXCEPTION in MKLRNG::eraseTable()"
                << std::endl;
    }
}

void
MKLRNG::regenTables()
//******************************************************************************
//  Regenerate random numbers in all tables
//! \author  J. Bakosi
//******************************************************************************
{
  for (tk::MKLRndTable* t : m_table) t->generate();
}

const tk::real*
MKLRNG::getRnd(tk::MKLRndTable* table)
//******************************************************************************
//  Constant accessor to random number table
//! \param[in]  table    Table id to access
//! \author  J. Bakosi
//******************************************************************************
{
  auto it = m_table.find(table);
  Assert(it != m_table.end(), tk::ExceptType::FATAL,
         "Cannot find random number table");
  return (*it)->getRnd();
}

tk::MKLRndStream*
MKLRNG::addStream(const int brng, const unsigned int seed)
//******************************************************************************
//  Add a random number stream
//! \param[in]  brng     Basic VSL generator type
//! \param[in]  seed     Seed
//! \author  J. Bakosi
//******************************************************************************
{
  // Create new stream
  tk::MKLRndStream* stream =
    new (std::nothrow) tk::MKLRndStream(m_nOMPthreads, brng, seed);
  ErrChk(stream != nullptr, tk::ExceptType::FATAL, "Cannot allocate memory");

  // Store new stream
  std::pair<Streams::iterator,bool> e = m_stream.insert(stream);
  if (!e.second) {
    if (stream) delete stream;
    Throw(tk::ExceptType::FATAL, "Cannot store random number stream");
  }

  // Return key to caller
  return stream;
}

void
MKLRNG::eraseStream(tk::MKLRndStream* stream) noexcept
//******************************************************************************
//  Erase a random number stream
//! \param[in]  stream    Pointer to stream to erase
//! \details    Exception safety: no-throw guarantee: never throws exceptions.
//! \author  J. Bakosi
//******************************************************************************
{
  try {

    auto it = m_stream.find(stream);
    Assert(it != m_stream.end(), tk::ExceptType::FATAL,
           "Cannot find random number stream in MKLRNG::eraseStream()");

    delete stream;
    m_stream.erase(it);

  } // emit only a warning on error
    catch (std::exception& e) {
      std::cout << "WARNING: " << e.what() << std::endl;
    }
    catch (...) {
      std::cout << "UNKNOWN EXCEPTION in MKLRNG::eraseStream()" << std::endl
                << "Continuing anyway..." << std::endl;
    }
}

const VSLStreamStatePtr*
MKLRNG::getStr(tk::MKLRndStream* stream)
//******************************************************************************
//  Constant accessor to VSL stream
//! \param[in]  stream    Stream id to access
//! \author  J. Bakosi
//******************************************************************************
{
  auto it = m_stream.find(stream);
  Assert(it != m_stream.end(), tk::ExceptType::FATAL,
         "Cannot find random number stream");
  return (*it)->getStr();
}
