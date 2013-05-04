//******************************************************************************
/*!
  \file      src/Random/MKLRandom.C
  \author    J. Bakosi
  \date      Fri 03 May 2013 07:05:49 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKL-based random number generator
  \details   MKL-based random number generator
*/
//******************************************************************************

#include <cassert>
#include <iostream>

#include <MKLRandom.h>
#include <Paradigm.h>

using namespace Quinoa;

MKLRandom::MKLRandom(Memory* const memory, Paradigm* const paradigm) :
  m_memory(memory),
  m_nOMPthreads(paradigm->getOpenMP()->nthread())
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \author  J. Bakosi
//******************************************************************************
{
}

MKLRandom::~MKLRandom() noexcept
//******************************************************************************
//  Destructor: Free all random number tables and streams
//! \details No-throw guarantee: this member function never throws exceptions.
//! \author  J. Bakosi
//******************************************************************************
{
  for (MKLRndTable* t : m_table) { delete t; }
  m_table.clear();
  for (MKLRndStream* s : m_stream) { delete s; }
  m_stream.clear();
}

MKLRndTable*
MKLRandom::addTable(const int brng,
                    const RndDist dist,
                    const int method,
                    const unsigned int seed,
                    const long long int number,
                    const string name)
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
  MKLRndTable* table = new (nothrow)
    MKLRndTable(m_memory, m_nOMPthreads, brng, dist, method, seed, number, name);
  if (table == nullptr) throw Exception(FATAL, "Cannot allocate memory");

  // Store new table
  pair<Tables::iterator,bool> e = m_table.insert(table);
  if (!e.second) {
    if (table) delete table;
    throw Exception(FATAL, "Cannot store random number table");
  }

  // Return key to caller
  return table;
}

void
MKLRandom::eraseTable(MKLRndTable* table) noexcept
//******************************************************************************
//  Erase a random number table
//! \param[in]  table    Pointer to table to erase
//! \details    Exception safety: no-throw guarantee: never throws exceptions.
//! \author  J. Bakosi
//******************************************************************************
{
  try {

    auto it = m_table.find(table);
    assert(it != m_table.end());

    delete table;
    m_table.erase(it);

  } // emit only a warning on error
    catch (exception& e) {
      cout << "WARNING: " << e.what() << endl;
    }
    catch (...) {
      cout << "UNKNOWN EXCEPTION in MKLRandom::eraseTable()" << endl
           << "Continuing anyway..." << endl;
    }
}

void
MKLRandom::regenTables()
//******************************************************************************
//  Regenerate random numbers in all tables
//! \author  J. Bakosi
//******************************************************************************
{
  for (MKLRndTable* t : m_table) t->generate();
}

const real*
MKLRandom::getRnd(MKLRndTable* table)
//******************************************************************************
//  Constant accessor to random number table
//! \param[in]  table    Table id to access
//! \author  J. Bakosi
//******************************************************************************
{
  auto it = m_table.find(table);
  assert(it != m_table.end());
  return (*it)->getRnd();
}

MKLRndStream*
MKLRandom::addStream(const int brng, const unsigned int seed)
//******************************************************************************
//  Add a random number stream
//! \param[in]  brng     Basic VSL generator type
//! \param[in]  seed     Seed
//! \author  J. Bakosi
//******************************************************************************
{
  // Create new stream
  MKLRndStream* stream = new (nothrow) MKLRndStream(m_nOMPthreads, brng, seed);
  if (stream == nullptr) throw Exception(FATAL, "Cannot allocate memory");

  // Store new stream
  pair<Streams::iterator,bool> e = m_stream.insert(stream);
  if (!e.second) {
    if (stream) delete stream;
    throw Exception(FATAL, "Cannot store random number stream");
  }

  // Return key to caller
  return stream;
}

void
MKLRandom::eraseStream(MKLRndStream* stream) noexcept
//******************************************************************************
//  Erase a random number stream
//! \param[in]  stream    Pointer to stream to erase
//! \details    Exception safety: no-throw guarantee: never throws exceptions.
//! \author  J. Bakosi
//******************************************************************************
{
  try {

    auto it = m_stream.find(stream);
    assert(it != m_stream.end());

    delete stream;
    m_stream.erase(it);

  } // emit only a warning on error
    catch (exception& e) {
      cout << "WARNING: " << e.what() << endl;
    }
    catch (...) {
      cout << "UNKNOWN EXCEPTION in MKLRandom::eraseStream()" << endl
           << "Continuing anyway..." << endl;
    }
}

const VSLStreamStatePtr*
MKLRandom::getStr(MKLRndStream* stream)
//******************************************************************************
//  Constant accessor to VSL stream
//! \param[in]  stream    Stream id to access
//! \author  J. Bakosi
//******************************************************************************
{
  auto it = m_stream.find(stream);
  assert(it != m_stream.end());
  return (*it)->getStr();
}
