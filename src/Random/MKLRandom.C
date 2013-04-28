//******************************************************************************
/*!
  \file      src/Random/MKLRandom.C
  \author    J. Bakosi
  \date      Sat 27 Apr 2013 08:36:59 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKL-based random number generator
  \details   MKL-based random number generator
*/
//******************************************************************************

#include <MKLRandom.h>
#include <MemoryException.h>
#include <MKLException.h>
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
  Assert(table != nullptr, MemoryException,FATAL,BAD_ALLOC);

  // Store new table
  pair<Tables::iterator,bool> e = m_table.insert(table);
  if (!e.second) {
    if (table) delete table;
    Throw(MemoryException,FATAL,BAD_INSERT);
  }

  // Return key to caller
  return table;
}

void
MKLRandom::eraseTable(MKLRndTable* table)
//******************************************************************************
//  Erase a random number table
//! \param[in]  table    Pointer to table to erase
//! \author  J. Bakosi
//******************************************************************************
{
  auto it = m_table.find(table);
  if (it != m_table.end()) {
    delete table;
    m_table.erase(it);
  } else {
    Throw(MKLException,WARNING,MKL_UNKNOWN_TABLE);
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
  if (it != m_table.end()) {
    return (*it)->getRnd();
  } else {
    Throw(MKLException,WARNING,MKL_UNKNOWN_TABLE);
    return nullptr;
  }
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
  Assert(stream != nullptr, MemoryException,FATAL,BAD_ALLOC);

  // Store new stream
  pair<Streams::iterator,bool> e = m_stream.insert(stream);
  if (!e.second) {
    if (stream) delete stream;
    Throw(MemoryException,FATAL,BAD_INSERT);
  }

  // Return key to caller
  return stream;
}

void
MKLRandom::eraseStream(MKLRndStream* stream)
//******************************************************************************
//  Erase a random number stream
//! \param[in]  stream    Pointer to stream to erase
//! \author  J. Bakosi
//******************************************************************************
{
  auto it = m_stream.find(stream);
  if (it != m_stream.end()) {
    delete stream;
    m_stream.erase(it);
  } else {
    Throw(MKLException,WARNING,MKL_UNKNOWN_STREAM);
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
  if (it != m_stream.end()) {
    return (*it)->getStr();
  } else {
    Throw(MKLException,WARNING,MKL_UNKNOWN_STREAM);
    return nullptr;
  }
}
