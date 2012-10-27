//******************************************************************************
/*!
  \file      src/Random/MKLRandom.C
  \author    J. Bakosi
  \date      Sat 27 Oct 2012 11:14:49 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKL-based random number generator
  \details   MKL-based random number generator
*/
//******************************************************************************

#include <MKLRandom.h>
#include <MemoryException.h>
#include <MKLException.h>

using namespace Quinoa;

MKLRandom::~MKLRandom()
//******************************************************************************
//  Destructor: Free all random number tables
//! \author  J. Bakosi
//******************************************************************************
{
  for (MKLRndTable* t : m_table) { delete t; }
  m_table.clear();
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
    MKLRndTable(m_memory, m_nthread, brng, dist, method, seed, number, name);
  if (table == nullptr) throw MemoryException(FATAL, BAD_ALLOC);

  // Store new table
  pair<Tables::iterator,bool> e = m_table.insert(table);
  if (!e.second) {
    if (table) delete table;
    throw MemoryException(FATAL, BAD_INSERT);
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
    throw MKLException(WARNING, MKLEXCEPT_UNKNOWN_TABLE);
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
    throw MKLException(WARNING, MKLEXCEPT_UNKNOWN_TABLE);
  }
}
