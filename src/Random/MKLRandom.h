//******************************************************************************
/*!
  \file      src/Random/MKLRandom.h
  \author    J. Bakosi
  \date      Fri 19 Oct 2012 07:40:27 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKL-based random number generator
  \details   MKL-based random number generator
*/
//******************************************************************************
#ifndef MKLRandom_h
#define MKLRandom_h

#include <unordered_set>

#include <Random.h>
#include <MKLRndTable.h>
#include <MKLRndStream.h>

using namespace std;

namespace Quinoa {

//! MKL-based random number generator
class MKLRandom : Random {

  public:
    //! Constructor
    MKLRandom(const long long int nthread, Memory* memory) :
      m_nthread(nthread), m_memory(memory) {}

    //! Destructor: Free all random number tables and streams
    ~MKLRandom();

    //! Add random table
    MKLRndTable* addTable(const int brng,
                          const Distribution dist,
                          const int method,
                          const unsigned int seed,
                          const long long int number,
                          const string name);

    //! Erase a random number table
    void eraseTable(MKLRndTable* table);

    //! Regenerate random numbers in all tables
    void regenTables();

    //! Constant accessor to random number table
    const real* getRnd(MKLRndTable* table);

  private:
    //! Don't permit copy constructor
    MKLRandom(const MKLRandom&) = delete;
    //! Don't permit copy assigment
    MKLRandom& operator=(const MKLRandom&) = delete;
    //! Don't permit move constructor
    MKLRandom(MKLRandom&&) = delete;
    //! Don't permit move assigment
    MKLRandom& operator=(MKLRandom&&) = delete;

    const int m_nthread;    //!< Number of threads to use
    Memory* m_memory;       //!< Memory object pointer

    //! Array of pointers to thread-streams (for sampling a few at a time)
    //VSLStreamStatePtr* m_stream;

    //! Type for a set of stream-tables to generate a large (and fixed) number
    //! of random numbers with fixed properties using several threads
    typedef unordered_set<MKLRndTable*> Tables;

    //! Stream tables to generate fixed numbers of random numbers with fixed
    //! properties using several threads
    Tables m_table;
};

} // namespace Quinoa

#endif // MKLRandom_h
