//******************************************************************************
/*!
  \file      src/Random/MKLRandom.h
  \author    J. Bakosi
  \date      Tue 13 Nov 2012 10:12:03 PM MST
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
class MKLRandom : private Random {

  public:
    //! Constructor
    MKLRandom(Memory* memory, Paradigm* paradigm);

    //! Destructor: Free all random number tables and streams
    ~MKLRandom();

    //! Add random table
    MKLRndTable* addTable(const int brng,
                          const RndDist dist,
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

    //! Add a random number stream
    MKLRndStream* addStream(const int brng, const unsigned int seed);

    //! Erase a random number stream
    void eraseStream(MKLRndStream* stream);

    //! Constant accessor to random number VSL stream
    const VSLStreamStatePtr* getStr(MKLRndStream* stream);

  private:
    //! Don't permit copy constructor
    MKLRandom(const MKLRandom&) = delete;
    //! Don't permit copy assigment
    MKLRandom& operator=(const MKLRandom&) = delete;
    //! Don't permit move constructor
    MKLRandom(MKLRandom&&) = delete;
    //! Don't permit move assigment
    MKLRandom& operator=(MKLRandom&&) = delete;

    Memory* m_memory;            //!< Memory object pointer
    Paradigm* m_paradigm;        //!< Pointer to Memory object
    int m_nOMPthreads;           //!< Number of OpenMP threads

    //! Type for a set of stream-tables to generate a large (and fixed) number
    //! of random numbers with fixed properties using several threads
    using Tables = unordered_set<MKLRndTable*>;

    //! Type for a set of streams to generate a few  random numbers with
    //! arbitrary properties using several threads with leap-frogging
    using Streams = unordered_set<MKLRndStream*>;

    //! Stream tables to generate fixed numbers of random numbers with fixed
    //! properties using several threads
    Tables m_table;

    //! Streams to generate a few random numbers at a time leap-frogging
    Streams m_stream;
};

} // namespace Quinoa

#endif // MKLRandom_h
