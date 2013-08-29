//******************************************************************************
/*!
  \file      src/Random/MKLRandom.h
  \author    J. Bakosi
  \date      Thu Aug 29 15:26:23 2013
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

namespace quinoa {

//! MKL-based random number generator
class MKLRandom : public Random {

  public:
    //! Constructor
    MKLRandom(Memory* const memory, Paradigm* const paradigm) noexcept;

    //! Destructor: Free all random number tables and streams
    virtual ~MKLRandom() noexcept;

    //! Add random table
    MKLRndTable* addTable(const int brng,
                          const RndDist dist,
                          const int method,
                          const unsigned int seed,
                          const long long int number,
                          const std::string name);

    //! Erase a random number table
    void eraseTable(MKLRndTable* table) noexcept;

    //! Regenerate random numbers in all tables
    void regenTables();

    //! Constant accessor to random number table
    const real* getRnd(MKLRndTable* table);

    //! Add a random number stream
    MKLRndStream* addStream(const int brng, const unsigned int seed);

    //! Erase a random number stream
    void eraseStream(MKLRndStream* stream) noexcept;

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

    Memory* const m_memory;            //!< Memory object pointer
    const int m_nOMPthreads;           //!< Number of OpenMP threads

    //! Type for a set of stream-tables to generate a large (and fixed) number
    //! of random numbers with fixed properties using several threads
    using Tables = std::unordered_set<MKLRndTable*>;

    //! Type for a set of streams to generate a few  random numbers with
    //! arbitrary properties using several threads with leap-frogging
    using Streams = std::unordered_set<MKLRndStream*>;

    //! Stream tables to generate fixed numbers of random numbers with fixed
    //! properties using several threads
    Tables m_table;

    //! Streams to generate a few random numbers at a time leap-frogging
    Streams m_stream;
};

} // namespace quinoa

#endif // MKLRandom_h
