//******************************************************************************
/*!
  \file      src/Random/MKLRandom.h
  \author    J. Bakosi
  \date      Mon Oct  7 10:27:52 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKL-based random number generator
  \details   MKL-based random number generator
*/
//******************************************************************************
#ifndef MKLRandom_h
#define MKLRandom_h

#include <unordered_set>

#include <Base.h>
#include <Random.h>
#include <MKLRndTable.h>
#include <MKLRndStream.h>

namespace quinoa {

//! MKL-based random number generator
class MKLRandom : public tk::Random {

  public:
    //! Constructor
    MKLRandom(const Base& base) noexcept;

    //! Destructor: Free all random number tables and streams
    virtual ~MKLRandom() noexcept;

    //! Add random table
    tk::MKLRndTable* addTable(const int brng,
                              const tk::RndDist dist,
                              const int method,
                              const unsigned int seed,
                              const long long int number,
                              const std::string name);

    //! Erase a random number table
    void eraseTable(tk::MKLRndTable* table) noexcept;

    //! Regenerate random numbers in all tables
    void regenTables();

    //! Constant accessor to random number table
    const tk::real* getRnd(tk::MKLRndTable* table);

    //! Add a random number stream
    tk::MKLRndStream* addStream(const int brng, const unsigned int seed);

    //! Erase a random number stream
    void eraseStream(tk::MKLRndStream* stream) noexcept;

    //! Constant accessor to random number VSL stream
    const VSLStreamStatePtr* getStr(tk::MKLRndStream* stream);

  private:
    //! Don't permit copy constructor
    MKLRandom(const MKLRandom&) = delete;
    //! Don't permit copy assigment
    MKLRandom& operator=(const MKLRandom&) = delete;
    //! Don't permit move constructor
    MKLRandom(MKLRandom&&) = delete;
    //! Don't permit move assigment
    MKLRandom& operator=(MKLRandom&&) = delete;

    const int m_nOMPthreads;           //!< Number of OpenMP threads

    //! Type for a set of stream-tables to generate a large (and fixed) number
    //! of random numbers with fixed properties using several threads
    using Tables = std::unordered_set<tk::MKLRndTable*>;

    //! Type for a set of streams to generate a few  random numbers with
    //! arbitrary properties using several threads with leap-frogging
    using Streams = std::unordered_set<tk::MKLRndStream*>;

    //! Stream tables to generate fixed numbers of random numbers with fixed
    //! properties using several threads
    Tables m_table;

    //! Streams to generate a few random numbers at a time leap-frogging
    Streams m_stream;
};

} // quinoa::

#endif // MKLRandom_h
