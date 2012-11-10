//******************************************************************************
/*!
  \file      src/Random/MKLRndTable.h
  \author    J. Bakosi
  \date      Fri 09 Nov 2012 06:48:40 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generation into tables using Intel's MKL
  \details   Tables are used to generate a fixed large number of fixed property
             random numbers by several threads using block-splitting.
*/
//******************************************************************************
#ifndef MKLRndTable_h
#define MKLRndTable_h

#include <QuinoaTypes.h>
#include <Memory.h>
#include <MKL.h>

namespace Quinoa {

//! Probability distributions for sampling into tables
enum RndDist { UNIFORM=0,        //!< Uniform
               GAUSSIAN,         //!< Gaussian
               GAMMA,            //!< Gamma
               NUM_DIST_TYPES
};

//! Constants for sampling the uniform distribution in tables
const real UNIFORM_LEFT_BOUND = 0.0;
const real UNIFORM_RIGHT_BOUND = 1.0;

//! Constants for sampling the Gaussian distribution in tables
const real GAUSSIAN_MEAN = 0.0;
const real GAUSSIAN_STD = 1.0;

//! Constants for sampling the gamma distribution in tables
const real GAMMA_MEAN = 1.0;
const real GAMMA_VAR = 0.25;
const real GAMMA_SHAPE = GAMMA_MEAN * GAMMA_MEAN / GAMMA_VAR;
const real GAMMA_DISPLACEMENT = 0.0;
const real GAMMA_SCALE = GAMMA_VAR / GAMMA_MEAN;

//! MKL-based random number generator into tables using block-splitting
class MKLRndTable : private MKL {

  public:
    //! Constructor: Create random number skip-ahead table
    MKLRndTable(Memory* memory,
                const int nthread,
                const int brng,
                const RndDist dist,
                const int method,
                const unsigned int seed,
                const long long int number,
                const string name);

    //! Destructor: Destroy random number skip-ahead table
    ~MKLRndTable();

    //! Regenerate random numbers in table
    void generate();

    //! Constant accessor to random number table
    const real* getRnd() { return m_rndPtr; }

  private:
    //! Don't permit copy constructor
    MKLRndTable(const MKLRndTable&) = delete;
    //! Don't permit copy assigment
    MKLRndTable& operator=(const MKLRndTable&) = delete;
    //! Don't permit move constructor
    MKLRndTable(MKLRndTable&&) = delete;
    //! Don't permit move assigment
    MKLRndTable& operator=(MKLRndTable&&) = delete;

    Memory* m_memory;                //!< Memory object pointer
    const int m_nthread;             //!< Number of threads to use
    const RndDist m_dist;            //!< RndDist (UNIFORM, GAUSSIAN, etc.)
    const int m_method;              //!< Generation method (dist. specific)
    const long long int m_chunk;     //!< Number of numbers generated per thread
    const long long int m_remainder; //!< Leftover
    VSLStreamStatePtr* m_stream;     //!< Array of pointers to thread-streams
    MemoryEntry* m_rnd;              //!< MemoryEntry Pointer to random numbers
    real* m_rndPtr;                  //!< Raw pointer to random numbers
};

} // namespace Quinoa

#endif // MKLRndTable_h
