//******************************************************************************
/*!
  \file      src/Random/MKLRndTable.h
  \author    J. Bakosi
  \date      Wed Sep 25 15:50:52 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generation into tables using Intel's MKL
  \details   Tables are used to generate a fixed large number of fixed property
             random numbers by several threads using block-splitting.
*/
//******************************************************************************
#ifndef MKLRndTable_h
#define MKLRndTable_h

#include <memory>

#include <QuinoaTypes.h>
#include <Memory.h>
#include <MKL.h>

namespace quinoa {

//! Probability distributions for sampling into tables
enum class RndDist : uint8_t { UNIFORM=0,        //!< Uniform
                               GAUSSIAN,         //!< Gaussian
                               GAMMA             //!< Gamma
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
class MKLRndTable : public MKL {

  public:
    //! Constructor: Create random number skip-ahead table
    explicit MKLRndTable(int nthread,
                         int brng,
                         RndDist dist,
                         int method,
                         unsigned int seed,
                         long long int number,
                         const std::string& name);

    //! Destructor: Destroy random number skip-ahead table
    ~MKLRndTable() noexcept override;

    //! Regenerate random numbers in table
    void generate() const;

    //! Constant accessor to random number table
    const real* getRnd() const noexcept { return m_rnd.get(); }

  private:
    //! Don't permit copy constructor
    MKLRndTable(const MKLRndTable&) = delete;
    //! Don't permit copy assigment
    MKLRndTable& operator=(const MKLRndTable&) = delete;
    //! Don't permit move constructor
    MKLRndTable(MKLRndTable&&) = delete;
    //! Don't permit move assigment
    MKLRndTable& operator=(MKLRndTable&&) = delete;

    //! Finalize, single exit point, called implicitly from destructor or
    //! explicitly from anywhere else
    void finalize() noexcept;

    const int m_nthread;             //!< Number of threads to use
    const RndDist m_dist;            //!< RndDist (UNIFORM, GAUSSIAN, etc.)
    const int m_method;              //!< Generation method (dist. specific)
    const long long int m_chunk;     //!< Number of numbers generated per thread

    const long long int m_remainder; //!< Leftover
    std::unique_ptr<VSLStreamStatePtr[]> m_stream; //!< Thread-streams
    std::unique_ptr<real[]> m_rnd;   //!< Random numbers
};

} // namespace quinoa

#endif // MKLRndTable_h
