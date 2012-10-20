// //******************************************************************************
// /*!
//   \file      src/Random/MKLRandom.h
//   \author    J. Bakosi
//   \date      Fri 19 Oct 2012 02:25:32 PM MDT
//   \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
//   \brief     MKL-based random number generator
//   \details   MKL-based random number generator
// */
// //******************************************************************************
// #ifndef MKLRandom_h
// #define MKLRandom_h
// 
// #include <unordered_set>
// 
// #include <mkl_vsl.h>
// 
// #include <Memory.h>
// #include <Random.h>
// 
// using namespace std;
// 
// namespace Quinoa {
// 
// //! Probability distributions for sampling into tables
// enum Distribution { UNIFORM=0,        //!< Uniform
//                     GAUSSIAN,         //!< Gaussian
//                     GAMMA,            //!< Gamma
//                     NUM_DIST_TYPES
// };
// 
// //! Constants for sampling the uniform distribution in tables
// const Real UNIFORM_LEFT_BOUND = 0.0;
// const Real UNIFORM_RIGHT_BOUND = 1.0;
// 
// //! Constants for sampling the Gaussian distribution in tables
// const Real GAUSSIAN_MEAN = 0.0;
// const Real GAUSSIAN_STD = 1.0;
// 
// //! Constants for sampling the gamma distribution in tables
// const Real GAMMA_MEAN = 1.0;
// const Real GAMMA_VAR = 0.25;
// const Real GAMMA_SHAPE = GAMMA_MEAN * GAMMA_MEAN / GAMMA_VAR;
// const Real GAMMA_DISPLACEMENT = 0.0;
// const Real GAMMA_SCALE = GAMMA_VAR / GAMMA_MEAN;
// 
// //! MKL-based random number generator
// class MKLRandom : Random {
// 
//   //! Type for a set of stream-tables to generate a large (and fixed) number of
//   //! random numbers with fixed properties using Random::m_nthreads
//   typedef unordered_set<RndTable*> Tables;
// 
//   public:
//     //! Constructor: Setup random number generator streams
//     MKLRandom(const int brng,
//               const long long int nthreads,
//               const uInt seed,
//               Memory* memory);
// 
//     //! Destructor: Destroy random number generator thread-streams and tables
//     ~MKLRandom();
// 
//     //! Add random table
//     void addTable(const int brng,
//                   const Distribution dist,
//                   const int method,
//                   const long long int number,
//                   const string name);
// 
//     //! Regenerate random numbers in all tables
//     void regenAllTables();
// 
//     //! Thread-streams accessor (for sampling a few numbers at a time)
//     const VSLStreamStatePtr* getStream() const { return m_stream; }
// 
//     //! Call MKL's vdRngUniform() and handle error
//     void uniform(const Int& method,
//                  VSLStreamStatePtr& stream,
//                  const Int& n,
//                  Real* r,
//                  const Real& a,
//                  const Real& b);
// 
//     //! Call MKL's vdRngGaussian() and handle error
//     void gaussian(const Int& method,
//                   VSLStreamStatePtr& stream,
//                   const Int& n,
//                   Real* r,
//                   const Real& a,
//                   const Real& b);
// 
//     //! Call MKL's vdRngGamma() and handle error
//     void gamma(const Int& method,
//                VSLStreamStatePtr& stream,
//                const Int& n,
//                Real* r,
//                const Real& alpha,
//                const Real& a,
//                const Real& beta);
// 
//   private:
//     //! Don't permit copy constructor
//     MKLRandom(const MKLRandom&) = delete;
//     //! Don't permit copy assigment
//     MKLRandom& operator=(const MKLRandom&) = delete;
//     //! Don't permit move constructor
//     MKLRandom(MKLRandom&&) = delete;
//     //! Don't permit move assigment
//     MKLRandom& operator=(MKLRandom&&) = delete;
// 
//     //! Call MKL's vslNewStream() and handle error
//     void newStream(VSLStreamStatePtr* stream,
//                    const int& brng,
//                    const unsigned int& seed);
// 
//     //! Call MKL's vslCopyStream() and handle error
//     void copyStream(VSLStreamStatePtr* newstream,
//                     const VSLStreamStatePtr& srcstream);
// 
//     //! Call MKL's vslSkipaheadStream() and handle error
//     void skipAheadStream(VSLStreamStatePtr& stream,
//                          const long long int& nskip);
// 
//     //! Call MKL's vslLeapfrogStream() and handle error
//     void leapfrogStream(VSLStreamStatePtr& stream,
//                         const int& k,
//                         const int& nstreams);
// 
//     //! Regenerate random numbers in a table
//     void regenTable(const RndTable* table);
// 
//     //! Memory object pointer
//     Memory* m_memory;
// 
//     //! Array of pointers to thread-streams (for sampling a few at a time)
//     VSLStreamStatePtr* m_stream;
// 
//     //! Stream tables to generate fixed numbers of random numbers with fixed
//     //! properties using Random::m_nthreads
//     Tables m_table;
// };
// 
// } // namespace Quinoa
// 
// #endif // MKLRandom_h
