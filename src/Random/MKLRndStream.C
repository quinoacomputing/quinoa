//******************************************************************************
/*!
  \file      src/Random/MKLRndStream.C
  \author    J. Bakosi
  \date      Fri 19 Oct 2012 03:14:34 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKL-based random number generator using leap frogging
  \details   MKL-based random number generator using leap forgging
*/
//******************************************************************************
// 
// MKLRandom::MKLRandom(const int brng,
//                      const long long int nthreads,
//                      const uInt seed,
//                      Memory* memory) :
//   Random(nthreads, seed), m_memory(memory)
// //******************************************************************************
// //  Constructor
// //! \param[in]   brng     Index of the basic generator to initialize the stream
// //! \param[in]   nthreads Initialize generators using nthreads threads
// //! \param[in]   seed     Initial condition of the stream
// //! \param[in]   memory   Memory store object pointer
// //! \details Initialize random number generator thread-streams for sampling a
// //!          few numbers at a time.
// //! \author  J. Bakosi
// //******************************************************************************
// {
//   // Allocate array of thread-stream pointers (used for a few at a time),
//   // initialize to zero
//   try {
//     m_stream = new VSLStreamStatePtr [m_nthreads]();
//   } catch (bad_alloc&) { throw MemoryException(FATAL, BAD_ALLOC); }
// 
//   // Create new thread-streams and initialize using the leapfrog method
//   for (Int t=0; t<m_nthreads; ++t) {
//     newStream(&m_stream[t], brng, m_seed);
//     leapfrogStream(m_stream[t], t, m_nthreads);
//   }
// }
// 
// MKLRandom::~MKLRandom()
// //******************************************************************************
// //  Destructor
// //! \details Destroy random number generator thread-streams and tables
// //! \author  J. Bakosi
// //******************************************************************************
// {
//   // Free thread-streams (used for a few at a time)
//   for (int t=0; t<m_nthreads; ++t) {
//     if (m_stream[t] != nullptr &&
//         vslDeleteStream(&m_stream[t]) != VSL_STATUS_OK)
//       cerr << "WARNING: Failed to delete MKL VSL stream" << endl;
//   }
//   // Free array of thread-stream pointers (used for a few at a time)
//   if (m_stream) delete [] m_stream;
// 
// }
