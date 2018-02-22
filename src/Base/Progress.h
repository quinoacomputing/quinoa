// *****************************************************************************
/*!
  \file      src/Base/Progress.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Simple progress indicator
  \details   Simple progress indicator.
*/
// *****************************************************************************
#ifndef Progress_h
#define Progress_h

#include <array>
#include <functional>

#include "TaggedTuple.h"
#include "Print.h"
#include "Tags.h"

namespace tk {

//! Simple progress class for outputing progress indicators during a task
//! \details This is a helper class to abstract away the details of using
//!   tk::Print::progress() used to output progress reports to the screen during
//!   a task consisting of multiple sub-tasks happening at the same time. The
//!   template argument is a compile-time integer which is the number of
//!   independent sub-tasks the progress indicator receives messages for and
//!   counts them independtly toward multiple maxima.
template< std::size_t N >
class Progress {

  public:
    //! Constructor
    //! \param[in] print Pretty printer object to use for printing progress
    //! \param[in] feedback Whether to send sub-task feedback to host    
    //! \param[in] prefix Strings to output prefixing the progress report
    //! \param[in] max Array of integers equaling the max number of items to be
    //!   expected per sub-task
    //! \details Note that prefix and max are swallowed, i.e., moved in,
    //!   treating them as rvalue references.
    explicit Progress( const tk::Print& print,
                       bool feedback,
                       std::array< std::string, N >&& prefix,
                       std::array< int, N >&& max = std::array< int, N >() )
    : m_print( print ),
      m_feedback( feedback ),
      m_prefix( std::move(prefix) ),
      m_finished( false ),
      m_progress_size( 0 ),
      m_max( std::move(max) )
    {
      m_done.fill( 0 ); 
    }

   //! Start counting sub-tasks outputing an intial task message
   //! \param[in] msg Message to output to screen. This message should be
   //!   descriptive of all the sub-tasks we are responsible for. I.e., this
   //!   is usually a list of multiple sub-tasks happening at the same time.
   void start( const std::string& msg ) {
     m_print.diagstart( msg );
     m_progress_size = 0;
   }

   //! \brief Start counting sub-tasks outputing an intial task message and set
   //!   max number of items to be expected per sub-task
   //! \param[in] msg Message to output to screen. This message should be
   //!   descriptive of all the sub-tasks we are responsible for. I.e., this
   //!   is usually a list of multiple sub-tasks happening at the same time.
   //! \param[in] max Array of integers equaling the max number of items to be
   //!   expected per sub-task
   //! \details This function can be used to do the same as start( msg ) and
   //!   update/reset the max number of items per sub-task in case they are not
   //!   all yet available when the constructor is called.
   void start( const std::string& msg, std::array< int, N >&& max ) {
     m_max = std::move(max);
     start( msg );
   }

   //! Receive an update to a sub-task counter and update progress report
   //! \details The template argument indexes the sub-task. A compile-time
   //!   assert emits an error in case it is out of bounds.
   template< std::size_t i > void inc() {
     static_assert( i < N, "Indexing out of bounds" );
     ++m_done[i];
     if (!m_finished) report();
   }

   //! Finish progress report updating it one last time
   //! \details When this function is called, all sub-tasks are assumed to be
   //!   finished, i.e., assumed to have reached their respective maximum
   //!   values. Thus we update our 'done' array to be equal to 'max' and output
   //!   the progress report one final time before outputing 'done'. When this
   //!   function is called it is possible that not all sub-task counters have
   //!   reached their maximum, which can happen if the final reduction (if
   //!   exists), signaling the absolute end of a task (consisting of multiple
   //!   sub-tasks we count counters for), is scheduled before (or is faster)
   //!   than as the individual sub-task counting messages arrive. Even if that
   //!   is the case, this call "officially" finishes all sub-tasks, and outputs
   //!   the progress report using the max values for all sub-tasks to leave a
   //!   consistent screen output finishing the task.
   void end() {
     m_finished = true;
     m_done = m_max;
     report();
     m_print.diagend( "done" );
   }

  private:
    const tk::Print& m_print;   //!< Pretty printer to use for screen output
    bool m_feedback;            //!< Whether to send sub-task feedback to host
    const std::array< std::string, N > m_prefix;        //!< Sub-task prefixes
    bool m_finished;            //!< Whether task has finished
    std::size_t m_progress_size;//!< Size of previous progress report
    std::array< int, N > m_max; //!< Max number of items per sub-task
    std::array< int, N > m_done;//!< Number of items done per sub-task

   //! Output progress report to screen
   //! \details This output contains a status on each of the multiple sub-task
   //!   counters as they all work towards their respective maxima.
   //! \see tk::Print::progress()
   void report() {
     if (m_feedback)
       m_print.progress< N >( m_prefix, m_done, m_max, m_progress_size );
   }
};

} // tk::

#endif // Progress_h
