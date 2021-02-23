// *****************************************************************************
/*!
  \file      src/Base/Progress.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Simple progress indicator
  \details   Simple progress indicator.
*/
// *****************************************************************************
#ifndef Progress_h
#define Progress_h

#include <array>
#include <functional>

#include "TaggedTuple.hpp"
#include "Print.hpp"
#include "Tags.hpp"

namespace tk {

//! Simple progress class for outputing progress indicators during a task
//! \details This is a helper class to abstract away the details of using
//!   Print::progress() used to output progress reports to the screen during
//!   a task consisting of multiple sub-tasks happening at the same time. The
//!   template argument is a compile-time integer which is the number of
//!   independent sub-tasks the progress indicator receives messages for and
//!   counts them independtly toward multiple maxima.
template< std::size_t N >
class Progress {

  public:
    //! Constructor
    //! \param[in] feedback Whether to send sub-task feedback to host    
    //! \param[in] prefix Strings to output prefixing the progress report
    //! \param[in] legend Legend for each prefix to output at start
    //! \param[in] max Array of integers equaling the max number of items to be
    //!   expected per sub-task
    explicit Progress( bool feedback,
                       const std::array< std::string, N >& prefix,
                       const std::array< std::string, N >& legend,
                       std::array< int, N >&& max = std::array< int, N >() )
    : m_feedback( feedback ),
      m_prefix( std::move(prefix) ),
      m_legend( std::move(legend) ),
      m_finished( false ),
      m_progress_size( 0 ),
      m_max( std::move(max) )
    {
      m_done.fill( 0 ); 
    }

   //! Start counting sub-tasks outputing an intial task message
   //! \param[in] print Pretty printer object to use for printing progress
   //! \param[in] msg Message to output to screen. This message should be
   //!   descriptive of all the sub-tasks we are responsible for. I.e., this
   //!   is usually a list of multiple sub-tasks happening at the same time.
   //!   Appending to msg we also output the legend of subtasks in parentheses.
   void start( const Print& print, const std::string& msg ) {
     std::string legend;
     if (m_feedback) {
       legend.append( " (" );
       for (const auto& l : m_legend) legend.append( l + ", " );
       legend.pop_back();
       legend.pop_back();
       legend.append( ")" );
     }
     legend.append( " ..." );
     print.diagstart( msg + legend );
     m_progress_size = 0;
   }

   //! \brief Start counting sub-tasks outputing an intial task message and set
   //!   max number of items to be expected per sub-task
   //! \param[in] print Pretty printer object to use for printing progress
   //! \param[in] msg Message to output to screen. This message should be
   //!   descriptive of all the sub-tasks we are responsible for. I.e., this
   //!   is usually a list of multiple sub-tasks happening at the same time.
   //! \param[in] max Array of integers equaling the max number of items to be
   //!   expected per sub-task
   //! \details This function can be used to do the same as start( msg ) and
   //!   update/reset the max number of items per sub-task in case they are not
   //!   all yet available when the constructor is called.
   void start( const Print& print,
               const std::string& msg,
               std::array< int, N >&& max )
   {
     m_max = std::move(max);
     start( print, msg );
   }

   //! Receive an update to a sub-task counter and update progress report
   //! \param[in] print Pretty printer object to use for printing progress
   //! \details The template argument indexes the sub-task. A compile-time
   //!   assert emits an error in case it is out of bounds.
   template< std::size_t i > void inc( const Print& print ) {
     static_assert( i < N, "Indexing out of bounds" );
     ++m_done[i];
     if (!m_finished) report( print );
   }

   //! Finish progress report updating it one last time
   //! \param[in] print Pretty printer object to use for printing progress
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
   void end( const Print& print ) {
     m_finished = true;
     m_done = m_max;
     report( print );
     print.diagend( "done" );
   }

  private:
    bool m_feedback;            //!< Whether to send sub-task feedback to host
    const std::array< std::string, N > m_prefix;        //!< Sub-task prefixes
    const std::array< std::string, N > m_legend;        //!< Sub-task legend
    bool m_finished;            //!< Whether task has finished
    std::size_t m_progress_size;//!< Size of previous progress report
    std::array< int, N > m_max; //!< Max number of items per sub-task
    std::array< int, N > m_done;//!< Number of items done per sub-task

   //! Output progress report to screen
   //! \param[in] print Pretty printer object to use for printing progress
   //! \details This output contains a status on each of the multiple sub-task
   //!   counters as they all work towards their respective maxima.
   //! \see Print::progress()
   void report( const Print& print ) {
     if (m_feedback)
       print.progress< N >( m_prefix, m_done, m_max, m_progress_size );
   }
};

} // tk::

#endif // Progress_h
