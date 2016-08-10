// *****************************************************************************
/*!
  \file      src/IO/ParticleWriter.h
  \author    J. Bakosi
  \date      Tue 09 Aug 2016 08:08:40 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Charm++ group for outputing particle data to file via H5Part
  \details   Charm++ group for outputing particle data to file via H5Part in
     parallel using MPI-IO.
*/
// *****************************************************************************
#ifndef ParticleWriter_h
#define ParticleWriter_h

#include "Exception.h"
#include "H5PartWriter.h"

#include "NoWarning/particlewriter.decl.h"
#include "NoWarning/conductor.decl.h"

namespace tk {

//! \brief Charm++ group used to output particle data to file in parallel using
//!   H5Part and MPI-IO
template< class HostProxy >
class ParticleWriter : public CBase_ParticleWriter< HostProxy > {

  public:
    //! Constructor
    //! \param[in] host Host proxy
    //! \param[in] filename Filename of particle output file
    explicit ParticleWriter( const HostProxy& host,
                             const std::string& filename ) :
      m_host( host ),
      m_writer( filename ),
      m_timestamped( 0 ),
      m_npar( 0 ),
      m_x(),
      m_y(),
      m_z() {}

    //! Chares contribute their number of particles they will output on my PE
    //! \param[in] npar Number of particles will be contributed
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void npar( std::size_t npar ) { m_npar += npar; }

    //! Write a new time stamp to the particle file
    //! \param[in] it Iteration number
    //! \details Since potentially this function may be called mutliple times
    //!   (by multiple chare array elements on our PE, we make sure the
    //!   function call to output the new time stamp is only called once per PE.
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    //! \author J. Bakosi
    void writeTimeStamp( uint64_t it ) {
      if (m_npar > 0 && !m_timestamped) m_writer.writeTimeStamp( it, m_npar );
      m_timestamped = true;
    }

    //! Receive, buffer, and write particle coordinates to file
    //! \param[in] x X coordinates of particles
    //! \param[in] y Y coordinates of particles
    //! \param[in] z Z coordinates of particles
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    //! \author J. Bakosi
    void writeCoords( const std::vector< tk::real >& x,
                      const std::vector< tk::real >& y,
                      const std::vector< tk::real >& z )
    {
      if (m_npar == 0) { signal2host_outcomplete( m_host ); return; }
      Assert( m_timestamped, "Outputing the time stamp must precede the "
                             "particle coordinates" );
      Assert( x.size() == y.size() && y.size() == z.size(),
              "Particle coordinates array sizes mismatch" );
      // buffer up coordinates
      m_x.insert( end(m_x), begin(x), end(x) );
      m_y.insert( end(m_y), begin(y), end(y) );
      m_z.insert( end(m_z), begin(z), end(z) );
      // if received from all chares on my PE, write to file
      if (m_x.size() == m_npar) {
        m_writer.writeCoords( m_x, m_y, m_z );
        signal2host_outcomplete( m_host );
        m_x.clear();        // prepare for next step
        m_y.clear();
        m_z.clear();
        m_npar = 0;
        m_timestamped = false;
      }
    }

  private:
    HostProxy m_host;
    tk::H5PartWriter m_writer;      //!< Particle file format writer
    bool m_timestamped;  //!< Whether time stamp has been written out on this PE
    uint64_t m_npar;               //!< Number of particles to be written
    std::vector< tk::real > m_x;   //!< Buffer collecting x coordinates
    std::vector< tk::real > m_y;   //!< Buffer collecting y coordinates
    std::vector< tk::real > m_z;   //!< Buffer collecting z coordinates

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wdocumentation"
    #endif
    /** @name Host signal calls
      * \brief These functions signal back to the host via a global reduction
      *   originating from each PE branch
      * \details Singal calls contribute to a reduction on all branches (PEs)
      *   of ParticleWriter to the host, e.g., inciter::CProxy_Conductor, given
      *   by the template argument HostProxy. The signal functions are overloads
      *   on the specialization, e.g., inciter::CProxy_Conductor, of the
      *   ParticleWriter template. They create Charm++ reduction targets via
      *   creating a callback that invokes the typed reduction client, where
      *   host is the proxy on which the reduction target method, given by the
      *   string followed by "redn_wrapper_", e.g., parcomplete(), is called
      *   upon completion of the reduction.
      *
      *   Note that we do not use Charm++'s CkReductionTarget macro here,
      *   but instead explicitly generate the code that that macro would
      *   generate. To explain why, here is Charm++'s CkReductionTarget macro's
      *   definition, given in ckreduction.h:
      *   \code{.cpp}
      *      #define CkReductionTarget(me, method) \
      *        CkIndex_##me::redn_wrapper_##method(NULL)
      *   \endcode
      *   This macro takes arguments 'me' (a class name) and 'method' a member
      *   function of class 'me' and generates the call
      *   'CkIndex_<class>::redn_wrapper_<method>(NULL)'. With the overloads the
      *   signal2* functions generate, we do the above macro's job for
      *   ParticleWriter specialized by HostProxy, hard-coded here, as well its
      *   reduction target. This is required since
      *    * Charm++'s CkReductionTarget macro's preprocessing happens earlier
      *      than type resolution and the string of the template argument would
      *      be substituted instead of the type specialized (which is not what
      *      we want here), and
      *    * the template argument class, e.g, CProxy_Conductor, is in a
      *      namespace different than that of ParticleWriter. When a new class
      *      is used to specialize ParticleWriter, the compiler will alert that
      *      a new overload needs to be defined.
      *
      * \note This simplifies client-code, e.g., inciter::Conductor, which now
      *   requires no explicit book-keeping with counters, etc. Also a reduction
      *   (instead of a direct call to the host) better utilizes the
      *   communication network as computational nodes can send their aggregated
      *   contribution to other nodes on a network instead of all chares sending
      *   their (smaller) contributions to the same host, (hopefully)
      *   implemented using a tree among the PEs.
      * \see http://charm.cs.illinois.edu/manuals/html/charm++/manual.html,
      *   Sections "Processor-Aware Chare Collections" and "Chare Arrays".
      * */
    ///@{
    //! \brief Signal back to host that the initialization of the row indices of
    //!   the linear system is complete
    void signal2host_outcomplete( const inciter::CProxy_Conductor& host )
    {
      using inciter::CkIndex_Conductor;
      Group::contribute(
        CkCallback( CkIndex_Conductor::redn_wrapper_outcomplete(NULL), host ) );
    }
    ///@}
    #if defined(__clang__)
      #pragma GCC diagnostic pop
    #endif
};

} // tk::

#define CK_TEMPLATES_ONLY
#include "NoWarning/particlewriter.def.h"
#undef CK_TEMPLATES_ONLY

#endif // ParticleWriter_h
