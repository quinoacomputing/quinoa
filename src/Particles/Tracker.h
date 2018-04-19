// *****************************************************************************
/*!
  \file      src/Particles/Tracker.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Tracker tracks Lagrangian particles in physical space
  \details   Tracker tracks Lagrangian particles in physical space. It works on
    a chunk of the Eulerian mesh, and tracks particles in elements and across
    mesh chunks held by different Charm++ chares.
*/
// *****************************************************************************
#ifndef Tracker_h
#define Tracker_h

#include <vector>
#include <array>
#include <set>
#include <unordered_map>

#include "NoWarning/pup.h"

#include "Keywords.h"
#include "Particles.h"
#include "DerivedData.h"
#include "ParticleWriter.h"
#include "ContainerUtil.h"
#include "PUPUtil.h"

#include "NoWarning/transporter.decl.h"

namespace tk {

//! Tracker advances Lagrangian particles in state space
class Tracker {

  public:
    //! Constructor
    //! \param[in] npar Number of particles per mesh element
    //! \param[in] inpoel Mesh element connectivity
    //! \param[in] feedback Whether to send sub-task feedback to host
    explicit Tracker( bool feedback = false,
                      std::size_t npar = 0,
                      const std::vector< std::size_t >& inpoel = {} ) :
      m_particles( npar * inpoel.size()/4, 3 ), // only the 3 spatial components
      m_elp( m_particles.nunk() ),
      m_parmiss(),
      m_parelse(),
      m_nchpar( 0 ),
      m_esupel( tk::genEsupel( inpoel, 4, tk::genEsup(inpoel,4) ) ),
      m_feedback( feedback )
    {}

    //! Generate particles to each of our mesh cells
    void
    genpar( const std::array< std::vector< tk::real >, 3 >& coord,
            const std::vector< std::size_t >& inpoel,
            std::size_t nchare,
            int chid );

    //! Output number of particles we will write to file in this step
    //! \param[in] hostproxy Charm++ host proxy to which address reductions
    //! \param[in] pw Charm++ particle writer proxy
    //! \param[in] array Charm++ array object pointer of the holder class
    template< class HostProxy, class ParticleWriterProxy, class ChareArray >
    void writeParticles( HostProxy& hostproxy,
                         const ParticleWriterProxy& pw,
                         ChareArray* const array )
    {
      // Send number of partciles we will contribute to particle writer
      pw.ckLocalBranch()->npar( m_particles.nunk() );
      // Tell the host that we are done with sending our number of particles
      signal2host_nparcomplete( hostproxy, array );
    }

    //! Output particles fields to file
    //! \param[in] pw Charm++ particle writer proxy
    //! \param[in] it Iteration count
    //! \param[in] nchare Number of chares that contribute
    template< class ParticleWriterProxy >
    void doWriteParticles( const ParticleWriterProxy& pw,
                           uint64_t it,
                           std::size_t nchare )
    {
      pw.ckLocalBranch()->writeCoords( nchare,
                                       it,
                                       m_particles.extract(0,0),
                                       m_particles.extract(1,0),
                                       m_particles.extract(2,0) );
    }

    //! Advance particle based on velocity from mesh cell
    //! \param[in] array Charm++ array object pointer of the holder class
    //! \param[in] i Particle index
    //! \param[in] e Mesh element index where the particle currently resides
    //! \param[in] dt Time step size
    //! \param[in] Np Four finite-element shapefunctions evaluated (as a result
    //!   of the particle search) at the particle location in element e
    template< class ChareArray >
    void advanceParticle( ChareArray* const array,
                          std::size_t i,
                          std::size_t e,
                          tk::real dt,
                          const std::array< tk::real, 4 >& Np )
    {
      // Extract the transport velocity at nodes
      auto v = array->velocity( e );
      // Advance particle coordinates using the interpolated velocity
      m_particles(i,0,0) +=
        dt*(Np[0]*v[0][0] + Np[1]*v[0][1] + Np[2]*v[0][2] + Np[3]*v[0][3]);
      m_particles(i,1,0) +=
        dt*(Np[0]*v[1][0] + Np[1]*v[1][1] + Np[2]*v[1][2] + Np[3]*v[1][3]);
      m_particles(i,2,0) +=
        dt*(Np[0]*v[2][0] + Np[1]*v[2][1] + Np[2]*v[2][2] + Np[3]*v[2][3]);
      // Apply boundary conditions to particle
      applyParBC( i );
    }

    //! Advance our particles and initiate search for their new mesh cells
    //! \param[in] hostproxy Charm++ host proxy to which address reductions
    //! \param[in] arrayProxy Charm++ array proxy to which address
    //!   point-to-point communications (this is the proxy that holds us)
    //! \param[in] coord Mesh node coordinates
    //! \param[in] inpoel Mesh element connectivity
    //! \param[in] msum Mesh chunks surrounding mesh chunks; we only use the
    //!   keys of this container to address fellow Charm++ chare array elements
    //!   via the arrayProxy
    //! \param[in] chid Charm++ array index (thisIndex of the holder class)
    //! \param[in] array Charm++ array object pointer of the holder class
    //! \param[in] dt Time step size
    template< class HostProxy, class ChareArrayProxy, class ChareArray >
    void track( HostProxy& hostproxy,
                const ChareArrayProxy& arrayProxy,
                const std::array< std::vector< tk::real >, 3 >& coord,
                const std::vector< std::size_t >& inpoel,
                const std::unordered_map<int, std::vector< std::size_t >> msum,
                int chid,
                ChareArray* const array,
                tk::real dt )
    {
      // Lambda to attempt to find and advance particle i in element e. Returns
      // true if the particle was found (and advanced), false if was not found.
      std::array< tk::real, 4 > N;
      auto adv = [ this, array, &N, &dt, &coord, &inpoel ]
                ( std::size_t i, std::size_t e ) -> bool
      {
        if (this->parinel( coord, inpoel, i, e, N )) {
          advanceParticle( array, i, e, dt, N );
          return true;
        }
        return false;
      };
      // Search cells of our mesh chunk for all particles
      for (std::size_t i=0; i<m_particles.nunk(); ++i) {
        // Get element ID where particle i has last been seen
        bool found = adv( i, m_elp[i] );
        // Next search in the elements surroundings the points of the element
        // where the particle has last been found
        if (!found) {
          auto last = m_esupel.second[m_elp[i]+1];
          for (auto j=m_esupel.second[m_elp[i]]+1; j<=last; ++j) {
            found = adv( i, m_esupel.first[j] );
            if (found) j = last+1;  // search for next particle
          }
        }
        // Next search all cells in our chunk of the mesh
        if (!found) {
          for (std::size_t e=0; e<inpoel.size()/4; ++e) {
            found = adv( i, e );
            if (found) e = inpoel.size()/4;  // search for next particle
          }
          // If the particle still has not been found, it left our chunk of the
          // mesh, mark as missing (will initiate communication to find it)
          if (!found) m_parmiss.insert( i );
        }
      }
      // If we have no missing particles, we are done, if we do, send out
      // requests to find them to those ChareArray chares which we neighbor
      // mesh cells with.
      if (m_parmiss.empty()) {
        signal2host_parcomcomplete( hostproxy, array );
      } else {
        std::vector< std::vector< tk::real > > pexp;
        for (auto i : m_parmiss) pexp.push_back( m_particles[i] );
        m_nchpar = 0;
        std::vector< std::size_t > miss( begin(m_parmiss), end(m_parmiss) );
        for (const auto& n : msum)
          arrayProxy[ n.first ].findpar( chid, miss, pexp );
      }
    }

    //! Find particles missing by the requestor and make those found ours
    //! \param[in] arrayProxy Charm++ array proxy to which address
    //!   point-to-point communications (this is the proxy that holds us)
    //! \param[in] coord Mesh node coordinates
    //! \param[in] inpoel Mesh element connectivity
    //! \param[in] fromch Chare ID the request originates from
    //! \param[in] miss Indices of particles to find
    //! \param[in] ps Particle data associated to those particle indices to find
    template< class ChareArrayProxy >
    void findpar( const ChareArrayProxy& arrayProxy,
                  const std::array< std::vector< tk::real >, 3 >& coord,
                  const std::vector< std::size_t >& inpoel,
                  int fromch,
                  const std::vector< std::size_t >& miss,
                  const std::vector< std::vector< tk::real > >& ps )
    {
      // Try to find particles missing by the requestor and own those found
      auto found = addpar( coord, inpoel, miss, ps );
      // Send the particle indices we found back to the requestor
      arrayProxy[ fromch ].foundpar( found );
    }

    //! Receive particle indices found elsewhere (by fellow neighbors)
    //! \param[in] hostproxy Charm++ host proxy to which address reductions
    //! \param[in] arrayProxy Charm++ array proxy to whose all elements te
    //!   address our desparate broadcast (this is the proxy that holds us)
    //! \param[in] msum Mesh chunks surrounding mesh chunks; we only use the
    //!   keys of this container to address fellow Charm++ chare array elements
    //!   via the arrayProxy
    //! \param[in] array Charm++ array object pointer of the holder class
    //! \param[in] chid Charm++ array index (thisIndex of the holder class)
    //! \param[in] found Indices of particles found
    template< class HostProxy, class ChareArrayProxy, class ChareArray >
    void
    foundpar( HostProxy& hostproxy,
              ChareArrayProxy& arrayProxy,
              const std::unordered_map< int, std::vector< std::size_t > > msum,
              ChareArray* const array,
              int chid,
              const std::vector< std::size_t >& found )
    {
      m_parelse.insert( begin(found), end(found) );
      if (++m_nchpar == msum.size()) {    // if we have heard from all neighbors
        remove( m_parelse );  // delete particles found elsewhere
        // find particle that are still have not been found (by close neighbors)
        std::set< std::size_t > far;
        std::set_difference( begin(m_parmiss), end(m_parmiss),
                             begin(m_parelse), end(m_parelse), 
                             std::inserter( far, begin(far) ) );
        m_parmiss = far;
        // if there are still missing particles (not found by close neighbors we
        // share mesh nodes with), we resort to requesting them to be searched by
        // all holder chares
        if (m_parmiss.empty()) {
          signal2host_parcomcomplete( hostproxy, array );
        } else {
          std::vector< std::vector< tk::real > > pexp;
          for (auto i : m_parmiss) pexp.push_back( m_particles[i] );
          m_nchpar = 0;
          std::vector< std::size_t > miss( begin(m_parmiss), end(m_parmiss) );
          m_parelse.clear();
          arrayProxy.collectpar( chid, miss, pexp ); // broadcast to everyone
        }
      }
    }

    //! Find particles missing by the requestor and make those found ours
    //! \param[in] arrayProxy Charm++ array proxy to which address
    //!   point-to-point communications (this is the proxy that holds us)
    //! \param[in] coord Mesh node coordinates
    //! \param[in] inpoel Mesh element connectivity
    //! \param[in] fromch Chare ID the request originates from
    //! \param[in] miss Indices of particles to find
    //! \param[in] ps Particle data associated to those particle indices to find
    template< class ChareArrayProxy >
    void collectpar( const ChareArrayProxy& arrayProxy,
                     const std::array< std::vector< tk::real >, 3 >& coord,
                     const std::vector< std::size_t >& inpoel,
                     int fromch,
                     const std::vector< std::size_t >& miss,
                     const std::vector< std::vector< tk::real > >& ps )
    {
      // Try to find particles missing by the requestor and own those found
      auto found = addpar( coord, inpoel, miss, ps );
      // Send the particle indices we found back to the requestor
      arrayProxy[ fromch ].collectedpar( found );
    }

    //! Collect particle indices found elsewhere (by far fellows)
    //! \param[in] hostproxy Charm++ host proxy to which address reductions
    //! \param[in] array Charm++ array object pointer of the holder class
    //! \param[in] found Indices of particles found
    //! \param[in] nchare Total number of holder array chares
    template< class HostProxy, class ChareArray >
    void collectedpar( HostProxy& hostproxy,
                       ChareArray* const array,
                       const std::vector< std::size_t >& found,
                       std::size_t nchare )
    {
      // Collect particle indices found elsewhere (by distant neighbors)
      m_parelse.insert( begin(found), end(found) );
      if (++m_nchpar == nchare) {  // if we have heard from everyone
        remove( m_parelse );  // delete particles found elsewhere
        Assert( m_parmiss == m_parelse, "Not all particles have been found" );
        signal2host_parcomcomplete( hostproxy, array );
      }
    }

    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er& p ) {
      p | m_particles;
      p | m_elp;
      p | m_parmiss;
      p | m_parelse;
      p | m_nchpar;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i Tracker object reference
    friend void operator|( PUP::er& p, Tracker& i ) { i.pup(p); }
    //@}

  private:
    //! Particle properties
    tk::Particles m_particles;
    //! Element ID in which a particle has last been found for all particles
    std::vector< std::size_t > m_elp;
    //! Indicies of particles not found here (missing)
    std::set< std::size_t > m_parmiss;
    //! Indicies of particles not found here but found by fellows
    std::set< std::size_t > m_parelse;
    //! Number of chares we received particles from
    std::size_t m_nchpar;
    //! Elements surrounding points of elements of mesh chunk we operate on
    std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
      m_esupel;
    //! Bool that determines whether to send sub-task feedback to host
    bool m_feedback;

    //! Try to find particles and add those found to the list of ours
    std::vector< std::size_t >
    addpar( const std::array< std::vector< tk::real >, 3 >& coord,
            const std::vector< std::size_t >& inpoel,
            const std::vector< std::size_t >& miss,
            const std::vector< std::vector< tk::real > >& ps );

    //! Search particle in a single mesh cell
    bool parinel( const std::array< std::vector< tk::real >, 3 >& coord,
                  const std::vector< std::size_t >& inpoel,
                  std::size_t p,
                  std::size_t e,
                  std::array< tk::real, 4 >& N );

     //! Apply boundary conditions to particles
    void applyParBC( std::size_t i );

    //! Remove a set of particles
    void remove( const std::set< std::size_t >& idx );

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wdocumentation"
    #endif
    /** @name Calls to host
      * \brief These functions signal back to the host via a global reduction
      *   originating from each chare we are held by
      * \details Singal calls contribute to a reduction on all holding chares
      *   to the host, e.g., inciter::CProxy_Transporter, given by the template
      *   argument HostProxy. The signal functions are overloads on the
      *   specialization, e.g., inciter::CProxy_Transporter, of the Tracker
      *   template. They create Charm++ reduction targets via creating a
      *   callback that invokes the typed reduction client, where host is the
      *   proxy on which the reduction target method, given by the string
      *   followed by "redn_wrapper_", e.g., nparcomplete(), is called upon
      *   completion of the reduction.
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
      *   below functions generate, we do the above macro's job for Tracker
      *   specialized by HostProxy, hard-coded here, as well its reduction
      *   target. This is required because
      *    * Charm++'s CkReductionTarget macro's preprocessing happens earlier
      *      than type resolution and the string of the template argument would
      *      be substituted instead of the type specialized (which is not what
      *      we want here), and
      *    * the template argument class, e.g, CProxy_Transporter, is in a
      *      namespace different than that of Tracker. When a new class is
      *      used to specialize Tracker, the compiler will alert that a new
      *      overload needs to be defined.
      *
      * \note This simplifies client-code, e.g., inciter::Transporter, which now
      *   requires no explicit book-keeping with counters, etc. Also a reduction
      *   (instead of a direct call to the host) better utilizes the
      *   communication network as computational nodes can send their aggregated
      *   contribution to other nodes on a network instead of all chares sending
      *   their (smaller) contributions to the same host, implemented using a
      *   tree among the chares and PEs.
      * \see http://charm.cs.illinois.edu/manuals/html/charm++/manual.html,
      *   Sections "Processor-Aware Chare Collections" and "Chare Arrays".
      * */
    ///@{
    //! Tell the host that we are done with sending our number of particles
    //! \param[in] host Host proxy to contribute to
    //! \param[in] array Charm++ array object pointer of the holder class
    template< class ChareArray >
    void signal2host_nparcomplete( const inciter::CProxy_Transporter& host,
                                   ChareArray* array )
    {
      using inciter::CkIndex_Transporter;
      array->contribute(
        CkCallback( CkIndex_Transporter::redn_wrapper_nparcomplete(NULL),
                    host ) );
    }
    //! Tell the host that we are done communicating particles
    //! \param[in] host Host proxy to contribute to
    //! \param[in] array Charm++ array object pointer of the holder class
    template< class ChareArray >
    void signal2host_parcomcomplete( inciter::CProxy_Transporter& host,
                                     ChareArray* array )
    {
      // send progress report to host
      if (m_feedback) host.chtrack();
      m_nchpar = 0;
      m_parmiss.clear();
      m_parelse.clear();
      using inciter::CkIndex_Transporter;
      array->contribute(
        CkCallback( CkIndex_Transporter::redn_wrapper_parcomcomplete(NULL),
                    host ) );
    }
    ///@}
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif
};

} // tk::

#endif // Tracker_h
