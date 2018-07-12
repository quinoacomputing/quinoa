// *****************************************************************************
/*!
  \file      src/Inciter/DistFCT.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Charm++ module interface for distributed flux-corrected transport
  \details   Charm++ module interface file for asynchronous distributed
    flux-corrected transport (FCT).

    There are a potentially large number of DistFCT Charm++ chares created by
    Transporter. Each DistFCT gets a chunk of the full load (part of the mesh)
    and performs flux-corrected transport.

    The implementation uses the Charm++ runtime system and is fully
    asynchronous, overlapping computation and communication. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality. The high-level
    overview of the algorithm structure and how it interfaces with Charm++ is
    discussed in the Charm++ interface file src/Inciter/distfct.ci.

    #### Call graph ####
    The following is a directed acyclic graph (DAG) that outlines the
    asynchronous algorithm implemented in this class The detailed discussion of
    the algorithm is given in the Charm++ interface file transporter.ci, which
    also repeats the graph below using ASCII graphics. On the DAG orange
    fills denote global synchronization points that contain or eventually lead
    to global reductions. Dashed lines are potential shortcuts that allow
    jumping over some of the task-graph under some circumstances or optional
    code paths (taken, e.g., only in DEBUG mode). See the detailed discussion in
    distfct.ci.
    \dot
    digraph "DistFCT SDAG" {
      rankdir="LR";
      node [shape=record, fontname=Helvetica, fontsize=10];
      Upd [ label="Upd" tooltip="update high-order solution"
                 style="solid"
                URL="\ref tk::Solver::updateSol"];
      LowUpd [ label="LowUpd" tooltip="update low-order solution"
               style="solid"
               URL="\ref tk::Solver::updateLowol"];
      OwnAEC [ label="OwnAEC"
               tooltip="own contributions to the antidiffusive element
                        contributions computed"
               URL="\ref inciter::DistFCT::aec"];
      ComAEC [ label="ComAEC"
               tooltip="contributions to the antidiffusive element contributions
                        communicated"
               URL="\ref inciter::DistFCT::comaec"];
      OwnALW [ label="OwnALW"
               tooltip="own contributions to the maximum and minimum unknowns of
                        elements surrounding nodes computed"
               URL="\ref inciter::DistFCT::alw"];
      ComALW [ label="ComALW"
               tooltip="contributions to the the maximum and minimum unknowns of
                        elements surrounding nodes communicated"
               URL="\ref inciter::DistFCT::comalw"];
      Ver [ label="Ver" tooltip="verify antidiffusive element contributions"
            URL="\ref inciter::DistFCT::verify"];
      OwnLim [ label="OwnLim"
               tooltip="compute limited antidiffusive element contributions"
               URL="\ref inciter::DistFCT::lim"];
      ComLim [ label="ComLim"
               tooltip="contributions to the limited antidiffusive element
                        contributions communicated"
               URL="\ref inciter::DistFCT::comlim"];
      Apply [ label="Apply"
              tooltip="apply limited antidiffusive element contributions"
              URL="\ref inciter::DistFCT::limit"];
      OwnAEC -> Ver [ style="dashed" ];
      OwnALW -> Ver [ style="dashed" ];
      Upd -> OwnAEC [ style="solid" ];
      Upd -> ComEC [ style="solid" ];
      LowUpd -> OwnALW [ style="solid" ];
      LowUpd -> ComALW [ style="solid" ];
      OwnAEC -> OwnLim [ style="solid" ];
      ComAEC -> OwnLim [ style="solid" ];
      OwnALW -> OwnLim [ style="solid" ];
      ComALW -> OwnLim [ style="solid" ];
      OwnAEC -> ComLim [ style="solid" ];
      ComAEC -> ComLim [ style="solid" ];
      OwnALW -> ComLim [ style="solid" ];
      ComALW -> ComLim [ style="solid" ];
      OwnLim -> Apply [ style="solid" ];
      ComLim -> Apply [ style="solid" ];
    }
    \enddot
    \include Inciter/distfct.ci
  \see       DistFCT.[Ch] and FluxCorrector.[Ch] for more info.
*/
// *****************************************************************************
#ifndef DistFCT_h
#define DistFCT_h

#include <cstddef>
#include <iosfwd>
#include <utility>
#include <vector>
#include <cstring>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <set>

#include "QuinoaConfig.h"
#include "Types.h"
#include "Fields.h"
#include "DerivedData.h"
#include "VectorReducer.h"
#include "FluxCorrector.h"
#include "Discretization.h"
#include "MatCG.h"
#include "DiagCG.h"
#include "Inciter/InputDeck/InputDeck.h"

#include "NoWarning/distfct.decl.h"

namespace inciter {

//! DistFCT Charm++ chare array used to advance PDEs in time with MatCG+LW+FCT
class DistFCT : public CBase_DistFCT {

  private:
    //! Variant storing one proxy storing the discretization proxy we work with
    using SchemeProxy = boost::variant< CProxy_MatCG, CProxy_DiagCG >;
    //! Variant type listing chare element proxy types of discretization proxy
    using ProxyElem =
      boost::variant< CProxy_MatCG::element_t, CProxy_DiagCG::element_t >;
  
    //! Functor to call the next() member function behind SchemeProxy
    struct Next : boost::static_visitor<> {
      Next( const tk::Fields& a ) : A(a) {}
      template< typename P >
        void operator()( const P& p ) const { p.ckLocal()->next( A ); }
      const tk::Fields& A;
    };
  
  public:
    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wunused-parameter"
      #pragma clang diagnostic ignored "-Wdeprecated-declarations"
    #elif defined(STRICT_GNUC)
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Wunused-parameter"
      #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    #elif defined(__INTEL_COMPILER)
      #pragma warning( push )
      #pragma warning( disable: 1478 )
    #endif
    // Include Charm++ SDAG code. See http://charm.cs.illinois.edu/manuals/html/
    // charm++/manual.html, Sec. "Structured Control Flow: Structured Dagger".
    DistFCT_SDAG_CODE
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #elif defined(STRICT_GNUC)
      #pragma GCC diagnostic pop
    #elif defined(__INTEL_COMPILER)
      #pragma warning( pop )
    #endif

    //! Constructor
    explicit
    DistFCT( const CProxy_Transporter& host,
             int nchare,
             std::size_t nu,
             std::size_t np,
             const std::unordered_map< int, std::vector< std::size_t > >& msum,
             const std::unordered_map< std::size_t, std::size_t >& bid,
             const std::unordered_map< std::size_t, std::size_t >& lid,
             const std::vector< std::size_t >& inpoel );

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! Migrate constructor
    explicit DistFCT( CkMigrateMessage* ) {}
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

    //! Compute lumped mass lhs required for the low order solution
    tk::Fields lump( const Discretization& d );

    //! \brief Compute mass diffusion rhs contribution required for the low
    //!   order solution
    tk::Fields diff( const Discretization& d, const tk::Fields& Un );

    //! Prepare for next time step stage
    void next();

    //! Receive sums of antidiffusive element contributions on chare-boundaries
    void comaec( const std::vector< std::size_t >& gid,
                 const std::vector< std::vector< tk::real > >& P );

    //! \brief Receive contributions to the maxima and minima of unknowns of all
    //!   elements surrounding mesh nodes on chare-boundaries
    void comalw( const std::vector< std::size_t >& gid,
                 const std::vector< std::vector< tk::real > >& Q );

    //! \brief Receive contributions of limited antidiffusive element
    //!   contributions on chare-boundaries
    void comlim( const std::vector< std::size_t >& gid,
                 const std::vector< std::vector< tk::real > >& A );

    //! Compute and sum antidiffusive element contributions (AEC) to mesh nodes
    void aec( const Discretization& d,
              const tk::Fields& dUh,
              const tk::Fields& Un,
              const std::unordered_map< std::size_t,
                      std::vector< std::pair< bool, tk::real > > >& bc );

    //! \brief Compute the maximum and minimum unknowns of all elements
    //!   surrounding nodes
    void alw( const tk::Fields& Un,
              const tk::Fields& Ul,
              const tk::Fields& dUl,
              const SchemeProxy& scheme );

    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) {
      CBase_DistFCT::pup(p);
      p | m_nhsol;
      p | m_nlsol;
      p | m_naec;
      p | m_nalw;
      p | m_nlim;
      p | m_nchare;
      p | m_msum;
      p | m_bid;
      p | m_lid;
      p | m_inpoel;
      p | m_fluxcorrector;
      p | m_p;
      p | m_q;
      p | m_a;
      p | m_pc;
      p | m_qc;
      p | m_ac;
      p | m_ul;
      p | m_dul;
      p | m_du;
      p | m_scheme;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i DistFCT object reference
    friend void operator|( PUP::er& p, DistFCT& i ) { i.pup(p); }
    //@}

  private:
    using ncomp_t = kw::ncomp::info::expect::type;

    //! Transporter (host) proxy
    CProxy_Transporter m_host;
    //! Counter for high order solution nodes updated
    std::size_t m_nhsol;
    //! Counter for low order solution nodes updated
    std::size_t m_nlsol;
    //! \brief Number of chares from which we received antidiffusive element
    //!   contributions on chare boundaries
    std::size_t m_naec;
    //! \brief Number of chares from which we received maximum and minimum
    //!   unknowns of elements surrounding nodes on chare boundaries
    std::size_t m_nalw;
    //! \brief Number of chares from which we received limited antidiffusion
    //!   element contributiones on chare boundaries
    std::size_t m_nlim;
    //! Total number of worker chares
    std::size_t m_nchare;
    //! \brief Global mesh node IDs bordering the mesh chunk held by fellow
    //!   chares associated to their chare IDs
    //! \details msum: mesh chunks surrounding mesh chunks and their neighbor
    //!   points
    //! \note This is a copy. Original in (bound) Discretization
    std::unordered_map< int, std::vector< std::size_t > > m_msum;
    //! \brief Local chare-boundary mesh node IDs at which we receive
    //!   contributions associated to global mesh node IDs of mesh elements we
    //!   contribute to
    //! \note This is a copy. Original in (bound) Discretization
    std::unordered_map< std::size_t, std::size_t > m_bid;
    //! Local mesh node ids associated to the global ones of owned elements
    //! \note This is a copy. Original in (bound) Discretization
    std::unordered_map< std::size_t, std::size_t > m_lid;
    //! Mesh connectivity of our chunk of the mesh
    //! \note This is a copy. Original in (bound) Discretization
    std::vector< std::size_t > m_inpoel;
    //! Flux corrector performing FCT
    FluxCorrector m_fluxcorrector;
    //! Flux-corrected transport data structures
    tk::Fields m_p, m_q, m_a;
    //! Receive buffers for FCT
    std::vector< std::vector< tk::real > > m_pc, m_qc, m_ac;
    //! Pointer to low order solution vector and increment
    //! \note These are copies. Original in (bound) Discretization
    tk::Fields m_ul, m_dul, m_du;
    //! Variant storing the discretization scheme class we interoperate with
    SchemeProxy m_scheme;

    //! \brief Verify antidiffusive element contributions up to linear solver
    //!   convergence
    void verify();

    //! Compute the limited antidiffusive element contributions
    void lim();

    //! Apply limited antidiffusive element contributions
    void apply();
};

} // inciter::

#endif // DistFCT_h
