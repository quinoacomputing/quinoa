// *****************************************************************************
/*!
  \file      src/Transfer/Worker.hpp
  \copyright 2020 Charmworks, Inc.
             All rights reserved. See the LICENSE file for details.
  \brief     Chare class declaration for workers holding part of a mesh
  \details   Chare class declaration for workers holding part of a mesh.
*/
// *****************************************************************************
#ifndef Worker_h
#define Worker_h

#include "Types.hpp"
#include "PUPUtil.hpp"
#include "UnsMesh.hpp"
#include "CommMap.hpp"
#include "Fields.hpp"

#include "NoWarning/worker.decl.h"

namespace exam2m {

class PotentialCollision {
  public:
    std::size_t source_index, dest_index;
    CkVector3d point;
    void pup(PUP::er& p) { p | source_index; p | dest_index; p | point; }
};

class SolutionData {
  public:
    std::size_t dest_index;
    std::vector< tk::real > solution;
    void pup(PUP::er& p) { p | dest_index; p | solution; }
};

//! Worker chare array holding part of a mesh
class Worker : public CBase_Worker {

  public:
    //! Constructor
    explicit Worker( CkArrayID p, MeshData d, CkCallback cb );

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! Migrate constructor
    // cppcheck-suppress uninitMemberVar
    explicit Worker( CkMigrateMessage* ) {}
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

    //! Set the source mesh data
    void setSourceTets( std::vector< std::size_t>* inpoel,
                        tk::UnsMesh::Coords* coords,
                        const tk::Fields& u );

    //! Set the destination mesh data
    void setDestPoints( tk::UnsMesh::Coords* coords,
                        tk::Fields& u,
                        CkCallback cb );

    //! Process potential collisions in the destination mesh
    void processCollisions( CProxy_Worker proxy,
                            int nchare,
                            int offset,
                            int nColls,
                            Collision* colls );

    //! Identify actual collisions in the source mesh
    void determineActualCollisions( CProxy_Worker proxy,
                                    int index,
                                    int nColls,
                                    PotentialCollision* colls ) const;

    //! Transfer the interpolated solution data back to destination mesh
    void transferSolution( const std::vector< SolutionData >& soln );

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) override {
      p | m_firstchunk;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i Worker object reference
    friend void operator|( PUP::er& p, Worker& i ) { i.pup(p); }
    //@}

  private:
    //! The ID of my first chunk (used for collision detection library)
    int m_firstchunk;
    //! Pointer to element connectivity
    std::vector< std::size_t >* m_inpoel;
    //! Pointer to point coordinates
    tk::UnsMesh::Coords* m_coord;
    //! Pointer to solution in mesh nodes
    tk::Fields* m_u;

    //! The number of messages sent by the dest mesh
    int m_numsent;
    //! The number of messages received by the dest mesh
    int m_numreceived;
    //! Called once the transfer is complete (m_numsent == m_numreceived)
    CkCallback m_donecb;

    //! Initialize dest mesh solution with background data
    void background();

    //! Contribute vertex information to the collsion detection library
    void collideVertices();

    //! Contribute tet information to the collision detection library
    void collideTets() const;

    //! Determine if a point is in a tet
    bool intet(const CkVector3d &point,
               std::size_t e,
               std::array< real, 4 >& N) const;
};

} // exam2m::

#endif // Worker_h
