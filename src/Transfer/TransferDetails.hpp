// *****************************************************************************
/*!
  \file      src/Transfer/TransferDetails.hpp
  \copyright 2020 Charmworks, Inc.
             All rights reserved. See the LICENSE file for details.
  \brief     Chare class declaration for mesh transfer workers holding part of a
    mesh
  \details   Chare class declaration for mesh transfer workers holding part of a
    mesh.
*/
// *****************************************************************************
#ifndef TransferDetails_h
#define TransferDetails_h

#include "Types.hpp"
#include "PUPUtil.hpp"
#include "UnsMesh.hpp"
#include "CommMap.hpp"
#include "Fields.hpp"

#include "NoWarning/transferdetails.decl.h"

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

//! TransferDetails chare array holding part of a mesh
class TransferDetails : public CBase_TransferDetails {

  public:
    //! Constructor
    explicit TransferDetails(
      CkArrayID p,
      MeshData d,
      CollideHandle ch,
      CkCallback cb );

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! Migrate constructor
    // cppcheck-suppress uninitMemberVar
    explicit TransferDetails( CkMigrateMessage* ) {}
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
    void processCollisions( CProxy_TransferDetails proxy,
                            int nchare,
                            int offset,
                            int nColls,
                            Collision* colls );

    //! Identify actual collisions in the source mesh
    void determineActualCollisions( CProxy_TransferDetails proxy,
                                    int index,
                                    int nColls,
                                    PotentialCollision* colls ) const;

    //! Transfer the interpolated solution data back to destination mesh
    void transferSolution( const std::vector< SolutionData >& soln );

  private:
    //! The ID of my first chunk (used for collision detection library)
    int m_firstchunk;
    //! Pointer to element connectivity
    std::vector< std::size_t >* m_inpoel;
    //! Pointer to point coordinates
    tk::UnsMesh::Coords* m_coord;
    //! Pointer to solution in mesh nodes
    tk::Fields* m_u;
    //! Collide handle
    CollideHandle collideHandle;

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
};

} // exam2m::

#endif // TransferDetails_h
