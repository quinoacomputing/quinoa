// Controller for the library

#include "Controller.hpp"
#include "Worker.hpp"

#include <cassert>

namespace exam2m {

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wmissing-prototypes"
  #pragma clang diagnostic ignored "-Wmissing-variable-declarations"
#endif

/* readonly */ CProxy_Controller controllerProxy;
//! \brief Charm handle to the collision detection library instance
/* readonly */ CollideHandle collideHandle;

void collisionHandler( [[maybe_unused]] void *param,
                        int nColl,
                        Collision *colls )
{
  controllerProxy.ckLocalBranch()->distributeCollisions( nColl, colls );
}

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
  #pragma clang diagnostic ignored "-Wunused-private-field"
  #pragma clang diagnostic ignored "-Wvla"
  #pragma clang diagnostic ignored "-Wvla-extension"
#endif

void addMesh(CkArrayID p, int elem, CkCallback cb) {
  controllerProxy[0].addMesh(p, elem, cb);
}

void setSourceTets(CkArrayID p, int index, std::vector< std::size_t >* inpoel, tk::UnsMesh::Coords* coords, const tk::Fields& u) {
  controllerProxy.ckLocalBranch()->setSourceTets(p, index, inpoel, coords, u);
}

void setDestPoints(CkArrayID p, int index, tk::UnsMesh::Coords* coords, tk::Fields& u, CkCallback cb) {
  controllerProxy.ckLocalBranch()->setDestPoints(p, index, coords, u, cb);
}

LibMain::LibMain(CkArgMsg* msg) {
  delete msg;
  controllerProxy = CProxy_Controller::ckNew();

  // TODO: Need to make sure this is actually correct
  CollideGrid3d gridMap(CkVector3d(0, 0, 0),CkVector3d(2, 100, 2));
  collideHandle = CollideCreate(gridMap,
      CollideSerialClient(collisionHandler, 0));
}

Controller::Controller() : current_chunk(0) {}

void Controller::addMesh(CkArrayID p, int elem, CkCallback cb) {
  auto id = static_cast<std::size_t>(CkGroupID(p).idx);
  if (proxyMap.count(id) == 0) {
    CkArrayOptions opts;
    opts.bindTo(p);
    opts.setNumInitial(elem);
    MeshData mesh;
    mesh.m_nchare = elem;
    mesh.m_firstchunk = current_chunk;
    mesh.m_proxy = CProxy_Worker::ckNew(p, mesh, cb, opts);
    proxyMap[id] = mesh;
    current_chunk += elem;
  } else {
    CkAbort("Uhoh...\n");
  }
}
void Controller::setMesh( CkArrayID p, MeshData d ) {
  proxyMap[static_cast<std::size_t>(CkGroupID(p).idx)] = d;
}

void Controller::setDestPoints(CkArrayID p, int index, tk::UnsMesh::Coords* coords, tk::Fields& u, CkCallback cb) {
  m_destmesh = static_cast<std::size_t>(CkGroupID(p).idx);
  Worker* w = proxyMap[m_destmesh].m_proxy[index].ckLocal();
  assert(w);
  w->setDestPoints(coords, u, cb);
}

void Controller::setSourceTets(CkArrayID p, int index, std::vector< std::size_t >* inpoel, tk::UnsMesh::Coords* coords, const tk::Fields& u) {
  m_sourcemesh = static_cast<std::size_t>(CkGroupID(p).idx);
  Worker* w = proxyMap[m_sourcemesh].m_proxy[index].ckLocal();
  assert(w);
  w->setSourceTets(inpoel, coords, u);
}

void
Controller::distributeCollisions(int nColl, Collision* colls)
// *****************************************************************************
//  Called when all potential collisions have been found, and now need to be
//  destributed to the chares in the destination mesh to determine actual
//  collisions.
//! \param[in] nColl Number of potential collisions found
//! \param[in] colls The list of potential collisions
// *****************************************************************************
{
  CkPrintf("Collisions found: %i\n", nColl);
  auto first = static_cast<int>(proxyMap[m_destmesh].m_firstchunk);
  auto nchare = static_cast<int>(proxyMap[m_destmesh].m_nchare);
  std::vector<Collision> separated[nchare];

  // Separate collisions based on the destination mesh chare they belong to
  for (int i = 0; i < nColl; i++) {
    if (colls[i].A.chunk >= first && colls[i].A.chunk < first + nchare) {
      separated[static_cast<std::size_t>(colls[i].A.chunk - first)].push_back(colls[i]);
    } else {
      separated[static_cast<std::size_t>(colls[i].B.chunk - first)].push_back(colls[i]);
    }
  }

  // Send out each list to the destination chares for further processing
  for (int i = 0; i < nchare; i++) {
    CkPrintf("Dest mesh chunk %i has %lu\n", i, separated[i].size());
    proxyMap[m_destmesh].m_proxy[i].processCollisions(
        proxyMap[m_sourcemesh].m_proxy,
        proxyMap[m_sourcemesh].m_nchare,
        proxyMap[m_sourcemesh].m_firstchunk,
        static_cast<int>(separated[i].size()),
        separated[i].data() );
  }
}

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

} // exam2m::

#include "NoWarning/controller.def.h"
