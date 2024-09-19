// Controller for the library
#ifndef M2MTransfer_hpp
#define M2MTransfer_hpp

#include "NoWarning/m2mtransfer.decl.h"

#include "collidecharm.h"
#include "Fields.hpp"

#include <iostream>

namespace exam2m {

//! External user interface functions to M2MTransfer
void collisionHandler( [[maybe_unused]] void *param,
                        int nColl,
                        Collision *colls );
void addMesh(CkArrayID p, int elem, CkCallback cb);
void setSourceTets(CkArrayID p, int index, std::vector< std::size_t >* inpoel, tk::UnsMesh::Coords* coords, const tk::Fields& u);
void setDestPoints(CkArrayID p, int index, tk::UnsMesh::Coords* coords, tk::Fields& u, CkCallback cb);

//! LibMain mainchare that creates collidecharm-proxies at startup
class LibMain : public CBase_LibMain {
public:
  LibMain(CkArgMsg* msg);
  explicit LibMain(CkMigrateMessage* msg) : CBase_LibMain(msg) {
    std::cout << "LibMain() migrate ctor cmplt." << std::endl;
  }
  void pup(PUP::er&) {}
  friend void operator|( PUP::er& p, LibMain& m ) { m.pup(p); }
};

//! MeshData class that contains the mesh
class MeshData {
  public:
    CProxy_TransferDetails m_proxy;
    int m_firstchunk;
    int m_nchare;
    void pup(PUP::er& p) {
      p | m_proxy;
      p | m_firstchunk;
      p | m_nchare;
    }
};

//! M2MTransfer chare-group which is inciter's interface to collidecharm
class M2MTransfer : public CBase_M2MTransfer {
  private:
    std::unordered_map<CmiUInt8, MeshData> proxyMap;
    int current_chunk;
    CmiUInt8 m_sourcemesh, m_destmesh;

  public:

    //! Constructor
    M2MTransfer();

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! Migrate constructor
    explicit M2MTransfer( CkMigrateMessage* m ) : CBase_M2MTransfer( m ) {}
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

    //! Register mesh with the mesh-to-mesh transfer library
    void addMesh(CkArrayID p, int elem, CkCallback cb);

    void setMesh(CkArrayID p, MeshData d);
    void setSourceTets(CkArrayID p, int index, std::vector< std::size_t >* inpoel,
                       tk::UnsMesh::Coords* coords, const tk::Fields& u);
    void setDestPoints(CkArrayID p, int index, tk::UnsMesh::Coords* coords,
                       tk::Fields& u, CkCallback cb);
    void distributeCollisions(int nColl, Collision* colls);

    //! Pack/Unpack serialize member function for Charm
    void pup(PUP::er&) {}
    //! Pack/Unpack serialize operator| for Charm
    friend void operator|( PUP::er& p, M2MTransfer& m ) { m.pup(p); }
};

}

#endif // M2MTransfer_hpp
