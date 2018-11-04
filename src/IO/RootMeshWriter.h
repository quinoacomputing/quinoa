// *****************************************************************************
/*!
  \file      src/IO/RootMeshWriter.h
  \copyright 2016-2017, Los Alamos National Security, LLC.
  \brief     Root mesh-based data writer
  \details   Root mesh-based data writer class declaration.
*/
// *****************************************************************************
#ifndef RootMeshWriter_h
#define RootMeshWriter_h

#include <cstddef>
#include <iosfwd>
#include <vector>
#include <array>

#include "NoWarning/TFile.h"
#include "NoWarning/TTree.h"
#include "NoWarning/TGraph2D.h"

#include "Types.h"

namespace tk {

class UnsMesh;

class RootMeshWriter {

  public:
    //! Constructor: create/open Root file
    explicit RootMeshWriter( const std::string filename, int mode );

    //! Destructor
    ~RootMeshWriter() noexcept;

    //! Write Root mesh to file
    void writeMesh( const UnsMesh& mesh );

    //! Write the names of nodal output variables to Root file
    void writeNodeVarNames( const std::vector< std::string >& nv ) const;

    //!  Write node scalar field to Root file
    void writeNodeScalar( uint64_t it,
                          int varid,
                          const std::vector< tk::real >& var ) const;

    //!  Write time stamp to Root file
    void writeTimeStamp( uint64_t it, tk::real time );

  private:
    //! Write Root header
    void writeHeader( const UnsMesh& mesh );

    //! Write nodes coordinates 
    void writeNodes( const UnsMesh& mesh );

    //! Write element conectivity to ROOT file
    void writeElements( const UnsMesh& mesh ) const;

    //! Write element block to ROOT file
    void writeElemBlock( int& elclass,
                         const std::vector< std::size_t >& inpoel ) const;

    //! Variables for ROOT files, tuples
    TFile* m_rfile = nullptr;
    TTree* m_tree_connect = nullptr;
    TTree* m_friendTree = nullptr;
    
    typedef struct mesh_data {
      int coordinates;
      int triangles;
      std::vector<tk::real> mx_root;
      std::vector<tk::real> my_root;
      std::vector<tk::real> mz_root;
      std::vector<std::size_t> tri_connect;
      std::vector<std::size_t> tet_connect;

      public : mesh_data () {
	coordinates = 0;
	triangles = 0;
      }

      public : mesh_data( std::size_t vertices, std::size_t tri_count ) {
	coordinates = static_cast< int >( vertices );
	triangles = static_cast< int >( tri_count );
      }

    } connect_store; 
    
    // declare the object for the connectivity
    connect_store* m_csobject = nullptr;

    const std::string m_filename;          //!< File name
};

} // tk::

#endif // RootMeshWriter_h
