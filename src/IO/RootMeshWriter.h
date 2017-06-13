// *****************************************************************************
/*!
  \file      src/IO/RootMeshWriter.h
  \author    A. Pakki 
  \copyright 2012-2015, Aditya Pakki, 2016, Los Alamos National Security, LLC.
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
#include <TGraph2D.h>
#include <TTree.h>

#include "Types.h"

namespace tk {

class UnsMesh;

class RootMeshWriter {

  public:
    //! Constructor: create/open Root file
    explicit RootMeshWriter( const std::string filename, int option);

    //! Destructor
    ~RootMeshWriter() noexcept;

    //! Write Root mesh to file
    void writeMesh( const UnsMesh& mesh ) const;

    //! Write the names of nodal output variables to Root file
    void writeNodeVarNames( const std::vector< std::string >& nv ) const;

    //!  Write node scalar field to Root file
    void writeNodeScalar( uint64_t it,
                          int varid,
                          const std::vector< tk::real >& var ) const;

    //!  Write time stamp to Root file
    void writeTimeStamp( uint64_t it, tk::real time ) const;

  private:
    //! Write Root header
    void writeHeader( const UnsMesh& mesh ) const;

    //! Write nodes coordinates 
    void writeNodes( const UnsMesh& mesh ) const;

    //! Write element conectivity to ROOT file
    void writeElements( const UnsMesh& mesh ) const;

    //! Write element block to ROOT file
    void writeElemBlock( int& elclass,
                         const std::vector< std::size_t >& inpoel ) const;

    //! Variables for ROOT files, tuples
    TFile *rfile = nullptr;
    TTree *tree_connect = nullptr;
    
    typedef struct mesh_data {
      int coordinates;
      int triangles;
      std::vector<float> mx_root;
      std::vector<float> my_root;
      std::vector<float> mz_root;
      std::vector<int> connectivity;

      public : mesh_data () {
	coordinates = 0;
	triangles = 0;
      }

      public : mesh_data( int vertices, int tri_count ) {
	coordinates = vertices;
	triangles = tri_count;
      }

    } connect_store; 
    
    // declare the object for the connectivity
    mutable connect_store *csobject = nullptr;

    const std::string m_filename;          //!< File name
    int m_outFile;                         //!< Root file handle
};

} // tk::

#endif // RootMeshWriter_h
