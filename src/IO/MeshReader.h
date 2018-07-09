// *****************************************************************************
/*!
  \file      src/IO/MeshReader.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Polymorphic mesh reader class for connecting to various readers
  \brief     Polymorphic mesh reader class for connecting to various lower
    level, specific mesh readers.
*/
// *****************************************************************************
#ifndef MeshReader_h
#define MeshReader_h

#include <vector>
#include <array>
#include <string>

#include "Types.h"
#include "MeshFactory.h"
#include "Make_unique.h"
#include "ExodusIIMeshReader.h"

#ifdef HAS_OMEGA_H
  #include "Omega_h_MeshReader.h"
#endif

namespace tk {

//! Polymorphic mesh reader class for connecting to various mesh readers
//! \details This class uses runtime polymorphism without client-side
//!   inheritance: inheritance is confined to the internals of the this class,
//!   invisible to client-code. The class exclusively deals with ownership
//!   enabling client-side value semantics. Credit goes to Sean Parent at Adobe.
//! \see http://sean-parent.stlab.cc/papers-and-presentations/#value-semantics-and-concept-based-polymorphism.
//! \see For example client code that models a MeshReader, see
//!   tk::ExodusIIMeshReader or tk::Omega_h_MeshReader.
class MeshReader {

  public:
    //! Constructor
    //! \param[in] filename Input mesh filename
    //! \details Dispatch constructor call to various low level mesh readers by
    //!    creating child class and assigning to base to be used in polymorphic
    //!    fashion.
    explicit MeshReader( const std::string& filename ) {
      auto meshtype = detectInput( filename );
      if (meshtype == MeshReaderType::EXODUSII) {
        using R = ExodusIIMeshReader;
        self = make_unique< Model<R> >( R(filename) );
      #ifdef HAS_OMEGA_H
      } else if (meshtype == MeshReaderType::OMEGA_H) {
        using R = Omega_h_MeshReader;
        self = make_unique< Model<R> >( R(filename) );
      #endif
      } else Throw( "Mesh type not implemented or not supported" );
    }

    //! Public interface to read part of the mesh (graph and coords) from file
    //! \details Total number of PEs defaults to 1 for a single-CPU read, this
    //!    PE defaults to 0 for a single-CPU read.
    void readMeshPart( std::vector< std::size_t >& ginpoel,
                       std::vector< std::size_t >& inpoel,
                       std::vector< std::size_t >& gid,
                       std::unordered_map< std::size_t, std::size_t >& lid,
                       tk::UnsMesh::Coords& coord, 
                       int numpes=1, int mype=0 )
    { self->readMeshPart( ginpoel, inpoel, gid, lid, coord, numpes, mype ); }

    //! Public interface to read face list of side sets from mesh file
    std::size_t
    readSidesetFaces( std::map< int, std::vector< std::size_t > >& belem,
                      std::map< int, std::vector< int > >& faceid )
    { return self->readSidesetFaces( belem, faceid ); }

    //! Public interface to read face connectivity of boundary faces from file
    void readFaces( std::size_t nbfac, std::vector< std::size_t >& conn )
    { self->readFaces( nbfac, conn ); }

    //! Public interfaces to read node list of all side sets from mesh file
    std::map< int, std::vector< std::size_t > > readSidesets()
    { return self->readSidesets(); }

    //! Copy assignment
    MeshReader& operator=( const MeshReader& x )
    { MeshReader tmp(x); *this = std::move(tmp); return *this; }
    //! Copy constructor
    MeshReader( const MeshReader& x ) : self( x.self->copy() ) {}
    //! Move assignment
    MeshReader& operator=( MeshReader&& ) noexcept = default;
    //! Move constructor
    MeshReader( MeshReader&& ) noexcept = default;

  private:
    //! \brief Concept is a pure virtual base class specifying the requirements
    //!   of polymorphic objects deriving from it
    struct Concept {
      Concept() = default;
      Concept( const Concept& ) = default;
      virtual ~Concept() = default;
      virtual Concept* copy() const = 0;
      virtual void readMeshPart(
                     std::vector< std::size_t >&,
                     std::vector< std::size_t >&,
                     std::vector< std::size_t >&,
                     std::unordered_map< std::size_t, std::size_t >&,
                     tk::UnsMesh::Coords&,
                     int, int ) = 0;
      virtual std::size_t
        readSidesetFaces( std::map< int, std::vector< std::size_t > >&,
                          std::map< int, std::vector< int > >& ) = 0;
      virtual void readFaces( std::size_t, std::vector< std::size_t >& )
        const = 0;
      virtual std::map< int, std::vector< std::size_t > > readSidesets() = 0;
    };

    //! \brief Model models the Concept above by deriving from it and overriding
    //!   the the virtual functions required by Concept
    template< typename T >
    struct Model : Concept {
      Model( T x ) : data( std::move(x) ) {}
      Concept* copy() const override { return new Model( *this ); }
      void readMeshPart( std::vector< std::size_t >& ginpoel,
                         std::vector< std::size_t >& inpoel,
                         std::vector< std::size_t >& gid,
                         std::unordered_map< std::size_t, std::size_t >& lid,
                         tk::UnsMesh::Coords& coord, 
                         int numpes, int mype ) override
        { data.readMeshPart( ginpoel, inpoel, gid, lid, coord, numpes, mype ); }
      std::size_t
        readSidesetFaces( std::map< int, std::vector< std::size_t > >& belem,
                          std::map< int, std::vector< int > >& faceid )
        override { return data.readSidesetFaces( belem, faceid ); }
      void readFaces( std::size_t nbfac, std::vector< std::size_t >& conn )
        const override { data.readFaces( nbfac, conn ); }
      std::map< int, std::vector< std::size_t > > readSidesets() override
        { return data.readSidesets(); }
      T data;
    };

    std::unique_ptr< Concept > self;    //!< Base pointer used polymorphically
};

} // tk::

#endif // MeshReader_h
