// *****************************************************************************
/*!
  \file      src/IO/SiloWriter.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Silo writer definition.
  \details   Silo writer definition. This class currently only supports writing
    an STL triangulation into a Silo file. See also
    https://wci.llnl.gov/codes/silo.
*/
// *****************************************************************************

#include <string>
#include <sstream>
#include <algorithm>
#include <cstring>

#include "Exception.h"
#include "SiloWriter.h"
#include "STLMesh.h"

void
tk::SiloError( char* msg )
// *****************************************************************************
//  Silo error handler
//! \param[in]  msg  Error message
//! \author J. Bakosi
// *****************************************************************************
{
  // Take out newlines from error message coming from library
  std::string str(msg);
  std::replace(str.begin(), str.end(), '\n', ' ');

  // Echo and throw
  std::stringstream ss;
  ss << "Silo library writer error: " << str;
  Throw( ss.str() );
}

using tk::SiloWriter;

SiloWriter::SiloWriter( const std::string& filename,
                        const STLMesh& mesh,
                        int errLevel ) :
  m_filename( filename ),
  m_mesh( mesh ),
  m_dbfile( nullptr )
// *****************************************************************************
//  Constructor
//  \param[in] filename Name of Silo file to be created
//  \param[in] mesh Mesh object
//  \param[in] errLevel Silo library error output level
//! \author J. Bakosi
// *****************************************************************************
{
  // Save Silo's error handler and reporting level
  m_errFunc = DBErrfunc();
  m_errLevel = DBErrlvl();

  // Set Silo library error level and handler
  DBShowErrors(errLevel, &SiloError);

  // Create Silo file
  m_dbfile = DBCreate(filename.c_str(), 0, DB_LOCAL, filename.c_str(), DB_HDF5);
  ErrChk( m_dbfile != NULL, "Cannot create Silo file" + filename );
}

SiloWriter::~SiloWriter() noexcept
// *****************************************************************************
//  Destructor
//! \details Exception safety: no-throw guarantee: never throws exceptions.
//! \author J. Bakosi
// *****************************************************************************
{
  // Close Silo file
  DBClose(m_dbfile);

  // Restore Silo library error level and handler
  DBShowErrors(m_errLevel, m_errFunc);
}

void
SiloWriter::write()
// *****************************************************************************
//  Write out Silo file
//! \author J. Bakosi
// *****************************************************************************
{
  tk::real* coords[] = { m_mesh.getx(), m_mesh.gety(), m_mesh.getz() };
  int nnodes = m_mesh.nnodes();
  int nfaces = nnodes/3;
  int zshapetype = DB_ZONETYPE_TRIANGLE;
  int zshapesize = 3;
  int zshapecnt = nfaces;

  // Write out STL face connectivity
  DBPutZonelist2(m_dbfile, "zonelist", nfaces, 3, m_mesh.nodelist(), nnodes,
                 0, 0, 0, &zshapetype, &zshapesize, &zshapecnt, 1, NULL);

  // Write out STL mesh: no zones, only faces with a simple face connectivity
  DBPutUcdmesh(m_dbfile, m_mesh.name().c_str(), 3, NULL, coords,
               nnodes, nfaces, "zonelist", NULL, DB_DOUBLE, NULL);
}
