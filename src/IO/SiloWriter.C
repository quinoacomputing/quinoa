//******************************************************************************
/*!
  \file      src/IO/SiloWriter.C
  \author    J. Bakosi
  \date      Sun 21 Jul 2013 03:20:47 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Silo (https://wci.llnl.gov/codes/silo) writer
  \details   Silo (https://wci.llnl.gov/codes/silo) writer
*/
//******************************************************************************

#include <iostream>

#include <Exception.h>
#include <SiloWriter.h>
#include <STLMesh.h>

using namespace Quinoa;

SiloWriter::SiloWriter(const std::string& filename,
                       STLMesh* const mesh,
                       const int errLevel) :
  m_filename(filename),
  m_mesh(mesh),
  m_dbfile(nullptr)
//******************************************************************************
//  Constructor
//  \param[in]  filename  Name of Silo file to be created
//  \param[in]  mesh      Mesh object
//  \param[in]  errLevel  Silo library error output level
//! \author J. Bakosi
//******************************************************************************
{
  // Set Silo library error level
  DBShowErrors(errLevel, NULL);

  // Create Silo file
  m_dbfile = DBCreate(filename.c_str(), 0, DB_LOCAL, filename.c_str(), DB_PDB);
  ErrChk(m_dbfile != NULL, ExceptType::FATAL,
        "Cannot create Silo file" + filename);
}

SiloWriter::~SiloWriter() noexcept
//******************************************************************************
//  Destructor: Close Silo file
//! \details    Exception safety: no-throw guarantee: never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  DBClose(m_dbfile);
}

void
SiloWriter::write()
//******************************************************************************
//  Write out Silo file
//! \author J. Bakosi
//******************************************************************************
{
  real* coords[] = { m_mesh->getx(), m_mesh->gety(), m_mesh->getz() };
  int nnodes = m_mesh->nnodes();
  int nfaces = nnodes/3;
  int zshapesize = 3;
  int zshapecnt = nfaces;
  int err;

  err = DBPutFacelist(m_dbfile, "facelist", nfaces, 3, m_mesh->nodelist(),
                      nnodes, 0, 0, &zshapesize, &zshapecnt, 1, NULL, NULL, 0);
  ErrChk(err == 0, ExceptType::FATAL, "Error in Silo::DBPutFacelist()");

  err = DBPutUcdmesh(m_dbfile, m_mesh->name().c_str(), 3, NULL, coords,
                     nnodes, nfaces, NULL, "facelist", DB_DOUBLE, NULL);
  ErrChk(err == 0, ExceptType::FATAL, "Error in Silo::DBPutUcdmesh()");
}
