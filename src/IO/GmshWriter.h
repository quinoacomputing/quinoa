//******************************************************************************
/*!
  \file      src/IO/GmshWriter.h
  \author    J. Bakosi
  \date      Thu 13 Sep 2012 04:12:16 AM KST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     GmshWriter base class declaration
  \details   GmshWriter base class declaration
*/
//******************************************************************************
#ifndef GmshWriter_h
#define GmshWriter_h

#include <string>

using namespace std;

#include <PlotWriter.h>

namespace Quinoa {

//! GmshWriter base class
class GmshWriter : PlotWriter {

  public:
    //! Constructor
    GmshWriter(string filename, UnsMesh* mesh, Memory* memory) :
      PlotWriter(filename, mesh, memory) {}

    //! Destructor
    ~GmshWriter();

  private:
    //! Don't permit copy operator
    GmshWriter(const GmshWriter&);
    //! Don't permit assigment operator
    GmshWriter& operator=(const GmshWriter&);
};

} // namespace Quinoa

#endif // GmshWriter_h
