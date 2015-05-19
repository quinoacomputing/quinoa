//******************************************************************************
/*!
  \file      src/Base/LinearMap.h
  \author    J. Bakosi
  \date      Tue 19 May 2015 02:46:25 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Advanced Charm++ array creation with a map in a linear fashion
  \details   Advanced Charm++ array creation refers to various ways arrays can
     be created with the Charm++ runtime system. See
     http://charm.cs.illinois.edu/manuals/html/charm++/manual.html, Sec.
     Advanced Array Creation. This class does a simple linear distribution.
*/
//******************************************************************************
#ifndef BlockMap_h
#define BlockMap_h

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <linearmap.decl.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#include <Exception.h>

namespace tk {

//! Charm++ array map for initial placement of array elements in linear fashion
//! \details The map object is used by the Charm++ array manager to determine
//!   the "home" PE of each element. The home PE is the PE upon which the array
//!   element is initially placed, which will retain responsibility for
//!   maintaining the location of the element.
class LinearMap : public CkArrayMap {

  public:
    //! Constructor
    //! \param[in] nelem Total number of array elements
    explicit LinearMap( int nelem ) :
      m_chunksize( nelem > CkNumPes() ? nelem/CkNumPes() : 1 )
    { Assert( nelem > 0, "Number of array elements must be positive" ); }

    //! \brief Return the home processor number for the array element for linear
    //!   distribution
    int procNum( int, const CkArrayIndex& idx ) override;

    //! Create initial set of array elements based on linear distribution
    void populateInitial( int, CkArrayIndex& idx, void *msg, CkArrMgr *mgr )
    override;

  private:
    int m_chunksize;            //!< Number of array elements per PE
};

} // tk::

#endif // LinearMap_h
