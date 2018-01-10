// *****************************************************************************
/*!
  \file      src/Inciter/AMR/Error.h
  \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
  \brief     Class for computing error estimates for mesh refinement
  \details   Class for computing error estimates for mesh refinement.
*/
// *****************************************************************************
#ifndef Error_h
#define Error_h

#include "PUPUtil.h"
#include "Fields.h"

namespace AMR {

//! Class for computing error estimates for mesh refinement
class Error {

  public:
    //! Constructor
    explicit Error();

    //! Compute error estimate for a scalar quantity
    void scalar( const tk::Fields& u );

    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) {
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] d Error object reference
    friend void operator|( PUP::er& p, Error& d ) { d.pup(p); }
    //@}
};

} // AMR::

#endif // Error_h
