// *****************************************************************************
/*!
  \file      src/Inciter/Transfer.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Types for solution transfer between solvers holding different
             meshes
  \details   Types for solution transfer between solvers holding different
             meshes.
*/
// *****************************************************************************
#ifndef Transfer_h
#define Transfer_h

namespace inciter {

//! \brief Description of solution transfer between two solvers holding
//!   different meshes
struct Transfer {
  std::size_t src;              //!< Source mesh id
  std::size_t dst;              //!< Destination mesh id
  std::vector< CkCallback > cb; //!< List of callbacks to continue with
  //std::vector< std::size_t > cbs;     // only for debugging

  //! Constructor: empty for Charm++
  explicit Transfer() = default;

  //! Constructor: initialize src and dest mesh ids
  explicit Transfer( std::size_t s, std::size_t d ) : src(s), dst(d) {}

  /** @name Pack/Unpack: Serialize Transfer object for Charm++ */
  ///@{
  //! Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  void pup( PUP::er& p ) {
    p | src;
    p | dst;
    p | cb;
  }
  //! Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] t Transfer object reference
  friend void operator|( PUP::er& p, Transfer& t ) { t.pup(p); }
  ///@}
};

} // inciter::

#endif // Transfer_h
