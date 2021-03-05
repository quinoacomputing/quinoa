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

#include "NoWarning/charm++.hpp"

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

  //! \brief Operator << for writing a Transfer object to output streams
  //! \param[in,out] os Output stream to write to
  //! \param[in] t Transfer object to write
  //! \return Updated output stream
  friend std::ostream& operator<<( std::ostream& os, const Transfer& t ) {
    os << t.src << '>' << t.dst;
    return os;
  }

  //! \brief Equal operator for, e.g., finding unique elements, used by, e.g.,
  //!    std::unique().
  //! \details Test on src and dst only.
  //! \param[in] transfer Transfer object to compare
  //! \return Boolean indicating if term equals 'this'
  bool operator== ( const Transfer& transfer) const {
    if (src == transfer.src && dst == transfer.dst)
      return true;
    else
      return false;
  }

  //! Less-than operator for ordering, used by, e.g., std::sort().
  //! \details Test on src and dst only.
  //! \param[in] transfer Transfer object to compare
  //! \return Boolean indicating if term is less than 'this'
  bool operator< ( const Transfer& transfer ) const {
    if (src < transfer.src)
      return true;
    else if (src == transfer.src && dst < transfer.dst)
      return true;
    else
      return false;
  }
};

} // inciter::

#endif // Transfer_h
