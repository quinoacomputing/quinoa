// *****************************************************************************
/*!
  \file      src/LinearSolver/CSR.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Compressed sparse row (CSR) storage for a sparse matrix
  \details   Compressed sparse row (CSR) storage for a sparse matrix.
*/
// *****************************************************************************
#ifndef CSR_h
#define CSR_h

#include <vector>
#include <ostream>

#include "NoWarning/pup_stl.hpp"

#include "Types.hpp"

namespace tk {

//! Compressed sparse row (CSR) storage for a sparse matrix
class CSR {

  public:
    //! \brief Constructor: Create a CSR symmetric matrix with DOF degrees of
    //!   freedom, storing only the upper triangular part.
    explicit CSR( std::size_t DOF,
                  const std::pair< std::vector< std::size_t >,
                                   std::vector< std::size_t > >& psup );

    //! Return const reference to sparse matrix entry at a position
    const tk::real&
    operator()( std::size_t row, std::size_t col, std::size_t pos=0 ) const;

    //! Return non-const reference to sparse matrix entry at a position
    //! \see "Avoid Duplication in const and Non-const Member Function," and
    //!   "Use const whenever possible," Scott Meyers, Effective C++, 3d ed.
    tk::real&
    operator()( std::size_t row, std::size_t col, std::size_t pos=0 ) {
      return const_cast< tk::real& >(
               static_cast< const CSR& >( *this ).operator()( row, col, pos ) );
    }

    //! Access real size of matrix
    std::size_t rsize() const { return rnz.size()*dof; }

    //! Write out CSR as stored
    std::ostream& write_as_stored( std::ostream &os ) const;
    //! Write out CSR nonzero structure
    std::ostream& write_as_structure( std::ostream &os ) const;
    //! Write out CSR as a real matrix
    std::ostream& write_as_matrix( std::ostream &os ) const;
    //! Write out CSR in Matlab/Octave format
    std::ostream& write_as_matlab( std::ostream &os ) const;

    /** @name Pack/unpack (Charm++ serialization) routines */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) {
      p | dof;
      p | rnz;
      p | ia;
      p | ja;
      p | a;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] c CSR object reference
    friend void operator|( PUP::er& p, CSR& c ) { c.pup(p); }
    ///@}

  private:
    std::size_t dof;                    //!< Number of degrees of freedom
    std::vector< std::size_t > rnz;     //!< Number of nonzeros of each row
    std::vector< std::size_t > ia;      //!< Row pointers
    std::vector< std::size_t > ja;      //!< Column indices
    std::vector< tk::real > a;          //!< Nonzero values
};

} // tk::

#endif // CSR_h
