// *****************************************************************************
/*!
  \file      src/LinSys/HypreVector.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Hypre IJ vector class
  \details   Hypre IJ vector class.
*/
// *****************************************************************************
#ifndef HypreVector_h
#define HypreVector_h

#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>

#include "Macro.h"

namespace tk {
namespace hypre {

//! Hypre vector class
class HypreVector {

  public:
    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wold-style-cast"
    #endif

    //! Create and initialize Hypre IJ vector
    void create( std::size_t lower, std::size_t upper ) {
      // Create Hypre IJ vector
      HYPRE_IJVectorCreate( MPI_COMM_WORLD,
                            static_cast< int >( lower+1 ),
                            static_cast< int >( upper ),
                            &m_v );
      // Choose parallel CSR format storage
      HYPRE_IJVectorSetObjectType( m_v, HYPRE_PARCSR );
      // Initialize before setting coefficients
      HYPRE_IJVectorInitialize( m_v );
    }

    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

    //! \brief Destructor: destroy Hypre IJ vector
    //! \note There is no problem with destroying a Hypre vector without
    //!   creating it.
    ~HypreVector() noexcept { HYPRE_IJVectorDestroy( m_v ); }

    //! Set values of vector
    void set( int nvalues, const int* indices, const double* values )
    { HYPRE_IJVectorSetValues( m_v, nvalues, indices, values ); }

    //! Get the local vector
    void get( int nvalues, const int* indices, double* values )
    { HYPRE_IJVectorGetValues( m_v, nvalues, indices, values ); }

    //! Assemble vector
    void assemble() { HYPRE_IJVectorAssemble( m_v ); }

    //! Print out vector to file (for debugging)
    //! \param[in] filename Base file name
    //! \details The files names will be <filename>.XXXXX, where XXXXX is the
    //!   processor id.
    void print( const std::string& filename )
    { HYPRE_IJVectorPrint( m_v, filename.c_str() ); }

    //! Hypre vector accessor
    HYPRE_ParVector get() const {
      HYPRE_ParVector v;
      #if defined(STRICT_GNUC)
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wstrict-aliasing"
      #endif
      HYPRE_IJVectorGetObject( m_v, reinterpret_cast<void**>(&v) );
      #if defined(STRICT_GNUC)
        #pragma GCC diagnostic pop
      #endif
      return v;
    }

  private:
    HYPRE_IJVector m_v;         //!< Hypre IJ vector
};

} // hypre::
} // tk::

#endif // HypreVector_h
