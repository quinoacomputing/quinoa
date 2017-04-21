// *****************************************************************************
/*!
  \file      src/LinSys/HypreMatrix.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Hypre IJ matrix class
  \details   Hypre IJ matrix class.
*/
// *****************************************************************************
#ifndef HypreMatrix_h
#define HypreMatrix_h

#include <HYPRE.h>

#include "Macro.h"
#include "NoWarning/HYPRE_parcsr_ls.h"

namespace tk {
namespace hypre {

//! Hypre matrix class
class HypreMatrix {

  public:
    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wold-style-cast"
    #endif

    //! Create and initialize Hypre IJ matrix
    void create( std::size_t lower, std::size_t upper ) {
      // Create Hypre IJ matrix
      HYPRE_IJMatrixCreate( MPI_COMM_WORLD,
                            static_cast< int >( lower+1 ),
                            static_cast< int >( upper ),
                            static_cast< int >( lower+1 ),
                            static_cast< int >( upper ),
                            &m_A );
      // Choose parallel CSR format storage
      HYPRE_IJMatrixSetObjectType( m_A, HYPRE_PARCSR );
      // Initialize before setting coefficients
      HYPRE_IJMatrixInitialize( m_A );
    }

    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

    //! \brief Destructor: destroy Hypre IJ matrix
    //! \note There is no problem with destroying a Hypre matrix without
    //!   creating it.
    ~HypreMatrix() noexcept { HYPRE_IJMatrixDestroy( m_A ); }

    //! Set values of matrix
    void set( int nrows, int* ncols, const int* rows, const int* cols,
              const double* values )
    { HYPRE_IJMatrixSetValues( m_A, nrows, ncols, rows, cols, values ); }

    //! Assemble matrix
    void assemble() { HYPRE_IJMatrixAssemble( m_A ); }

    //! Print out matrix to file (for debugging)
    //! \param[in] filename Base file name
    //! \details The files names will be <filename>.XXXXX, where XXXXX is the
    //!   processor id.
    void print( const std::string& filename )
    { HYPRE_IJMatrixPrint( m_A, filename.c_str() ); }

    //! Hypre matrix accessor
    HYPRE_ParCSRMatrix get() const {
      HYPRE_ParCSRMatrix m;
      #if defined(STRICT_GNUC)
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wstrict-aliasing"
      #endif
      HYPRE_IJMatrixGetObject( m_A, reinterpret_cast<void**>(&m) );
      #if defined(STRICT_GNUC)
        #pragma GCC diagnostic pop
      #endif
      return m;
    }

  private:
    HYPRE_IJMatrix m_A;         //!< Hypre IJ matrix
};

} // hypre::
} // tk::

#endif // HypreMatrix_h
