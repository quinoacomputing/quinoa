//******************************************************************************
/*!
  \file      src/LinSys/HypreVector.h
  \author    J. Bakosi
  \date      Tue 03 May 2016 09:19:53 AM MDT
  \copyright 2012-2016, Jozsef Bakosi.
  \brief     Hypre IJ vector class
  \details   Hypre IJ vector class.
*/
//******************************************************************************
#ifndef HypreVector_h
#define HypreVector_h

#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>

namespace tk {
namespace hypre {

//! Hypre vector class
class HypreVector {

  public:
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
      HYPRE_IJVectorGetObject( m_v, reinterpret_cast<void**>(&v) );
      return v;
    }

  private:
    HYPRE_IJVector m_v;         //!< Hypre IJ vector
};

} // hypre::
} // tk::

#endif // HypreVector_h
