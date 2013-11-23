//******************************************************************************
/*!
  \file      src/RNG/MatrixRank.h
  \author    J. Bakosi
  \date      Fri 22 Nov 2013 05:30:31 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Statistical tests suggested by George Marsaglia
  \details   Statistical tests suggested by George Marsaglia
*/
//******************************************************************************
#ifndef MatrixRank_h
#define MatrixRank_h

#include <Marsaglia.h>

namespace rngtest {

//! MatrixRank : Marsaglia
class MatrixRank : public Marsaglia {

  public:
    //! Constructor
    explicit MatrixRank() = default;

    //! Destructor
    ~MatrixRank() noexcept override = default;

  private:
    //! Don't permit copy constructor
    MatrixRank(const MatrixRank&) = delete;
    //! Don't permit copy assigment
    MatrixRank& operator=(const MatrixRank&) = delete;
    //! Don't permit move constructor
    MatrixRank(MatrixRank&&) = delete;
    //! Don't permit move assigment
    MatrixRank& operator=(MatrixRank&&) = delete;
};

} // rngtest::

#endif // MatrixRank_h
