//******************************************************************************
/*!
  \file      src/Statistics/Statistics.h
  \author    J. Bakosi
  \date      Wed 28 Aug 2013 08:53:53 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Statistics
  \details   Statistics
*/
//******************************************************************************
#ifndef Statistics_h
#define Statistics_h

#include <QuinoaTypes.h>
#include <Memory.h>
#include <Distribution.h>
#include <QuinoaControlTypes.h>

namespace Quinoa {

class Memory;
class Paradigm;
class QuinoaControl;
class Physics;

//! Statistics estimator
class Statistics {

  public:
    //! Constructor
    explicit Statistics(Memory* const memory,
                        Paradigm* const paradigm,
                        QuinoaControl* const control,
                        Physics* const physics);

    //! Destructor
    virtual ~Statistics() noexcept;

    //! Accumulate statistics
    void accumulate();

    //! Number of ordinary moments accessor
    int nord() const noexcept { return m_nord; }

    //! Number of central moments accessor
    int ncen() const noexcept { return m_ncen; }

    //! Ordinary moments accessor
    const real* ordinary() const noexcept { return m_ordinary.ptr; }

    //! Central moments accessor
    const real* central() const noexcept { return m_central.ptr; }

    //! Find out whether ordinary moment is to be plotted
    bool plotOrdinary(const int m) const;

    //! Return the name of ordinary moment
    const std::string& nameOrdinary(const int m) const;

    //! Return the name of central moment
    const std::string& nameCentral(const int m) const;

  private:
    //! Don't permit copy constructor
    Statistics(const Statistics&) = delete;
    //! Don't permit copy assigment
    Statistics& operator=(const Statistics&) = delete;
    //! Don't permit move constructor
    Statistics(Statistics&&) = delete;
    //! Don't permit move assigment
    Statistics& operator=(Statistics&&) = delete;

    //! Finalize, single exit point, called implicitly from destructor or
    //! explicitly from anywhere else
    void finalize() noexcept;

    //! Estimate ordinary moments
    void estimateOrdinary();

    //! Estimate central moments
    void estimateCentral();

    //! Find out whether product only contains ordinary moment terms
    bool ordinary(const std::vector<control::Term>& product) const;

    //! Return mean for fluctuation
    int mean(const control::Term& term) const;

    //! Convert string to upper case
    std::string toUpper(const std::string& s) const;

    //! Return true if string is all lower case
    bool isLower(const std::string&s) const;

    Memory* const m_memory;                   //!< Memory object
    const uint64_t m_nthread;                 //!< Number of threads
    const uint64_t m_npar;                    //!< Number of particles
    Physics* const m_physics;                 //!< Physics object
    const int m_nprop;                        //!< Number of particle properties
    const std::vector<control::Product> m_statistics;//!< Requested tatistics

    //! Instantaneous variable pointers for computing ordinary moments
    std::vector<std::vector<const real*>> m_instOrd;
    Data<real> m_ordinary;                    //!< Ordinary moments
    std::vector<bool> m_plotOrdinary;         //!< Whether to plot ord moments
    //! Ordinary moment field names
    std::vector<control::FieldName> m_ordFieldName;
    std::vector<std::string> m_nameOrdinary;  //!< Ordinary moment names
    int m_nord;                               //!< Number of ordinary moments

    //! Instantaneous variable pointers for computing central moments
    std::vector<std::vector<const real*>> m_instCen;
    Data<real> m_central;                     //!< Central moments
    //! Ordinary moments about which to compute central moments
    std::vector<std::vector<const real*>> m_center;
    std::vector<std::string> m_nameCentral;   //!< Central moment names
    int m_ncen;                               //!< Number of central moments
};

} // namespace Quinoa

#endif // Statistics_h
