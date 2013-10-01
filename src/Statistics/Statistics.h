//******************************************************************************
/*!
  \file      src/Statistics/Statistics.h
  \author    J. Bakosi
  \date      Mon 30 Sep 2013 10:13:49 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Statistics
  \details   Statistics
*/
//******************************************************************************
#ifndef Statistics_h
#define Statistics_h

#include <QuinoaTypes.h>
#include <Base.h>
#include <Distribution.h>
#include <Quinoa/InputDeck/Types.h>

namespace quinoa {

class Physics;

//! Statistics estimator
class Statistics {

  public:
    //! Constructor
    explicit Statistics(const Base& base, const real* const particles);

    //! Destructor
    virtual ~Statistics() = default;

    //! Accumulate statistics
    void accumulate();

    //! Number of ordinary moments accessor
    int nord() const noexcept { return m_nord; }

    //! Number of central moments accessor
    int ncen() const noexcept { return m_ncen; }

    //! Ordinary moments accessor
    const real* ordinary() const noexcept { return m_ordinary.get(); }

    //! Central moments accessor
    const real* central() const noexcept { return m_central.get(); }

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

    //! Estimate ordinary moments
    void estimateOrdinary();

    //! Estimate central moments
    void estimateCentral();

    //! Find out whether product only contains ordinary moment terms
    bool ordinary(const std::vector<ctr::Term>& product) const;

    //! Return mean for fluctuation
    int mean(const ctr::Term& term) const;

    //! Convert string to upper case
    std::string toUpper(const std::string& s) const;

    //! Return true if string is all lower case
    bool isLower(const std::string&s) const;

    const Base& m_base;                       //!< Essentials
    const uint64_t m_nthread;                 //!< Number of threads
    const uint64_t m_npar;                    //!< Number of particles
    const real* const m_particles;            //!< Particles
    const int m_nprop;                        //!< Number of particle properties
    const std::vector<ctr::Product> m_statistics;//!< Requested tatistics

    //! Instantaneous variable pointers for computing ordinary moments
    std::vector<std::vector<const real*>> m_instOrd;
    std::unique_ptr<real[]> m_ordinary;       //!< Ordinary moments
    std::vector<bool> m_plotOrdinary;         //!< Whether to plot ord moments
    //! Ordinary moment field names
    std::vector<ctr::FieldName> m_ordFieldName;
    std::vector<std::string> m_nameOrdinary;  //!< Ordinary moment names
    int m_nord;                               //!< Number of ordinary moments

    //! Instantaneous variable pointers for computing central moments
    std::vector<std::vector<const real*>> m_instCen;
    std::unique_ptr<real[]> m_central;        //!< Central moments
    //! Ordinary moments about which to compute central moments
    std::vector<std::vector<const real*>> m_center;
    std::vector<std::string> m_nameCentral;   //!< Central moment names
    int m_ncen;                               //!< Number of central moments
};

} // namespace quinoa

#endif // Statistics_h
