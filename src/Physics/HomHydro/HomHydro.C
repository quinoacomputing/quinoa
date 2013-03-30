//******************************************************************************
/*!
  \file      src/Physics/HomHydro/HomHydro.C
  \author    J. Bakosi
  \date      Sat 30 Mar 2013 12:53:26 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous hydrodynamics
  \details   Homogeneous hydrodynamics
*/
//******************************************************************************

#include <sstream>
#include <iomanip>

#include <Memory.h>
#include <MemoryException.h>
#include <ControlTypes.h>
#include <Control.h>
#include <HomHydro.h>
#include <HydroException.h>
#include <PDFWriter.h>
#include <GlobWriter.h>
#include <TxtPlotWriter.h>
#include <SimplifiedLangevin.h>
#include <GeneralizedLangevin.h>
#include <Timer.h>
#include <Statistics.h>

using namespace Quinoa;

HomHydro::HomHydro(Memory* const memory,
                   Paradigm* const paradigm,
                   Control* const control,
                   Timer* const timer) :
  Physics(memory, paradigm, control, timer),
  m_term(control->get<control::TERM>()),
  m_totalTime(timer->create("Total solution"))
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  control  Control object pointer
//! \param[in]  timer    Timer object pointer
//! \author  J. Bakosi
//******************************************************************************
{
  // Instantiate selected hydrodynamics model
  switch (control->get<control::HYDRO>()) {

    case control::HydroType::NO_HYDRO :
      Throw(HydroException,FATAL,HydroExceptType::NO_SUCH_HYDRO);
      break;

    case control::HydroType::SLM :
      m_hydro = new (nothrow) SimplifiedLangevin(memory, paradigm, control);
      Assert(m_hydro != nullptr, MemoryException,FATAL,BAD_ALLOC);
      break;

    case control::HydroType::GLM :
      m_hydro = new (nothrow) GeneralizedLangevin(memory, paradigm, control);
      Assert(m_hydro != nullptr, MemoryException,FATAL,BAD_ALLOC);
      break;

    default :
      Throw(HydroException,FATAL,HYDRO_UNIMPLEMENTED);
  }

  // Instantiate statistics estimator
  m_statistics = new (nothrow) Statistics(memory, paradigm, control, m_hydro);
  Assert(m_statistics != nullptr, MemoryException,FATAL,BAD_ALLOC);

  // Instantiate glob file writer
  m_glob = new (nothrow) GlobWriter(m_control->get<control::GLOBNAME>());
  Assert(m_glob != nullptr, MemoryException,FATAL,BAD_ALLOC);

  // Instantiate plot file writer
  m_plot = new (nothrow) TxtPlotWriter(m_control->get<control::PLOTNAME>(),
                                       m_statistics);
  Assert(m_plot != nullptr, MemoryException,FATAL,BAD_ALLOC);
}

HomHydro::~HomHydro()
//******************************************************************************
//  Destructor
//! \author  J. Bakosi
//******************************************************************************
{
  if (m_plot) { delete m_plot; m_plot = nullptr; }
  if (m_glob) { delete m_glob; m_glob = nullptr; }
  if (m_statistics) { delete m_statistics; m_statistics = nullptr; }
  if (m_hydro) { delete m_hydro; m_hydro = nullptr; }
}

void
HomHydro::solve()
//******************************************************************************
//  Solve
//! \author  J. Bakosi
//******************************************************************************
{
  int it = 0;
  real t = 0.0;
  bool wroteJpdf = false;
  bool wroteGlob = false;
  bool wrotePlot = false;

  const auto nstep = m_control->get<control::NSTEP>();
  const auto ttyi  = m_control->get<control::TTYI>();
  const auto pdfi  = m_control->get<control::PDFI>();
  const auto glob  = m_control->get<control::GLOB>();
  const auto plti  = m_control->get<control::PLTI>();
  const auto dt    = m_control->get<control::DT>();

  m_timer->start(m_totalTime);

  // Echo headers
  if (nstep) {
    reportHeader();
    m_plot->header();
  }

  // Time stepping loop
  while (fabs(t-m_term) > numeric_limits<real>::epsilon() && it < nstep) {

    // Advance particles
    m_hydro->advance(dt);

    // Accumulate statistics
    m_statistics->accumulate();

    // Output pdf at selected times
    if (!(it % pdfi)) { outJpdf(t); wroteJpdf = true; }

    // Append glob file at selected times
    if (!(it % glob)) { m_glob->write(it,t); wroteGlob = true; }

    // Append plot file at selected times
    if (!(it % plti)) { m_plot->write(it,t); wrotePlot = true; }

    // Echo one-liner info
    if (!(it % ttyi)) {
      report(it, nstep, t, dt, wroteJpdf, wroteGlob, wrotePlot);
      wroteJpdf = wroteGlob = wrotePlot = false;
    }

    // Increase timestep and iteration counter
    t += dt;
    ++it;
    if (t > m_term) t = m_term;
  } // Time stepping loop
}

void
HomHydro::reportHeader()
//******************************************************************************
//  Echo report header
//! \author  J. Bakosi
//******************************************************************************
{
  cout << "      it             t            dt"
          "        ETE        ETA   out\n"
       << "------------------------------------"
          "----------------------------" << endl;
}

void
HomHydro::report(const int it,
                 const int nstep,
                 const real t,
                 const real dt,
                 const bool wroteJpdf,
                 const bool wroteGlob,
                 const bool wrotePlot)
//******************************************************************************
//  One-liner report
//! \param[in]  it         Iteration counter
//! \param[in]  nstep      Terminate time
//! \param[in]  t          Time
//! \param[in]  dt         Time step size
//! \param[in]  wroteJpdf  True if joint PDF was output
//! \param[in]  wroteGlob  True if glob was output
//! \param[in]  wrotePlot  True if plot was output
//! \author  J. Bakosi
//******************************************************************************
{
  Watch ete, eta;       // estimated time elapsed and to accomplishment
  m_timer->eta(m_totalTime, m_term, t, nstep, it, ete, eta);

  cout << setfill(' ') << setw(8) << it << "  "
       << scientific << setprecision(6) << setw(12) << t << "  " << dt << "  "
       << setfill('0') << setw(3) << ete.h.count() << ":"
                       << setw(2) << ete.m.count() << ":"
                       << setw(2) << ete.s.count() << "  "
                       << setw(3) << eta.h.count() << ":"
                       << setw(2) << eta.m.count() << ":"
                       << setw(2) << eta.s.count() << "  ";

  if (wroteGlob) cout << "G";
  if (wroteJpdf) cout << "J";
  if (wrotePlot) cout << "P";

  cout << endl;
}

void
HomHydro::outJpdf(const real t)
//******************************************************************************
//  Output joint PDF
//! \param[in]  t    Time stamp
//! \author  J. Bakosi
//******************************************************************************
{
//   // Contruct filename
//   stringstream ss;
//   ss << m_control->get<control::PDFNAME>() << "." << t << ".msh";
//   string filename = ss.str();
// 
//   // Create joint PDF
//   JPDF jpdf(m_nscalar, 0.02);
// 
//   // Estimate joint PDF
//   m_mix->jpdf(jpdf);
// 
//   // Output joint PDF
//   PDFWriter jpdfFile(filename);
//   jpdfFile.writeGmsh(&jpdf);
}

void
HomHydro::echo()
//******************************************************************************
//  Echo informaion on homogeneous hydrodynamics
//! \author  J. Bakosi
//******************************************************************************
{
}

void
HomHydro::init()
//******************************************************************************
//  Initialize homogeneous hydrodynamics
//! \author  J. Bakosi
//******************************************************************************
{
  m_hydro->init();
}
