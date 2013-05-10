//******************************************************************************
/*!
  \file      src/Model/Mix/Mix.C
  \author    J. Bakosi
  \date      Fri May 10 09:00:47 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mix model base
  \details   Mix model base
*/
//******************************************************************************

// #include <Mix.h>
// #include <Control.h>
// #include <Exception.h>
// 
// using namespace Quinoa;
// 
// Mix::Mix(Memory* const memory,
//          Paradigm* const paradigm,
//          Control* const control,
//          const string& name) :
//   Model(memory, paradigm, control, name, control->get<control::NPAR>()),
//   m_nscalar(control->get<control::NSCALAR>())
// //******************************************************************************
// //  Constructor
// //! \param[in]  memory   Memory object pointer
// //! \param[in]  paradigm Parallel programming object pointer
// //! \param[in]  control  Control object pointer
// //! \param[in]  name     Mix model name
// //! \author  J. Bakosi
// //******************************************************************************
// {
//   ErrChk(m_nscalar > 0, FATAL, "Wrong number of scalars");
// }
