// *****************************************************************************
/*!
  \file      src/PDE/EoS/MatBlock.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Material block class
  \details   This file defines a material block vector that stores EOS parameters
    for all defined materials.
*/
// *****************************************************************************
#ifndef MatBlock_h
#define MatBlock_h

#include <cmath>
#include "Data.hpp"
#include "EoS/StiffenedGas.hpp"

namespace inciter {

using ncomp_t = kw::ncomp::info::expect::type;

  std::vector< EoS_Base* > m_mat_blk; // EOS material block

  void initialize_material_blk(tk::real g, tk::real ps, int k)
    {

     m_mat_blk.push_back(new StiffenedGas(g, ps, k));

    }
  
    
} //inciter::

#endif // MatBlock_h
