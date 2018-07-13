// *****************************************************************************
/*!
  \file      src/Base/CartesianProduct.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Cartesian product using brigand
  \details   Cartesian product using brigand
*/
// *****************************************************************************
#ifndef CartesianProduct_h
#define CartesianProduct_h

#include <brigand/algorithms/fold.hpp>
#include <brigand/algorithms/transform.hpp>

namespace tk {

//! Cartesian product of two brigand lists
//! \see <brigand>/test/apply.cpp
template< class li, class lo >
using cartesian_product = brigand::reverse_fold<brigand::list<li, lo>, brigand::list<brigand::list<>>, brigand::bind<brigand::join, brigand::bind<brigand::transform, brigand::_2, brigand::defer<brigand::bind<brigand::join, brigand::bind<brigand::transform, brigand::parent<brigand::_1>, brigand::defer<brigand::bind<brigand::list, brigand::bind<brigand::push_front, brigand::_1, brigand::parent<brigand::_1>>>>>>>>>>;

} // tk::

#endif // CartesianProduct_h
