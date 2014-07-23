//******************************************************************************
/*!
  \file      src/Control/tkTags.h
  \author    J. Bakosi
  \date      Wed 23 Jul 2014 07:18:46 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Toolkit Tags
  \details   Toolkit Tags
*/
//******************************************************************************
#ifndef tkTags_h
#define tkTags_h

namespace tk {
namespace tag {

struct seed {};
struct uniform_method {};
struct gaussian_method {};
struct rng {};
struct rngmkl {};
struct rngsse {};
struct seqlen {};
struct verbose {};

} // tag::
} // tk::

#endif // tkTags_h
