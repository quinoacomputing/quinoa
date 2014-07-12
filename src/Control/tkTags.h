//******************************************************************************
/*!
  \file      src/Control/tkTags.h
  \author    J. Bakosi
  \date      Sat 12 Jul 2014 09:02:13 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     tk Tags
  \details   tk Tags
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
