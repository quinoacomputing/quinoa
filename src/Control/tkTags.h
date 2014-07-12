//******************************************************************************
/*!
  \file      src/Control/tkTags.h
  \author    J. Bakosi
  \date      Sun 01 Jun 2014 11:46:02 AM MDT
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

} // tag::
} // tk::

#endif // tkTags_h
