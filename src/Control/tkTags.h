//******************************************************************************
/*!
  \file      src/Control/tkTags.h
  \author    J. Bakosi
  \date      Thu 16 Jan 2014 09:33:15 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
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
struct mklrng {};
struct rngsse {};
struct seqlen {};

} // tag::
} // tk::

#endif // tkTags_h
