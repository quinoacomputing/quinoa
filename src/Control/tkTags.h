//******************************************************************************
/*!
  \file      src/Control/tkTags.h
  \author    J. Bakosi
  \date      Fri 22 Aug 2014 10:35:06 PM MDT
  \copyright 2012-2014, Jozsef Bakosi.
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
struct error {};

} // tag::
} // tk::

#endif // tkTags_h
