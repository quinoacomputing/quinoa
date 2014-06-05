//******************************************************************************
/*!
  \file      src/Base/TestU01PUP.h
  \author    J. Bakosi
  \date      Wed 21 May 2014 10:18:25 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     TestU01 Charm++ Pack/UnPack utilities
  \details   TestU01 Charm++ Pack/UnPack utilities
*/
//******************************************************************************
#ifndef TestU01PUP_h
#define TestU01PUP_h

inline void operator|( PUP::er& p, sres_Poisson ) {
  std::cout << "pup struct sres_Poisson\n";
}

inline void operator|( PUP::er& p, sknuth_Res2 ) {
  std::cout << "pup struct sknuth_Res2\n";
}

#endif // PUPUtil_h
