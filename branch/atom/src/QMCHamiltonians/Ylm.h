//////////////////////////////////////////////////////////////////
// (c) Copyright 2007-  by Ken Esler
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_YLM_H
#define QMCPLUSPLUS_YLM_H

#include <cmath>
#include <complex>
#include "OhmmsPETE/TinyVector.h"

namespace qmcplusplus
{

template<typename T>
inline T LegendrePll (int l, T x)
{
  if (l==0)
    return 1.0;
  else
  {
    T sqt = std::sqrt(1.0-x)*std::sqrt(1.0+x);
    T val = 1.0;
    T dblfact = 1.0;
    for (int i=1; i<=l; i++)
    {
      val *= -sqt;
      val *= dblfact;
      dblfact += 2.0;
    }
    return val;
  }
}

template<typename T>
inline T LegendrePlm (int l, int m, T x)
{
  if (m < 0)
  {
    m = std::abs (m);
    T posval = LegendrePlm (l, m, x);
    T sign = (m%2==0) ? 1.0 : -1.0;
    T mfact = 1.0;
    for (int i=2; i<=(l-m); i++)
      mfact *= static_cast<T>(i);
    T pfact = 1.0;
    for (int i=2; i<=(l+m); i++)
      pfact *= static_cast<T>(i);
    return posval * sign*mfact/pfact;
  }
  // Now we can assume that m is not negative
  T pmm = LegendrePll (m, x);
  T pmp1m = x*(2*m+1)*pmm;
  if (m == l)
    return pmm;
  else
    if (l==(m+1))
      return pmp1m;
    else
      // Use recursive formula
    {
      T Plm2m = pmm;
      T Plm1m = pmp1m;
      T Pl;
      for (int i=m+2; i<=l;  i++)
      {
        Pl = (1.0/static_cast<T>(i-m)) * (x*(2*i-1)*Plm1m - (i+m-1)*Plm2m);
        Plm2m = Plm1m;
        Plm1m = Pl;
      }
      return Pl;
    }
}

template<typename T>
inline std::complex<T> Ylm(int l, int m, const TinyVector<T,3>& r)
{
  T costheta = r[0];
  T phi   = std::atan2(r[2],r[1]);
  int lmm = l - m;
  int lpm = l + m;
  T mfact = 1.0;
  T pfact = 1.0;
  for (int i=lmm; i>0; i--)
    mfact *= static_cast<T>(i);
  for (int i=lpm; i>0; i--)
    pfact *= static_cast<T>(i);
  T prefactor = std::sqrt (static_cast<T>(2*l+1)*mfact/(4.0*M_PI*pfact));
  prefactor *= LegendrePlm (l, m, costheta);
  return std::complex<T>(prefactor*std::cos(m*phi), prefactor*std::sin(m*phi));
  //T prefactor = std::sqrt (static_cast<T>(2*l+1)*mfact/(4.0*M_PI*pfact));
  //T Plm = LegendrePlm (l, m, costheta);
  //std::complex<T> e2imphi (std::cos(m*phi), std::sin(m*phi));
  //return prefactor * Plm * e2imphi;
}
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 19:14:53 -0500 (Thu, 25 Apr 2013) $
 * $Id: Ylm.h 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/
