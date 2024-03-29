//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim and Jordan Vincent
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include <math.h>
#include "SQD/SphericalPotential/ZOverRPotential.h"
#include "Numerics/RadialFunctorUtility.h"
namespace ohmmshf
{

/** constructor for the Nuclear Potential \f$V(r)=frac{Z}{r}\f$
 *
 * \param z the charge of the Nuclear Potential
 */
ZOverRPotential::ZOverRPotential(value_type z): Z(z)
{
  Qinfty=Z;
}

RadialPotentialBase::value_type
ZOverRPotential::evaluate(const BasisSetType& psi,
                          RadialOrbitalSet_t& V, int norb)
{
  if(!Vext)
  {
    Vext = new RadialOrbital_t(psi(0));
    integrand=new RadialOrbital_t(psi(0));
    for(int ig=0; ig < psi.m_grid->size(); ig++)
    {
      (*Vext)(ig) = -Z/psi.m_grid->r(ig);
    }
  }
  for(int ig=0; ig < psi.m_grid->size(); ig++)
  {
    value_type t = (*Vext)(ig);
    value_type sum = 0.0;
    for(int o=0; o < norb; o++)
    {
      V[o](ig) += t;
      sum += pow(psi(o,ig),2);
    }
    (*integrand)(ig) = t*sum;
  }
  return integrate_RK2(*integrand);
}

int ZOverRPotential::getNumOfNodes(int n, int l)
{
  MinEigenValue = -2.0*Z*Z/static_cast<value_type>(n*n);
  return n-l-1;
}
}
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 19:14:53 -0500 (Thu, 25 Apr 2013) $
 * $Id: ZOverRPotential.cpp 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/
