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
#include "SQD/SphericalPotential/SJPseudoPotential.h"
#include "Numerics/RadialFunctorUtility.h"
namespace ohmmshf
{

/**
 * \param vreg a registry for the optimizable variables \f$r_c\f$
 * and \f$\lambda\f$
 * \param zeff \f$Z_{Eff}\f$ the effective charge of the core
 * \param r_c \f$r_{c}\f$ the core-radius
 * \param lambda \f$\lambda\f$ the decay parameter
 * \brief The Constructor for the Starkloff-Joannopoulos
 * pseudopotential
 *
 */

SJPseudoPotential::SJPseudoPotential(VarRegistry<value_type>& vreg,
                                     value_type zeff, value_type r_c,
                                     value_type lambda):
  Zeff(zeff), rc(r_c), SJ_lambda(lambda)
{
  vreg["SJ_lambda"]=SJ_lambda;
  vreg["r_core"]=rc;
  //vreg.add("SJ_lambda",&SJ_lambda);
  //vreg.add("r_core",&rc);
}

/**
 * \param zeff \f$Z_{Eff}\f$ the effective charge of the core
 * \param r_c \f$r_{c}\f$ the core-radius
 * \param lambda \f$\lambda\f$ the decay parameter
 * \brief The Constructor for the Starkloff-Joannopoulos
 pseudopotential
 *
 */

SJPseudoPotential::SJPseudoPotential(value_type zeff, value_type r_c,
                                     value_type lambda):
  Zeff(zeff), rc(r_c), SJ_lambda(lambda) { }

RadialPotentialBase::value_type
SJPseudoPotential::evaluate(const BasisSetType& psi,
                            RadialOrbitalSet_t& V, int norb)
{
  RadialOrbital_t integrand(psi(0));
  for(int ig=0; ig < psi.m_grid->size(); ig++)
  {
    value_type r = psi.m_grid->r(ig);
    value_type SJ_num = 1.0-exp(-SJ_lambda*r);
    value_type SJ_den = 1.0+exp(-SJ_lambda*(r-rc));
    value_type v = -Zeff/r*SJ_num/SJ_den;
    value_type sum = 0.0;
    for(int o=0; o < norb; o++)
    {
      V[o](ig) += v;
      sum += pow(psi(o,ig),2);
    }
    integrand(ig) = v*sum;
  }
  return integrate_RK2(integrand);
}

int SJPseudoPotential::getNumOfNodes(int n, int l)
{
  //warning, the number of nodes should always be 1
  MinEigenValue = -(2.0*Zeff)/rc;
  return n-l-1;
}

}
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 19:14:53 -0500 (Thu, 25 Apr 2013) $
 * $Id: SJPseudoPotential.cpp 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/
