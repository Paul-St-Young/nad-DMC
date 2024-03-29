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
#ifndef OHMMSHF_RADIALPOTENTIALBASE_H
#define OHMMSHF_RADIALPOTENTIALBASE_H

#include "AtomicHF/HFAtomicOrbitals.h"
#include "Optimize/VarList.h"
/*! \author Jeongnim Kim
 *  \author Jordan Vincent
 *  \note  The original Prim was written in F90 by Tim Wilkens.
 */

class Clebsch_Gordan;

namespace ohmmshf
{

/**class RadialPotentialBase
 *\brief An abstract class to build grid potentials
 */
struct RadialPotentialBase
{

  typedef HFAtomicOrbitals::value_type         value_type;
  typedef HFAtomicOrbitals::RadialGrid_t       RadialGrid_t;
  typedef HFAtomicOrbitals::RadialOrbital_t    RadialOrbital_t;
  typedef HFAtomicOrbitals::RadialOrbitalSet_t RadialOrbitalSet_t;

  ///constructor
  RadialPotentialBase() {}

  ///destructor
  virtual ~RadialPotentialBase() { }

  /*! \fn virtual
    value_type evaluate(const HFAtomicOrbitals& mo,
    RadialOrbitalSet_t& V, int norb)
    * \param psi the wavefunction
    * \param V the potential
    * \param norb the number of orbitals

   */
  virtual
  value_type evaluate(const HFAtomicOrbitals& mo,
                      RadialOrbitalSet_t& V, int norb) = 0;
};

/**class ZOverRFunctor
 *\brief Implements the Nuclear potential
 */
struct ZOverRFunctor: public RadialPotentialBase
{
  value_type Z;
  ZOverRFunctor(value_type z);
  value_type evaluate(const HFAtomicOrbitals& mo,
                      RadialOrbitalSet_t& V, int norb);
};

/**class HarmonicFunctor
 *\brief Implements the Harmonic potential
 */
struct HarmonicFunctor: public RadialPotentialBase
{
  value_type Omega;
  HarmonicFunctor(value_type omega);
  value_type evaluate(const HFAtomicOrbitals& mo,
                      RadialOrbitalSet_t& V, int norb);
};

/**class PseudoPotential
  *\brief Implements the Starkloff-Joanappoulos
  pseudopotential
  */
struct PseudoPotential: public RadialPotentialBase
{
  value_type Zeff, SJ_lambda, rc;
  PseudoPotential(VarRegistry<value_type>&,value_type,
                  value_type,value_type);
  PseudoPotential(value_type,value_type,value_type);
  value_type evaluate(const HFAtomicOrbitals& mo,
                      RadialOrbitalSet_t&V, int norb);
};


/**class HartreePotential
 *\brief Implements the Hartree potential
 */
struct HartreePotential: public RadialPotentialBase
{
  Clebsch_Gordan *CG_coeff;
  HartreePotential(Clebsch_Gordan*);
  value_type evaluate(const HFAtomicOrbitals& mo,
                      RadialOrbitalSet_t& V, int norb);
};

/**class ExchangePotential
  *\brief Implements the exchange potential
  */
struct ExchangePotential: public RadialPotentialBase
{
  Clebsch_Gordan *CG_coeff;
  ExchangePotential(Clebsch_Gordan*);
  value_type evaluate(const HFAtomicOrbitals& mo,
                      RadialOrbitalSet_t& V, int norb);
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 19:14:53 -0500 (Thu, 25 Apr 2013) $
 * $Id: RadialPotential.h 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/


