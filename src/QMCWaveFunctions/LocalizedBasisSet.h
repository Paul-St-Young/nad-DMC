//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
/** @file LocalizedBasisSet.h
 * @brief A derived class from BasisSetBase
 *
 * This is intended as a replacement for MolecularOrbitalBase and
 * any other localized basis set.
 */
#ifndef QMCPLUSPLUS_LOCALIZEDBASISSET_H
#define QMCPLUSPLUS_LOCALIZEDBASISSET_H

#include "QMCWaveFunctions/BasisSetBase.h"
#include "Particle/DistanceTable.h"

namespace qmcplusplus
{

/** A localized basis set derived from BasisSetBase<typename COT::value_type>
 *
 * This class performs the evaluation of the basis functions and their
 * derivatives for each of the N-particles in a configuration.
 * The template parameter COT denotes Centered-Orbital-Type which provides
 * a set of localized orbitals associated with a center.
 */
template<class COT>
struct LocalizedBasisSet: public BasisSetBase<typename COT::value_type>
{
  typedef COT                     ThisCOT_t;
  typedef typename COT::RadialOrbital_t    ThisRadialOrbital_t;
  typedef BasisSetBase<typename COT::value_type> BasisSetType;
  typedef typename BasisSetType::RealType      RealType;
  typedef typename BasisSetType::ValueType     ValueType;
  typedef typename BasisSetType::IndexType     IndexType;
  typedef typename BasisSetType::IndexVector_t IndexVector_t;
  typedef typename BasisSetType::ValueVector_t ValueVector_t;
  typedef typename BasisSetType::ValueMatrix_t ValueMatrix_t;
  typedef typename BasisSetType::GradVector_t  GradVector_t;
  typedef typename BasisSetType::GradMatrix_t  GradMatrix_t;

  using BasisSetType::ActivePtcl;
  using BasisSetType::Counter;
  using BasisSetType::BasisSetSize;
  using BasisSetType::Phi;
  using BasisSetType::dPhi;
  using BasisSetType::d2Phi;
  using BasisSetType::grad_grad_Phi;
  using BasisSetType::grad_grad_grad_Phi;
  using BasisSetType::Y;
  using BasisSetType::dY;
  using BasisSetType::d2Y;
  ///Reference to the center
  const ParticleSet& CenterSys;
  ///number of centers, e.g., ions
  int NumCenters;
  ///number of quantum particles
  int NumTargets;

  /** container to store the offsets of the basis functions
   *
   * the number of basis states for center J is BasisOffset[J+1]-Basis[J]
   */
  vector<int>  BasisOffset;

  /** container of the pointers to the Atomic Orbitals
   *
   * size of LOBasis =  number  of centers (e.g., ions)
   * AO[i] returns a Centered Orbital for an ion i
   */
  vector<COT*> LOBasis;

  /** container of the unique pointers to the Atomic Orbitals
   *
   * size of LOBasisSet = number  of unique centers
   */
  vector<COT*> LOBasisSet;

  /** distance table, e.g., ion-electron
   *
   * Localized basis sets require a pair relationship between CenterSys
   * and the quantum particle set.
   */
  const DistanceTableData* myTable;

  /** constructor
   * @param ions ionic system
   * @param els electronic system
   */
  LocalizedBasisSet(ParticleSet& ions, ParticleSet& els): CenterSys(ions), myTable(0)
  {
    myTable = DistanceTable::add(ions,els);
    NumCenters=CenterSys.getTotalNum();
    NumTargets=els.getTotalNum();
    LOBasis.resize(NumCenters,0);
    LOBasisSet.resize(CenterSys.getSpeciesSet().getTotalNum(),0);
    BasisOffset.resize(NumCenters+1);
  }

  LocalizedBasisSet<COT>* makeClone() const
  {
    LocalizedBasisSet<COT>* myclone = new LocalizedBasisSet<COT>(*this);
    for(int i=0; i<LOBasisSet.size(); ++i)
    {
      COT* cc=LOBasisSet[i]->makeClone();
      myclone->LOBasisSet[i]=cc;
      for(int j=0; j<CenterSys.getTotalNum(); ++j)
      {
        if(CenterSys.GroupID[j]==i)
          myclone->LOBasis[j]=cc;
      }
    }
    return myclone;
  }

  /**
   @param atable the distance table (ion-electron)
   @brief Assign the distance table (ion-electron) and
   *determine the total number of basis states.
  */
  void setBasisSetSize(int nbs)
  {
    if(nbs == BasisSetSize)
      return;
    if(myTable ==0)
    {
      app_error() << "LocalizedBasisSet cannot function without a distance table. Abort" << endl;
    }
    //reset the distance table for the atomic orbitals
    for(int i=0; i<LOBasisSet.size(); i++)
      LOBasisSet[i]->setTable(myTable);
    //evaluate the total basis dimension and offset for each center
    BasisOffset[0] = 0;
    for(int c=0; c<NumCenters; c++)
      BasisOffset[c+1] = BasisOffset[c]+LOBasis[c]->getBasisSetSize();
    BasisSetSize = BasisOffset[NumCenters];
    this->resize(NumTargets);
  }

  void resetParameters(const opt_variables_type& active)
  {
    //reset each unique basis functions
    for(int i=0; i<LOBasisSet.size(); i++)
      LOBasisSet[i]->resetParameters(active);
  }

  /** reset the distance table with a new target P
   */
  void resetTargetParticleSet(ParticleSet& P)
  {
    myTable = DistanceTable::add(CenterSys,P);
    for(int i=0; i<LOBasisSet.size(); i++)
      LOBasisSet[i]->setTable(myTable);
  }

  inline void
  evaluateWithHessian(const ParticleSet& P, int iat)
  {
    for(int c=0; c<NumCenters; c++)
      LOBasis[c]->evaluateForWalkerMove(c,iat,BasisOffset[c],Phi,dPhi,grad_grad_Phi);
    Counter++; // increment a conter
  }

  inline void
  evaluateWithThirdDeriv(const ParticleSet& P, int iat)
  {
    // should only work for s,p
    for(int c=0; c<NumCenters; c++)
      LOBasis[c]->evaluateForWalkerMove(c,iat,BasisOffset[c],Phi,dPhi,grad_grad_Phi,grad_grad_grad_Phi);
    Counter++; // increment a conter
  }

  inline void
  evaluateThirdDerivOnly(const ParticleSet& P, int iat)
  {
    // should only work for s,p
    for(int c=0; c<NumCenters; c++)
      LOBasis[c]->evaluateThirdDerivOnly(c,iat,BasisOffset[c],grad_grad_grad_Phi);
    Counter++; // increment a conter
  }

  inline void
  evaluateForWalkerMove(const ParticleSet& P)
  {
    for(int c=0; c<NumCenters; c++)
      LOBasis[c]->evaluateForWalkerMove(c,0,NumTargets,BasisOffset[c],Y,dY,d2Y);
    Counter++; // increment a conter
  }

  inline void
  evaluateForWalkerMove(const ParticleSet& P, int iat)
  {
    for(int c=0; c<NumCenters; c++)
      LOBasis[c]->evaluateForWalkerMove(c,iat,BasisOffset[c],Phi,dPhi,d2Phi);
    Counter++;
  }

  inline void
  evaluateForPtclMove(const ParticleSet& P, int iat)
  {
    for(int c=0; c<NumCenters; c++)
      LOBasis[c]->evaluateForPtclMove(c,iat,BasisOffset[c],Phi);
    Counter++;
    ActivePtcl=iat;
  }

  inline void
  evaluateAllForPtclMove(const ParticleSet& P, int iat)
  {
    for(int c=0; c<NumCenters; c++)
      LOBasis[c]->evaluateAllForPtclMove(c,iat,BasisOffset[c],Phi,dPhi,d2Phi);
    Counter++;
    ActivePtcl=iat;
  }

  inline void
  evaluateForPtclMoveWithHessian(const ParticleSet& P, int iat)
  {
    for(int c=0; c<NumCenters; c++)
      LOBasis[c]->evaluateAllForPtclMove(c,iat,BasisOffset[c],Phi,dPhi,grad_grad_Phi);
    Counter++;
    ActivePtcl=iat;
  }

  /** add a new set of Centered Atomic Orbitals
   * @param icenter the index of the center
   * @param aos a set of Centered Atomic Orbitals
   */
  void add(int icenter, COT* aos)
  {
    aos->setTable(myTable);
    LOBasisSet[icenter]=aos;
    for(int i=0; i<NumCenters; i++)
    {
      if(CenterSys.GroupID[i] == icenter)
        LOBasis[i]=aos;
    }
  }
};
}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 19:14:53 -0500 (Thu, 25 Apr 2013) $
 * $Id: LocalizedBasisSet.h 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/

