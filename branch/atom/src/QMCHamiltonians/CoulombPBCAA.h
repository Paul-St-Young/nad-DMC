//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim and Kris Delaney
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
#ifndef QMCPLUSPLUS_COULOMBPBCAA_TEMP_H
#define QMCPLUSPLUS_COULOMBPBCAA_TEMP_H
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMCHamiltonians/ForceBase.h"
#include "LongRange/LRCoulombSingleton.h"

namespace qmcplusplus
{

/** @ingroup hamiltonian
 *\brief Calculates the AA Coulomb potential using PBCs
 *
 * Functionally identical to CoulombPBCAA but uses a templated version of
 * LRHandler.
 */
struct CoulombPBCAA: public QMCHamiltonianBase, public ForceBase
{

  typedef LRCoulombSingleton::LRHandlerType LRHandlerType;
  typedef LRCoulombSingleton::GridType       GridType;
  typedef LRCoulombSingleton::RadFunctorType RadFunctorType;
  LRHandlerType* AA;
  GridType* myGrid;
  RadFunctorType* rVs;

  bool is_active;
  bool FirstTime;
  int SourceID;
  int NumSpecies;
  int ChargeAttribIndx;
  int MemberAttribIndx;
  int NumCenters;
  RealType myConst;
  RealType myRcut;
  string PtclRefName;
  vector<RealType> Zat,Zspec;
  vector<int> NofSpecies;
  vector<int> SpeciesID;

  Matrix<RealType> SR2;
  Vector<RealType> dSR;
  Vector<ComplexType> del_eikr;
  /// Flag for whether to compute forces or not
  bool ComputeForces;
//     madelung constant
  RealType MC0;


  /** constructor */
  CoulombPBCAA(ParticleSet& ref, bool active,
                   bool computeForces=false);

  ~CoulombPBCAA();

  void resetTargetParticleSet(ParticleSet& P);

  Return_t evaluate(ParticleSet& P);

  void update_source(ParticleSet& s);

  inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }
  Return_t registerData(ParticleSet& P, BufferType& buffer);
  Return_t updateBuffer(ParticleSet& P, BufferType& buffer);
  void copyFromBuffer(ParticleSet& P, BufferType& buf);
  void copyToBuffer(ParticleSet& P, BufferType& buf);
  Return_t evaluatePbyP(ParticleSet& P, int iat);
  void acceptMove(int iat);

  /** Do nothing */
  bool put(xmlNodePtr cur)
  {
    return true;
  }

  bool get(std::ostream& os) const
  {
    os << "CoulombPBCAA potential: " << PtclRefName;
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

  void initBreakup(ParticleSet& P);

  Return_t evalSR(ParticleSet& P);
  Return_t evalLR(ParticleSet& P);
  Return_t evalSRwithForces(ParticleSet& P);
  Return_t evalLRwithForces(ParticleSet& P);
  Return_t evalConsts();
  Return_t evaluateForPbyP(ParticleSet& P);

  void addObservables(PropertySetType& plist, BufferType& collectables);

  void setObservables(PropertySetType& plist)
  {
    QMCHamiltonianBase::setObservables(plist);
    if (ComputeForces)
      setObservablesF(plist);
  }

  void setParticlePropertyList(PropertySetType& plist, int offset)
  {
    QMCHamiltonianBase::setParticlePropertyList(plist, offset);
    if (ComputeForces)
      setParticleSetF(plist, offset);
  }

};

}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 5915 $   $Date: 2013-08-05 08:11:42 -0500 (Mon, 05 Aug 2013) $
 * $Id: CoulombPBCAA.h 5915 2013-08-05 13:11:42Z jnkim $
 ***************************************************************************/

