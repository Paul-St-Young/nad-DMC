//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeongnim Kim
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

/** @file DistanceEstimator.h
 * @brief Declare DistanceEstimator, indistinguishable particle distance 
 * is NOT yet supported. Also, no hdf5 support yet
 */
#ifndef QMCPLUSPLUS_DISTANCE_ESTIMATOR_H
#define QMCPLUSPLUS_DISTANCE_ESTIMATOR_H
#include <QMCHamiltonians/QMCHamiltonianBase.h>
namespace qmcplusplus
{

/** DistanceEstimator evaluate the distance between two particles
 *
 * <estimator name="distance" particle1="C" particle2="H" minExp=-1 maxExp=2/>
 */
class DistanceEstimator: public QMCHamiltonianBase
{
public:

  DistanceEstimator(ParticleSet& elns, ParticleSet& ions);

  void resetTargetParticleSet(ParticleSet& P);

  Return_t evaluate(ParticleSet& P);

  inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  void addObservables(PropertySetType& plist);
  void addObservables(PropertySetType& plist, BufferType& collectables);
  void registerCollectables(vector<observable_helper*>& h5desc, hid_t gid) const ;
  void setObservables(PropertySetType& plist);
  void setParticlePropertyList(PropertySetType& plist, int offset);
  bool put(xmlNodePtr cur);
  bool get(std::ostream& os) const;
  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

protected:
  ParticleSet *ionPtcl; // ion particle set
  ParticleSet *elnPtcl; // electron particle set

  /** number of species */
  int elnSpecies;
  int ionSpecies;

  string particle1,particle2; // r = distance(particle1,particle2)
  int particle1Idx,particle2Idx; // ionPtcl.R[particle1Idx] should be particle1's position
  int minExp,maxExp; // calculates r^{minExp}, r^{minExp+1},...,r^{maxExp}
  Vector<RealType> values;
  bool hdf5_out;
};

}
#endif
