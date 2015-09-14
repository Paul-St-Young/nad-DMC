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
#include <QMCHamiltonians/DistanceEstimator.h>
#include <Particle/DistanceTableData.h>
#include <Particle/DistanceTable.h>
#include <OhmmsData/AttributeSet.h>
#include <Utilities/SimpleParser.h>
#include <set>

namespace qmcplusplus
{

DistanceEstimator::DistanceEstimator(ParticleSet& elns, ParticleSet& srcs): ions(srcs)
{
  // set maximum distance use the simulation cell radius if any direction is periodic
  // otherwise should be given with rcut="Dmax" in initialization of estimator
  if(elns.Lattice.SuperCellEnum)
    Dmax=elns.Lattice.SimulationCellRadius;

  // get the ion distance table, add names particle pairs
  d_table = DistanceTable::add(ions,ions);
  const SpeciesSet& species(srcs.getSpeciesSet());
  int ng=species.size();
  nag=srcs.getTotalNum();
  for(int i=0; i<ng; ++i)
  {
    for(int j=i+1; j<ng; j++)
    {
      std::stringstream nm;
      nm<<species.speciesName[i]<<"_"<<species.speciesName[j];
      names.push_back(nm.str());
    }
  }
  dist.resize(names.size());
}

void DistanceEstimator::resetTargetParticleSet(ParticleSet& P)
{
  d_table = DistanceTable::add(ions,ions);
}

DistanceEstimator::Return_t DistanceEstimator::evaluate(ParticleSet& P)
{
  for(int iat=0; iat<nag; ++iat)
  {
    int j(0);
    for(int nn=d_table->M[iat]; nn<d_table->M[iat+1]; ++nn)
    {
      RealType r=d_table->r(nn);
      if(r>=Dmax)
        continue;
      dist[j++]=r;
    }
  }
  return 0.0;
}

void DistanceEstimator::registerCollectables(vector<observable_helper*>& h5list
    , hid_t gid) const
{
}


void DistanceEstimator::addObservables(PropertySetType& plist, BufferType& collectables)
{
  addObservables(plist);
}


bool DistanceEstimator::put(xmlNodePtr cur)
{
  OhmmsAttributeSet attrib;
  attrib.add(Dmax,"rcut");
  attrib.put(cur);
  return true;
}

bool DistanceEstimator::get(std::ostream& os) const
{
  os << myName << " rcut=" << Dmax << endl;
  return true;
}

QMCHamiltonianBase* DistanceEstimator::makeClone(ParticleSet& qp
    , TrialWaveFunction& psi)
{
  //default constructor is sufficient
  DistanceEstimator* myClone = new DistanceEstimator(*this);
  myClone->Dmax=Dmax;
  myClone->resetTargetParticleSet(qp);
  return myClone;
}

}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2945 $   $Date: 2008-08-05 10:21:33 -0500 (Tue, 05 Aug 2008) $
 * $Id: ForceBase.h 2945 2008-08-05 15:21:33Z jnkim $
 ***************************************************************************/
