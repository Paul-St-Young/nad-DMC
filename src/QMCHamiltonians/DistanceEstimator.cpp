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
#include <Utilities/IteratorUtility.h>
#include <OhmmsData/AttributeSet.h>

namespace qmcplusplus
{

DistanceEstimator::DistanceEstimator(ParticleSet& elns,ParticleSet& ions)
  : elnPtcl(&elns), ionPtcl(&ions), hdf5_out(false), minExp(-1), maxExp(2)
{
  UpdateMode.set(COLLECTABLE,1);
  elnSpecies=elns.getSpeciesSet().getTotalNum();
  ionSpecies=ions.getSpeciesSet().getTotalNum();
  values.resize(maxExp-minExp);
  hdf5_out=false;
}

void DistanceEstimator::resetTargetParticleSet(ParticleSet& P)
{ // !!!! only allow reseting electron particle set for now
  elnPtcl = &P;
}

DistanceEstimator::Return_t DistanceEstimator::evaluate(ParticleSet& P)
{
  PosType vecr=ionPtcl->R[particle1Idx]-ionPtcl->R[particle2Idx];
  RealType r=std::sqrt( dot(vecr,vecr) );
  int exp=minExp;
  for (int i=0;i<values.size();i++){
    if (exp==0) exp++;
    values[i]=std::pow(r,(RealType) exp++);
  } 
  if (exp!=maxExp+1) APP_ABORT("DistanceEstimator::evaluate");
  return 0.0;
}

void DistanceEstimator::addObservables(PropertySetType& plist, BufferType& collectables)
{
  if(hdf5_out)
  {
  }
  else
  {
    myIndex=plist.size();
    for (int i=minExp; i<=maxExp; i++)
    {
      if (i!=0){
        std::stringstream sstr;
        sstr << "r" << particle1 << particle2 << "^" <<i;
        int id=plist.add(sstr.str());
      }
    }
  }
}

void DistanceEstimator::addObservables(PropertySetType& plist)
{
  myIndex=plist.size();
  for (int i=minExp; i<=maxExp; i++)
  {
    if (i!=0){
      std::stringstream sstr;
      sstr << "r^" <<i;
      int id=plist.add(sstr.str());
    }
  }
}

void DistanceEstimator::setObservables(PropertySetType& plist)
{
  if (!hdf5_out)
    std::copy(values.begin(),values.end(),plist.begin()+myIndex);
}

void DistanceEstimator::setParticlePropertyList(PropertySetType& plist, int offset)
{
  if (!hdf5_out)
    std::copy(values.begin(),values.end(),plist.begin()+myIndex+offset);
}

void DistanceEstimator::registerCollectables(vector<observable_helper*>& h5desc
                                       , hid_t gid) const
{
  if (hdf5_out)
  {
  /* save for reference
    vector<int> ndim(1,NumK);
    observable_helper* h5o=new observable_helper(myName);
    h5o->set_dimensions(ndim,myIndex);
    h5o->open(gid);
    h5desc.push_back(h5o);
    hsize_t kdims[2];
    kdims[0] = NumK;
    kdims[1] = OHMMS_DIM;
    string kpath = myName + "/kpoints";
    hid_t k_space = H5Screate_simple(2,kdims, NULL);
    hid_t k_set   = H5Dcreate (gid, kpath.c_str(), H5T_NATIVE_DOUBLE, k_space, H5P_DEFAULT);
    hid_t mem_space = H5Screate_simple (2, kdims, NULL);
    double *ptr = &(sourcePtcl->SK->KLists.kpts_cart[0][0]);
    herr_t ret = H5Dwrite(k_set, H5T_NATIVE_DOUBLE, mem_space, k_space, H5P_DEFAULT, ptr);
    H5Dclose (k_set);
    H5Sclose (mem_space);
    H5Sclose (k_space);
    H5Fflush(gid, H5F_SCOPE_GLOBAL);
  */
  } else {
  }
}

bool DistanceEstimator::put(xmlNodePtr cur)
{ // cur should look something like: 
// <estimator name="distance" particle1="C" particle2="H" minExp=-1 maxExp=2/>

  // process input line
  OhmmsAttributeSet attrib;
  attrib.add(particle1,"particle1");
  attrib.add(particle2,"particle2");
  attrib.add(minExp,"minExp");
  attrib.add(maxExp,"maxExp");
  attrib.put(cur);
  if (particle1.size()==0 or particle2.size()==0){
    APP_ABORT("DistanceEstimator::put - Please specify particle1 and particle2 in input xml");
  }// end if

  // find the specified particles
  SpeciesSet myions=ionPtcl->getSpeciesSet();
  for (int i=0;i<myions.size();i++){
    if (myions.speciesName[i]==particle1) particle1Idx=i;
    if (myions.speciesName[i]==particle2) particle2Idx=i;
  }

  return true;
}

bool DistanceEstimator::get(std::ostream& os) const
{
  return true;
}

QMCHamiltonianBase* DistanceEstimator::makeClone(ParticleSet& qp
    , TrialWaveFunction& psi)
{ // no cumulative variable, nothing needs to be saved
  DistanceEstimator* myclone = new DistanceEstimator(*this);
  myclone->hdf5_out=hdf5_out;
  myclone->myIndex=myIndex;
  return myclone;
}
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2945 $   $Date: 2008-08-05 10:21:33 -0500 (Tue, 05 Aug 2008) $
 * $Id: ForceBase.h 2945 2008-08-05 15:21:33Z jnkim $
 ***************************************************************************/
