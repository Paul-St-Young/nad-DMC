//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
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
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "Utilities/IteratorUtility.h"
#include "LongRange/StructFact.h"
#include "Particle/HDFWalkerOutput.h"
#include <map>

#ifdef QMC_CUDA
#include "Particle/accept_kernel.h"
#endif

namespace qmcplusplus
{

/** store minimum Walker data for the next section
 */
struct MCSample
{
  typedef ParticleSet::Walker_t Walker_t;

  ParticleSet::ParticlePos_t R;
  ParticleSet::ParticlePos_t ionPos;
  ParticleSet::ParticleGradient_t G;
  ParticleSet::ParticleLaplacian_t L;
  ParticleSet::RealType LogPsi, Sign, PE, KE;

  inline MCSample(const Walker_t& w)
  :R(w.R), ionPos(w.ionPos), G(w.G), L(w.L)
  {
    LogPsi=w.Properties(LOGPSI);
    Sign=w.Properties(SIGN);
    PE=w.Properties(LOCALPOTENTIAL);
    KE=w.Properties(LOCALENERGY)-PE;
  }

  inline MCSample(int n)
  {
    R.resize(n);
    G.resize(n);
    L.resize(n);
    ionPos.resize(n);
  }

  inline MCSample(int n, int m)
  {
    R.resize(n);
    G.resize(n);
    L.resize(n);
    ionPos.resize(m);
  }

  inline void put(const Walker_t& w)
  {
    R=w.R;
    G=w.G;
    L=w.L;
    ionPos.resize(w.ionPos.size());
    ionPos=w.ionPos;
    LogPsi=w.Properties(LOGPSI);
    Sign=w.Properties(SIGN);
    PE=w.Properties(LOCALPOTENTIAL);
    KE=w.Properties(LOCALENERGY)-PE;
  }

  inline void get(Walker_t& w) const
  {
    w.R=R;
    w.G=G;
    w.L=L;
    w.ionPos.resize(ionPos.size());
    w.ionPos=ionPos;
    w.Properties(LOGPSI)=LogPsi;
    w.Properties(SIGN)=Sign;
    w.Properties(LOCALPOTENTIAL)=PE;
    w.Properties(LOCALENERGY)=PE+KE;
  }
};

MCWalkerConfiguration::MCWalkerConfiguration():
  OwnWalkers(true),ReadyForPbyP(false),UpdateMode(Update_Walker),Polymer(0),

  MaxSamples(10),CurSampleCount(0),GlobalNumWalkers(0),reptile(0)
#ifdef QMC_CUDA
  ,RList_GPU("MCWalkerConfiguration::RList_GPU"),
  GradList_GPU("MCWalkerConfiguration::GradList_GPU"),
  LapList_GPU("MCWalkerConfiguration::LapList_GPU"),
  Rnew_GPU("MCWalkerConfiguration::Rnew_GPU"),
  NLlist_GPU ("MCWalkerConfiguration::NLlist_GPU"),
  AcceptList_GPU("MCWalkerConfiguration::AcceptList_GPU"),
  iatList_GPU("iatList_GPU")
#endif
{
  //move to ParticleSet
  //initPropertyList();
}

MCWalkerConfiguration::MCWalkerConfiguration(const MCWalkerConfiguration& mcw)
  : ParticleSet(mcw), OwnWalkers(true), GlobalNumWalkers(mcw.GlobalNumWalkers),
    UpdateMode(Update_Walker), ReadyForPbyP(false), Polymer(0),
    MaxSamples(mcw.MaxSamples), CurSampleCount(0)
#ifdef QMC_CUDA
    ,RList_GPU("MCWalkerConfiguration::RList_GPU"),
    GradList_GPU("MCWalkerConfiguration::GradList_GPU"),
    LapList_GPU("MCWalkerConfiguration::LapList_GPU"),
    Rnew_GPU("MCWalkerConfiguration::Rnew_GPU"),
    NLlist_GPU ("MCWalkerConfiguration::NLlist_GPU"),
    AcceptList_GPU("MCWalkerConfiguration::AcceptList_GPU"),
    iatList_GPU("iatList_GPU")
#endif
{
  GlobalNumWalkers=mcw.GlobalNumWalkers;
  WalkerOffsets=mcw.WalkerOffsets;
  Properties=mcw.Properties;
  //initPropertyList();
}

///default destructor
MCWalkerConfiguration::~MCWalkerConfiguration()
{
  if(OwnWalkers)
    destroyWalkers(WalkerList.begin(), WalkerList.end());
}


void MCWalkerConfiguration::createWalkers(int n, int nions)
{
  if(WalkerList.empty())
  {
    while(n)
    {
      Walker_t* awalker=new Walker_t(GlobalNum,nions);
      awalker->R = R;

      WalkerList.push_back(awalker);
      --n;
    }
  }
  else
  {
    if(WalkerList.size()>=n)
    {
      int iw=WalkerList.size();//copy from the back
      for(int i=0; i<n; ++i)
      {
        WalkerList.push_back(new Walker_t(*WalkerList[--iw]));
      }
    }
    else
    {
      int nc=n/WalkerList.size();
      int nw0=WalkerList.size();
      for(int iw=0; iw<nw0; ++iw)
      {
        for(int ic=0; ic<nc; ++ic)
          WalkerList.push_back(new Walker_t(*WalkerList[iw]));
      }
      n-=nc*nw0;
      while(n>0)
      {
        WalkerList.push_back(new Walker_t(*WalkerList[--nw0]));
        --n;
      }
    }
  }
  resizeWalkerHistories();
}


void MCWalkerConfiguration::createWalkers(int n)
{
  if(WalkerList.empty())
  {
    while(n)
    {
      Walker_t* awalker=new Walker_t(GlobalNum);
      awalker->R = R;

      WalkerList.push_back(awalker);
      --n;
    }
  }
  else
  {
    if(WalkerList.size()>=n)
    {
      int iw=WalkerList.size();//copy from the back
      for(int i=0; i<n; ++i)
      {
        WalkerList.push_back(new Walker_t(*WalkerList[--iw]));
      }
    }
    else
    {
      int nc=n/WalkerList.size();
      int nw0=WalkerList.size();
      for(int iw=0; iw<nw0; ++iw)
      {
        for(int ic=0; ic<nc; ++ic)
          WalkerList.push_back(new Walker_t(*WalkerList[iw]));
      }
      n-=nc*nw0;
      while(n>0)
      {
        WalkerList.push_back(new Walker_t(*WalkerList[--nw0]));
        --n;
      }
    }
  }
  resizeWalkerHistories();
}



void MCWalkerConfiguration::resize(int numWalkers, int numPtcls)
{
  if(GlobalNum && WalkerList.size())
    app_warning() << "MCWalkerConfiguration::resize cleans up the walker list." << endl;
  ParticleSet::resize(unsigned(numPtcls));
  int dn=numWalkers-WalkerList.size();
  if(dn>0)
    createWalkers(dn);
  if(dn<0)
  {
    int nw=-dn;
    if(nw<WalkerList.size())
    {
      iterator it = WalkerList.begin();
      while(nw)
      {
        delete *it;
        ++it;
        --nw;
      }
      WalkerList.erase(WalkerList.begin(),WalkerList.begin()-dn);
    }
  }
  //iterator it = WalkerList.begin();
  //while(it != WalkerList.end()) {
  //  delete *it++;
  //}
  //WalkerList.erase(WalkerList.begin(),WalkerList.end());
  //R.resize(np);
  //GlobalNum = np;
  //createWalkers(nw);
}

///returns the next valid iterator
MCWalkerConfiguration::iterator
MCWalkerConfiguration::destroyWalkers(iterator first, iterator last)
{
  if(OwnWalkers)
  {
    iterator it = first;
    while(it != last)
    {
      delete *it++;
    }
  }
  return WalkerList.erase(first,last);
}

void MCWalkerConfiguration::createWalkers(iterator first, iterator last)
{
  destroyWalkers(WalkerList.begin(),WalkerList.end());
  OwnWalkers=true;
  while(first != last)
  {
    WalkerList.push_back(new Walker_t(**first));
    ++first;
  }
}

void
MCWalkerConfiguration::destroyWalkers(int nw)
{
  if(nw > WalkerList.size())
  {
    app_warning() << "  Cannot remove walkers. Current Walkers = " << WalkerList.size() << endl;
    return;
  }
  nw=WalkerList.size()-nw;
  int iw=nw;
  while(iw<WalkerList.size())
  {
    delete WalkerList[iw++];
  }
  //iterator it(WalkerList.begin()+nw),it_end(WalkerList.end());
  //while(it != it_end)
  //{
  //  delete *it++;
  //}
  WalkerList.erase(WalkerList.begin()+nw,WalkerList.end());
}

void MCWalkerConfiguration::copyWalkers(iterator first, iterator last, iterator it)
{
  while(first != last)
  {
    (*it++)->makeCopy(**first++);
  }
}


void
MCWalkerConfiguration::copyWalkerRefs(Walker_t* head, Walker_t* tail)
{
  if(OwnWalkers)
    //destroy the current walkers
  {
    destroyWalkers(WalkerList.begin(), WalkerList.end());
    WalkerList.clear();
    OwnWalkers=false;//set to false to prevent deleting the Walkers
  }
  if(WalkerList.size()<2)
  {
    WalkerList.push_back(0);
    WalkerList.push_back(0);
  }
  WalkerList[0]=head;
  WalkerList[1]=tail;
}

/** Make Metropolis move to the walkers and save in a temporary array.
 * @param it the iterator of the first walker to work on
 * @param tauinv  inverse of the time step
 *
 * R + D + X
 */
void MCWalkerConfiguration::sample(iterator it, RealType tauinv)
{
  APP_ABORT("MCWalkerConfiguration::sample obsolete");
//  makeGaussRandom(R);
//  R *= tauinv;
//  R += (*it)->R + (*it)->Drift;
}

void MCWalkerConfiguration::reset()
{
  iterator it(WalkerList.begin()), it_end(WalkerList.end());
  while(it != it_end)
    //(*it)->reset();++it;}
  {
    (*it)->Weight=1.0;
    (*it)->Multiplicity=1.0;
    ++it;
  }
}

//void MCWalkerConfiguration::clearAuxDataSet() {
//  UpdateMode=Update_Particle;
//  int nbytes=128*GlobalNum*sizeof(RealType);//could be pagesize
//  if(WalkerList.size())//check if capacity is bigger than the estimated one
//    nbytes = (WalkerList[0]->DataSet.capacity()>nbytes)?WalkerList[0]->DataSet.capacity():nbytes;
//  iterator it(WalkerList.begin());
//  iterator it_end(WalkerList.end());
//  while(it!=it_end) {
//    (*it)->DataSet.clear();
//    //CHECK THIS WITH INTEL 10.1
//    //(*it)->DataSet.reserve(nbytes);
//    ++it;
//  }
//  ReadyForPbyP = true;
//}
//
//bool MCWalkerConfiguration::createAuxDataSet(int nfield) {
//
//  if(ReadyForPbyP) return false;
//
//  ReadyForPbyP=true;
//  UpdateMode=Update_Particle;
//  iterator it(WalkerList.begin());
//  iterator it_end(WalkerList.end());
//  while(it!=it_end) {
//    (*it)->DataSet.reserve(nfield); ++it;
//  }
//
//  return true;
//}


/** reset the Property container of all the walkers
 */
void MCWalkerConfiguration::resetWalkerProperty(int ncopy)
{
  int m(PropertyList.size());
  app_log() << "  Resetting Properties of the walkers " << ncopy << " x " << m << endl;
  Properties.resize(ncopy,m);
  iterator it(WalkerList.begin()),it_end(WalkerList.end());
  while(it != it_end)
  {
    (*it)->resizeProperty(ncopy,m);
    (*it)->Weight=1;
    ++it;
  }
  resizeWalkerHistories();
}

void MCWalkerConfiguration::resizeWalkerHistories()
{
  //using vector<vector<RealType> > is too costly. 
  int np=PropertyHistory.size();
  if(np)
    for(int iw=0; iw<WalkerList.size(); ++iw)
      WalkerList[iw]->PropertyHistory=PropertyHistory;

  np=PHindex.size();
  if(np)
    for(int iw=0; iw<WalkerList.size(); ++iw)
      WalkerList[iw]->PHindex=PHindex;;
}

/** allocate the SampleStack
 * @param n number of samples per thread
 */
void MCWalkerConfiguration::setNumSamples(int n)
{
  clearEnsemble();
  MaxSamples=n;
  //do not add anything
  if(n==0)
    return;
  //SampleStack.resize(n,0);
  SampleStack.reserve(n);
  int nadd=n-SampleStack.size();
  while(nadd>0)
  {
    SampleStack.push_back(new MCSample(GlobalNum));
    --nadd;
  }
}


void MCWalkerConfiguration::setNumSamples(int n, int nions)
{
  clearEnsemble();
  MaxSamples=n;
  //do not add anything
  if(n==0)
    return;
  //SampleStack.resize(n,0);
  SampleStack.reserve(n);
  int nadd=n-SampleStack.size();
  while(nadd>0)
  {
    SampleStack.push_back(new MCSample(GlobalNum,nions));
    --nadd;
  }
}


/** save the current walkers to SampleStack
 */
void MCWalkerConfiguration::saveEnsemble()
{
  saveEnsemble(WalkerList.begin(),WalkerList.end());
}

/** save the [first,last) walkers to SampleStack
 */
void MCWalkerConfiguration::saveEnsemble(iterator first, iterator last)
{
  //safety check
  if(MaxSamples==0)
    return;
  while((first != last) && (CurSampleCount<MaxSamples))
  {
    SampleStack[CurSampleCount]->put(**first);
    ++first;
    ++CurSampleCount;
  }
}

/** load SampleStack to WalkerList
 */
void MCWalkerConfiguration::loadEnsemble()
{
  int nsamples=std::min(MaxSamples,CurSampleCount);
  if(SampleStack.empty() || nsamples==0)
    return;
  Walker_t::PropertyContainer_t prop(1,PropertyList.size());
  delete_iter(WalkerList.begin(),WalkerList.end());
  WalkerList.resize(nsamples);
  for(int i=0; i<nsamples; ++i)
  {
    Walker_t* awalker=new Walker_t(GlobalNum);
    awalker->Properties.copy(prop);
    SampleStack[i]->get(*awalker);
    WalkerList[i]=awalker;
  }
  resizeWalkerHistories();
  clearEnsemble();
}

///** load SampleStack to WalkerList
// */
//void MCWalkerConfiguration::loadEnsemble(const Walker_t& wcopy)
//{
//  int nsamples=std::min(MaxSamples,CurSampleCount);
//  if(SampleStack.empty() || nsamples==0) return;
//
//  Walker_t::PropertyContainer_t prop(1,PropertyList.size());
//
//  while(WalkerList.size()) pop_back();
//  WalkerList.resize(nsamples);
//
//  for(int i=0; i<nsamples; ++i)
//  {
//    Walker_t* awalker=new Walker_t(GlobalNum);
//    awalker->Properties.copy(prop);
//    SampleStack[i]->get(*awalker);
//    WalkerList[i]=awalker;
//  }
//  resizeWalkerHistories();
//  clearEnsemble();
//}
//
//void MCWalkerConfiguration::loadEnsemble(MCWalkerConfiguration& other)
//{
//  if(SampleStack.empty()) return;
//
//  Walker_t twalker(*WalkerList[0]);
//  for(int i=0; i<MaxSamples; ++i)
//  {
//    Walker_t* awalker=new Walker_t(twalker);
//    SampleStack[i]->get(*awalker);
//    other.WalkerList.push_back(awalker);
//  }
//
//  clearEnsemble();
//}

bool
MCWalkerConfiguration::dumpEnsemble(std::vector<MCWalkerConfiguration*>& others
                                    , HDFWalkerOutput* out, int np)
{
#if !(defined(__bgp__)||(__bgq__))
  MCWalkerConfiguration wtemp;
  wtemp.resize(0,GlobalNum);
  wtemp.loadEnsemble(others,false);
  int w=wtemp.getActiveWalkers();
  if(w==0)
    return false;
  vector<int> nwoff(np+1,0);
  for(int ip=0; ip<np; ++ip)
    nwoff[ip+1]=nwoff[ip]+w;
  wtemp.setGlobalNumWalkers(nwoff[np]);
  wtemp.setWalkerOffsets(nwoff);
  out->dump(wtemp);
#endif
  return true;
}

void MCWalkerConfiguration::loadEnsemble(std::vector<MCWalkerConfiguration*>& others, bool doclean)
{
  vector<int> off(others.size()+1,0);
  for(int i=0; i<others.size(); ++i)
  {
    off[i+1]=off[i]+std::min(others[i]->MaxSamples,others[i]->CurSampleCount);
  }
  int nw_tot=off.back();
  if(nw_tot)
  {
    Walker_t::PropertyContainer_t prop(1,PropertyList.size());
    while(WalkerList.size())
      pop_back();
    WalkerList.resize(nw_tot);
    for(int i=0; i<others.size(); ++i)
    {
      vector<MCSample*>& astack(others[i]->SampleStack);
      for(int j=0, iw=off[i]; iw<off[i+1]; ++j, ++iw)
      {
        Walker_t* awalker=new Walker_t(GlobalNum);
        awalker->Properties.copy(prop);
        astack[j]->get(*awalker);
        WalkerList[iw]=awalker;
      }
      if(doclean)
        others[i]->clearEnsemble();
    }
  }
  if(doclean)
    resizeWalkerHistories();
}

void MCWalkerConfiguration::loadEnsemble(std::vector<MCWalkerConfiguration*>& others, int nions, bool doclean)
{
  vector<int> off(others.size()+1,0);
  for(int i=0; i<others.size(); ++i)
  {
    off[i+1]=off[i]+std::min(others[i]->MaxSamples,others[i]->CurSampleCount);
  }
  int nw_tot=off.back();
  if(nw_tot)
  {
    Walker_t::PropertyContainer_t prop(1,PropertyList.size());
    while(WalkerList.size())
      pop_back();
    WalkerList.resize(nw_tot);
    for(int i=0; i<others.size(); ++i)
    {
      vector<MCSample*>& astack(others[i]->SampleStack);
      for(int j=0, iw=off[i]; iw<off[i+1]; ++j, ++iw)
      {
        Walker_t* awalker=new Walker_t(GlobalNum,nions);
        awalker->Properties.copy(prop);
        astack[j]->get(*awalker);
        WalkerList[iw]=awalker;
      }
      if(doclean)
        others[i]->clearEnsemble();
    }
  }
  if(doclean)
    resizeWalkerHistories();
}


void MCWalkerConfiguration::clearEnsemble()
{
  //delete_iter(SampleStack.begin(),SampleStack.end());
  for(int i=0; i<SampleStack.size(); ++i)
    if(SampleStack[i])
      delete SampleStack[i];
  SampleStack.clear();
  MaxSamples=0;
  CurSampleCount=0;
}

bool MCWalkerConfiguration::swap_pos(iterator it, iterator it_end)
{ // !!!!!!!!!!!!!!! hack !!!!!!!!!!!!!!!!!
  PosType ip,ep;
  for(; it != it_end; ++it)
  { // for entire walker population, trick QMCPACK into dumping ion positions by
      // swapping the first few electron positions with ions
       // This will work as long as there are no fewer electrons than ions

    Walker_t& thisWalker(**it);

    for (int ionIdx=0;ionIdx<ionPos.size();ionIdx++){
        ip =  thisWalker.ionPos[ionIdx];
        ep =  thisWalker.R[ionIdx];

        thisWalker.ionPos[ionIdx] = ep;
        thisWalker.R[ionIdx]  = ip;
    }
  }
   return 1;
}

#ifdef QMC_CUDA
void MCWalkerConfiguration::updateLists_GPU()
{
  int nw = WalkerList.size();
  int NumSpecies = getSpeciesSet().TotalNum;
  if (Rnew_GPU.size() != nw)
  {
    Rnew_GPU.resize(nw);
    RhokLists_GPU.resize(NumSpecies);
    for (int isp=0; isp<NumSpecies; isp++)
      RhokLists_GPU[isp].resize(nw);
    Rnew_host.resize(nw);
    Rnew.resize(nw);
    AcceptList_GPU.resize(nw);
    AcceptList_host.resize(nw);
    RList_GPU.resize(nw);
    GradList_GPU.resize(nw);
    LapList_GPU.resize(nw);
    DataList_GPU.resize(nw);
  }
  hostlist.resize(nw);
  for (int iw=0; iw<nw; iw++)
  {
    if (WalkerList[iw]->R_GPU.size() != R.size())
      cerr << "Error in R_GPU size for iw = " << iw << "!\n";
    hostlist[iw] = (CUDA_PRECISION*)WalkerList[iw]->R_GPU.data();
  }
  RList_GPU = hostlist;
  for (int iw=0; iw<nw; iw++)
  {
    if (WalkerList[iw]->Grad_GPU.size() != R.size())
      cerr << "Error in Grad_GPU size for iw = " << iw << "!\n";
    hostlist[iw] = (CUDA_PRECISION*)WalkerList[iw]->Grad_GPU.data();
  }
  GradList_GPU = hostlist;
  for (int iw=0; iw<nw; iw++)
  {
    if (WalkerList[iw]->Lap_GPU.size() != R.size())
      cerr << "Error in Lap_GPU size for iw = " << iw << "!\n";
    hostlist[iw] = (CUDA_PRECISION*)WalkerList[iw]->Lap_GPU.data();
  }
  LapList_GPU = hostlist;
  for (int iw=0; iw<nw; iw++)
    hostlist[iw] = WalkerList[iw]->cuda_DataSet.data();
  DataList_GPU = hostlist;
  for (int isp=0; isp<NumSpecies; isp++)
  {
    for (int iw=0; iw<nw; iw++)
      hostlist[iw] = WalkerList[iw]->get_rhok_ptr(isp);
    RhokLists_GPU[isp] = hostlist;
  }
}

void
MCWalkerConfiguration::allocateGPU(size_t buffersize)
{
  int N = WalkerList[0]->R.size();
  int Numk = 0;
  if (SK)
    Numk = SK->KLists.numk;
  int NumSpecies = getSpeciesSet().TotalNum;
  for (int iw=0; iw<WalkerList.size(); iw++)
  {
    Walker_t &walker = *(WalkerList[iw]);
    walker.resizeCuda(buffersize, NumSpecies, Numk);
  }
}



void MCWalkerConfiguration::copyWalkersToGPU(bool copyGrad)
{
  R_host.resize(WalkerList[0]->R.size());
  for (int iw=0; iw<WalkerList.size(); iw++)
  {
    for (int i=0; i<WalkerList[iw]->size(); i++)
      for (int dim=0; dim<OHMMS_DIM; dim++)
        R_host[i][dim] = WalkerList[iw]->R[i][dim];
    WalkerList[iw]->R_GPU = R_host;
  }
  if (copyGrad)
    copyWalkerGradToGPU();
}


void MCWalkerConfiguration::copyWalkerGradToGPU()
{
  R_host.resize(WalkerList[0]->R.size());
  for (int iw=0; iw<WalkerList.size(); iw++)
  {
    for (int i=0; i<WalkerList[iw]->size(); i++)
      for (int dim=0; dim<OHMMS_DIM; dim++)
        R_host[i][dim] = WalkerList[iw]->G[i][dim];
    WalkerList[iw]->Grad_GPU = R_host;
  }
}


void MCWalkerConfiguration::proposeMove_GPU
(vector<PosType> &newPos, int iat)
{
  if (Rnew_host.size() < newPos.size())
    Rnew_host.resize(newPos.size());
  for (int i=0; i<newPos.size(); i++)
    for (int dim=0; dim<OHMMS_DIM; dim++)
      Rnew_host[i][dim] = newPos[i][dim];
  Rnew_GPU.asyncCopy(Rnew_host);
  Rnew = newPos;
  CurrentParticle = iat;
}


void MCWalkerConfiguration::acceptMove_GPU(vector<bool> &toAccept)
{
  if (AcceptList_host.size() < toAccept.size())
    AcceptList_host.resize(toAccept.size());
  for (int i=0; i<toAccept.size(); i++)
    AcceptList_host[i] = (int)toAccept[i];
  AcceptList_GPU.asyncCopy(AcceptList_host);
//   app_log() << "toAccept.size()        = " << toAccept.size() << endl;
//   app_log() << "AcceptList_host.size() = " << AcceptList_host.size() << endl;
//   app_log() << "AcceptList_GPU.size()  = " << AcceptList_GPU.size() << endl;
//   app_log() << "WalkerList.size()      = " << WalkerList.size() << endl;
//   app_log() << "Rnew_GPU.size()        = " << Rnew_GPU.size() << endl;
//   app_log() << "RList_GPU.size()       = " << RList_GPU.size() << endl;
  if (RList_GPU.size() != WalkerList.size())
    cerr << "Error in RList_GPU size.\n";
  if (Rnew_GPU.size() != WalkerList.size())
    cerr << "Error in Rnew_GPU size.\n";
  if (AcceptList_GPU.size() != WalkerList.size())
    cerr << "Error in AcceptList_GPU_GPU size.\n";
  accept_move_GPU_cuda
  (RList_GPU.data(), (CUDA_PRECISION*)Rnew_GPU.data(),
   AcceptList_GPU.data(), CurrentParticle, WalkerList.size());
}



void MCWalkerConfiguration::NLMove_GPU(vector<Walker_t*> &walkers,
                                       vector<PosType> &newpos,
                                       vector<int> &iat)
{
  int N = walkers.size();
  if (NLlist_GPU.size() < N)
  {
    NLlist_GPU.resize(N);
    NLlist_host.resize(N);
  }
  if (Rnew_GPU.size() < N)
  {
    Rnew_host.resize(N);
    Rnew_GPU.resize(N);
  }
  for (int iw=0; iw<N; iw++)
  {
    Rnew_host[iw]  = newpos[iw];
    NLlist_host[iw] =
      (CUDA_PRECISION*)(walkers[iw]->R_GPU.data()) + OHMMS_DIM*iat[iw];
  }
  Rnew_GPU   = Rnew_host;
  NLlist_GPU = NLlist_host;
  NL_move_cuda (NLlist_GPU.data(), (CUDA_PRECISION*)Rnew_GPU.data(), N);
}
#endif

}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 5826 $   $Date: 2013-05-06 20:02:09 -0500 (Mon, 06 May 2013) $
 * $Id: MCWalkerConfiguration.cpp 5826 2013-05-07 01:02:09Z jnkim $
 ***************************************************************************/
