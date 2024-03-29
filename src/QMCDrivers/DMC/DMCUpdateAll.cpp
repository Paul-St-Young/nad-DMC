//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCDrivers/DMC/DMCUpdateAll.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCDrivers/DriftOperators.h"

namespace qmcplusplus
{

/// Constructor.
DMCUpdateAllWithRejection::DMCUpdateAllWithRejection(MCWalkerConfiguration& w,
    TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator_t& rg):
  QMCUpdateBase(w,psi,h,rg)
{ 
  UpdatePbyP=false;
}

/// destructor
DMCUpdateAllWithRejection::~DMCUpdateAllWithRejection() { }

//void DMCUpdateAllWithRejection::initWalkers(WalkerIter_t it, WalkerIter_t it_end){
//}

//void DMCUpdateAllWithRejection::updateWalkers(WalkerIter_t it, WalkerIter_t it_end){
//}

/** advance all the walkers with killnode==no
 * @param nat number of particles to move
 *
 * When killnode==no, any move resulting in node-crossing is treated
 * as a normal rejection.
 */
void DMCUpdateAllWithRejection::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end,
    bool measure)
{
  for(; it != it_end; ++it)
  {
    Walker_t& thisWalker(**it);
    W.loadWalker(thisWalker,false);
    //create a 3N-Dimensional Gaussian with variance=1
    RealType nodecorr=setScaledDriftPbyPandNodeCorr(Tau,MassInvP,W.G,drift);
    //RealType nodecorr = setScaledDriftPbyPandNodeCorr(m_tauovermass,W.G,drift);
    makeGaussRandomWithEngine(deltaR,RandomGen);
    //if(!W.makeMoveWithDrift(thisWalker,drift,deltaR, m_sqrttau))
    if (!W.makeMoveWithDrift(thisWalker,drift ,deltaR,SqrtTauOverMass))
    {
      H.rejectedMove(W,thisWalker);
      continue;
    }
    //save old local energy
    RealType eold    = thisWalker.Properties(LOCALENERGY);
    RealType signold = thisWalker.Properties(SIGN);
    RealType enew  = eold;
    //evaluate wave functior
    RealType logpsi(Psi.evaluateLog(W));
    if(UseTMove)
      nonLocalOps.reset();
    bool accepted=false;
    RealType rr_accepted = 0.0;
    nodecorr=0.0;
    if(branchEngine->phaseChanged(Psi.getPhaseDiff()))
    {
      thisWalker.Age++;
      H.rejectedMove(W,thisWalker);
    }
    else
    {
      if(UseTMove)
        enew=H.evaluate(W,nonLocalOps.Txy);
      else
        enew=H.evaluate(W);

      RealType logGf = -0.5*Dot(deltaR,deltaR);
      //RealType nodecorr = setScaledDriftPbyPandNodeCorr(m_tauovermass,W.G,drift);
      RealType nodecorr=setScaledDriftPbyPandNodeCorr(Tau,MassInvP,W.G,drift);
      deltaR = thisWalker.R - W.R - drift;
      RealType logGb=logBackwardGF(deltaR);
      //RealType logGb = -m_oneover2tau*Dot(deltaR,deltaR);
      RealType prob= std::min(std::exp(logGb-logGf +2.0*(logpsi-thisWalker.Properties(LOGPSI))),1.0);
      //calculate rr_proposed here
      deltaR = W.R-thisWalker.R;
      RealType rr_proposed = Dot(deltaR,deltaR);
      if(RandomGen() > prob)
      {
        thisWalker.Age++;
        enew=eold;
        thisWalker.Properties(R2ACCEPTED)=0.0;
        thisWalker.Properties(R2PROPOSED)=rr_proposed;
        H.rejectedMove(W,thisWalker);
      }
      else
      {
        accepted=true;
        thisWalker.Age=0;
        W.saveWalker(thisWalker);
        rr_accepted = rr_proposed;
        thisWalker.resetProperty(logpsi,Psi.getPhase(),enew,rr_accepted,rr_proposed,nodecorr);
        H.auxHevaluate(W,thisWalker);
        H.saveProperty(thisWalker.getPropertyBase());
      }
    }
    if(UseTMove)
    {
      int ibar=nonLocalOps.selectMove(RandomGen());
      //make a non-local move
      if(ibar)
      {
        int iat=nonLocalOps.id(ibar);
        W.R[iat] += nonLocalOps.delta(ibar);
        W.update();
        logpsi=Psi.evaluateLog(W);
        setScaledDrift(Tau,W.G,drift);
        thisWalker.resetProperty(logpsi,Psi.getPhase(),eold);
        thisWalker.R[iat] = W.R[iat];
        ++NonLocalMoveAccepted;
      }
    }
    thisWalker.Weight *= branchEngine->branchWeight(enew,eold);
    //branchEngine->accumulate(eold,1);
    if(accepted)
      ++nAccept;
    else
      ++nReject;
  }
}

/*
 * DMCUpdateAllWithKill member functions
 */
/// Constructor.
DMCUpdateAllWithKill::DMCUpdateAllWithKill(MCWalkerConfiguration& w,
    TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator_t& rg): QMCUpdateBase(w,psi,h,rg)
{ 
  UpdatePbyP=false;
}

/// destructor
DMCUpdateAllWithKill::~DMCUpdateAllWithKill() { }

//void DMCUpdateAllWithKill::initWalkers(WalkerIter_t it, WalkerIter_t it_end){
//}

//void DMCUpdateAllWithKill::updateWalkers(WalkerIter_t it, WalkerIter_t it_end){
//}
/** advance all the walkers with killnode==yes
 */
void DMCUpdateAllWithKill::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end,
    bool measure)
{
  for(; it != it_end; ++it)
  {
    Walker_t& thisWalker(**it);
    W.loadWalker(thisWalker,false);
    //RealType nodecorr = setScaledDriftPbyPandNodeCorr(m_tauovermass,W.G,drift);
    RealType nodecorr=setScaledDriftPbyPandNodeCorr(Tau,MassInvP,W.G,drift);
    //create a 3N-Dimensional Gaussian with variance=1
    makeGaussRandomWithEngine(deltaR,RandomGen);
    //if(!W.makeMoveWithDrift(thisWalker,drift,deltaR, m_sqrttau))
    if (!W.makeMoveWithDrift(thisWalker,drift ,deltaR,SqrtTauOverMass))
    {
      H.rejectedMove(W,thisWalker);
      continue;
    }
    //save old local energy
    RealType eold = thisWalker.Properties(LOCALENERGY);
    RealType enew = eold;
    RealType signold = thisWalker.Properties(SIGN);
    RealType logpsi(Psi.evaluateLog(W));
    bool accepted=false;
    RealType rr_accepted = 0.0;
    nodecorr=0.0;
    if(branchEngine->phaseChanged(Psi.getPhaseDiff()))
    {
      thisWalker.Age++;
      thisWalker.willDie();
    }
    else
    {
      enew=H.evaluate(W);
      RealType logGf = -0.5*Dot(deltaR,deltaR);
      //nodecorr = setScaledDriftPbyPandNodeCorr(m_tauovermass,W.G,drift);
      nodecorr=setScaledDriftPbyPandNodeCorr(Tau,MassInvP,W.G,drift);
      deltaR = thisWalker.R - W.R - drift;
      RealType logGb=logBackwardGF(deltaR);
      //RealType logGb = -m_oneover2tau*Dot(deltaR,deltaR);
      RealType prob= std::min(std::exp(logGb-logGf +2.0*(logpsi-thisWalker.Properties(LOGPSI))),1.0);
      //calculate rr_proposed here
      deltaR = W.R-thisWalker.R;
      RealType rr_proposed = Dot(deltaR,deltaR);
      if(RandomGen() > prob)
      {
        enew=eold;
        thisWalker.Age++;
        thisWalker.Properties(R2ACCEPTED)=0.0;
        thisWalker.Properties(R2PROPOSED)=rr_proposed;
        H.rejectedMove(W,thisWalker);
      }
      else
      {
        thisWalker.Age=0;
        accepted=true;
        W.saveWalker(thisWalker);
        rr_accepted = rr_proposed;
        thisWalker.resetProperty(logpsi,Psi.getPhase(),enew,rr_accepted,rr_proposed,nodecorr);
        H.auxHevaluate(W,thisWalker);
        H.saveProperty(thisWalker.getPropertyBase());
      }
//         cout<<logpsi<<"  "<<Psi.getPhase()<<"  "<<enew<<"  "<<rr_accepted<<"  "<<rr_proposed<<"  "<<nodecorr<<endl;
      thisWalker.Weight *= branchEngine->branchWeight(enew,eold);
    }
    if(accepted)
      ++nAccept;
    else
      ++nReject;
  }
}
}

/***************************************************************************
 * $RCSfile: DMCUpdateAll.cpp,v $   $Author: jnkim $
 * $Revision: 5884 $   $Date: 2013-06-10 08:50:16 -0500 (Mon, 10 Jun 2013) $
 * $Id: DMCUpdateAll.cpp 5884 2013-06-10 13:50:16Z jnkim $
 ***************************************************************************/
