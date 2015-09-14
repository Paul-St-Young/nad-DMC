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
DMCUpdateAllWithIons::DMCUpdateAllWithIons(MCWalkerConfiguration& w,
					   TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator_t& rg, ParticleSet& ions,std::vector<int> ionsToMove):
  QMCUpdateBase(w,psi,h,rg), DMCIons(ions)
{ 

  UpdatePbyP=false;

  /*

  MCWalkerConfiguration::iterator
    wit(W.begin()), wit_end(W.end());
  
  for(; wit != wit_end; ++wit)
    {
      Walker_t& thisWalker(**wit);
      Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
      W.loadWalker(thisWalker,true);

      for (int i=0;i<DMCIons.R.size();++i)
        {
          // Mass==0 is default for fixed nuclei (not used)                                     
          // i.e. if defined Mass>1 then those ions are moved                                   
	  if (DMCIons.Mass[i]>0.0)
	    {	      
	      //W.ionPos[i]=DMCIons.R[i];
	      //thisWalker.ionPos[i]=DMCIons.R[i];
	    }
        }
      W.saveWalker(thisWalker);
    }
  */

}

/// destructor
DMCUpdateAllWithIons::~DMCUpdateAllWithIons() { }

void DMCUpdateAllWithIons::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end,
    bool measure)
{
  for(; it != it_end; ++it)
  {
    Walker_t& thisWalker(**it);
    W.loadWalker(thisWalker,true);


    //create a 3N-Dimensional Gaussian with variance=1
    RealType nodecorr=setScaledDriftPbyPandNodeCorr(Tau,MassInvP,W.G,drift);
    //RealType nodecorr = setScaledDriftPbyPandNodeCorr(m_tauovermass,W.G,drift);
    makeGaussRandomWithEngine(deltaR,RandomGen);
    //if(!W.makeMoveWithDrift(thisWalker,drift,deltaR, m_sqrttau))

    RealType GridSpacing=0.000001;

    // move ions

    PosType ionR[thisWalker.ionPos.size()];
    PosType ionRorig[thisWalker.ionPos.size()];
    for (int i=0;i<thisWalker.ionPos.size();++i)
      {
	ionR[i]=thisWalker.ionPos[i];
	ionRorig[i]=ionR[i];
	DMCIons.R[i] = ionR[i]; 
      }

    RealType coeff = 10.0;
    RealType rc = 2.68;

    RealType f2old=testGaussian(ionR[0],ionR[1],coeff,rc);

    RealType wfs[3][2];
    RealType wfs2[3][2];
    RealType f=std::exp(Psi.evaluateLog(W));
    RealType h=GridSpacing;
    
    PosType ionsDrift[thisWalker.ionPos.size()];
    
    for (int j=0;j<thisWalker.ionPos.size();++j)
      {
	for (int i=0;i<3;++i)
	  {
	    DMCIons.R[j][i] = ionR[j][i]-h;
	    DMCIons.update();
	    W.update();	  
	    wfs[i][0]=std::exp(Psi.evaluateLog(W));
	    wfs2[i][0]=testGaussian(DMCIons.R[0],DMCIons.R[1],coeff,rc);
	    
	    DMCIons.R[j][i] = ionR[j][i]+h;
	    DMCIons.update();
	    W.update();	  
	    wfs[i][1]=std::exp(Psi.evaluateLog(W));
	    wfs2[i][1]=testGaussian(DMCIons.R[0],DMCIons.R[1],coeff,rc);
	    
	    DMCIons.R[j][i] = ionR[j][i];
	    DMCIons.update();
	    W.update();	  	  
	  }
	
	ionsDrift[j] = ionDrift(wfs,h,f,Tau/DMCIons.Mass[j]/2.0)+ionDrift(wfs2,h,f2old,Tau/DMCIons.Mass[j]/2.0);
      }
    
    RealType ilogGf=0.0;
    for (int i=0;i<thisWalker.ionPos.size();++i)
      {
	ilogGf += -0.5*dot(deltaR[i],deltaR[i]);
	ionR[i]=ionR[i]+std::sqrt(Tau/DMCIons.Mass[i])*deltaR[i] + ionsDrift[i];
	DMCIons.R[i] = ionR[i]; 
      }
    

    makeGaussRandomWithEngine(deltaR,RandomGen);


    if (!W.makeMoveWithDrift(thisWalker,drift ,deltaR,SqrtTauOverMass))
    {
      
      for (int i=0;i<thisWalker.ionPos.size();++i)
	{
	  DMCIons.R[i] = thisWalker.ionPos[i];
	}
      DMCIons.update();
      H.update_source(DMCIons);
      W.update();      
      
      H.rejectedMove(W,thisWalker);
      continue;
    }
        
    
    //save old local energy
    RealType eold    = thisWalker.Properties(LOCALENERGY);
    RealType signold = thisWalker.Properties(SIGN);
    RealType enew  = eold;
    
    //makeGaussRandomWithEngine(deltaR,RandomGen);

    
    DMCIons.update(); // update ions
    H.update_source(DMCIons); // update ions in Hamiltonian 
    W.update();  // update "wave function" i.e. distances

    //evaluate wave function
    RealType f2=testGaussian(ionR[0],ionR[1],coeff,rc);
    RealType g2 = f2*f2/f2old/f2old;
    RealType logpsi(Psi.evaluateLog(W));


    if(UseTMove)
      nonLocalOps.reset();
    bool accepted=false;
    RealType rr_accepted = 0.0;
    nodecorr=0.0;
    if(branchEngine->phaseChanged(Psi.getPhaseDiff()))
    {
      thisWalker.Age++;
      
      for (int i=0;i<thisWalker.ionPos.size();++i)
	{
	  DMCIons.R[i] = thisWalker.ionPos[i];
	}
      DMCIons.update();
      H.update_source(DMCIons);
      W.update();
      
      H.rejectedMove(W,thisWalker);      
    }
    else
    {
      if(UseTMove)
        enew=H.evaluate(W,nonLocalOps.Txy);
      else
	{	  
	  RealType wfs[3][2];
	  RealType wfs2[3][2];
	  RealType f=std::exp(logpsi);
	  RealType h=GridSpacing;
	  
	  RealType ionsKineticE=0.0;
	  
	  for (int j=0;j<thisWalker.ionPos.size();++j)
	    {
	      for (int i=0;i<3;++i)
		{
		  DMCIons.R[j][i] = ionR[j][i]-h;
		  DMCIons.update();
		  W.update();	  
		  wfs[i][0]=std::exp(Psi.evaluateLog(W));
		  wfs2[i][0]=testGaussian(DMCIons.R[0],DMCIons.R[1],coeff,rc);
		  
		  DMCIons.R[j][i] = ionR[j][i]+h;
		  DMCIons.update();
		  W.update();	  
		  wfs[i][1]=std::exp(Psi.evaluateLog(W));
		  wfs2[i][1]=testGaussian(DMCIons.R[0],DMCIons.R[1],coeff,rc);
		  
		  DMCIons.R[j][i] = ionR[j][i];
		  DMCIons.update();
		  W.update();	  	  
		}
	      //ionsKineticE += IonKineticEnergy(wfs,f,h,DMCIons.Mass[j]);
	      ionsKineticE += IonKineticEnergy3(wfs,wfs2,f,f2,h,DMCIons.Mass[j]);
	      ionsDrift[j] = ionDrift(wfs,h,f,Tau/DMCIons.Mass[j]/2.0)+ionDrift(wfs2,h,f2,Tau/DMCIons.Mass[j]/2.0);
	    }
	  
	  
	  enew=H.evaluate(W)+ionsKineticE;
	  
	}
      
      RealType logGf = -0.5*Dot(deltaR,deltaR)+ilogGf;
      //RealType nodecorr = setScaledDriftPbyPandNodeCorr(m_tauovermass,W.G,drift);
      RealType nodecorr=setScaledDriftPbyPandNodeCorr(Tau,MassInvP,W.G,drift);
      deltaR = thisWalker.R - W.R - drift;
      PosType r_temp[thisWalker.ionPos.size()];
      RealType ilogGb=0.0;
      for (int i=0;i<thisWalker.ionPos.size();++i)
	{
	  for (int j=0;j<3;++j)
	    r_temp[i][j]=ionRorig[i][j]-ionR[i][j]-ionsDrift[i][j];
	  ilogGb += -dot(r_temp[i],r_temp[i])/2.0/Tau*DMCIons.Mass[i];
	}
      RealType logGb=logBackwardGF(deltaR)+ilogGb;
      //RealType logGb = -m_oneover2tau*Dot(deltaR,deltaR);
      RealType prob= std::min(g2*std::exp(logGb-logGf +2.0*(logpsi-thisWalker.Properties(LOGPSI))),1.0);
      //calculate rr_proposed here
      deltaR = W.R-thisWalker.R;
      RealType rr_proposed = Dot(deltaR,deltaR);
      if(RandomGen() > prob)
      {
        thisWalker.Age++;	

	for (int i=0;i<thisWalker.ionPos.size();++i)
	  {
	    DMCIons.R[i] = thisWalker.ionPos[i];
	  }
	DMCIons.update();
	H.update_source(DMCIons);
	W.update();
	
	
        enew=eold;
        thisWalker.Properties(R2ACCEPTED)=0.0;
        thisWalker.Properties(R2PROPOSED)=rr_proposed;
        H.rejectedMove(W,thisWalker);
      }
      else
	{
	  accepted=true;
	  for (int i=0;i<thisWalker.ionPos.size();++i)
	    {
	      W.ionPos[i]=ionR[i];
	      thisWalker.ionPos[i]=ionR[i];
	    }	  
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

}

/***************************************************************************
 * $RCSfile: DMCUpdateAll.cpp,v $   $Author: jnkim $
 * $Revision: 5884 $   $Date: 2013-06-10 08:50:16 -0500 (Mon, 10 Jun 2013) $
 * $Id: DMCUpdateAll.cpp 5884 2013-06-10 13:50:16Z jnkim $
 ***************************************************************************/
