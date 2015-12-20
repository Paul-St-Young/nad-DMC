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

void DMCUpdateAllWithIons::updateCoeff(){
    // update determinant coeffients with ion separation
    if (Interpolate){
      RealType CHdistance=std::sqrt( dot(DMCIons.R[1],DMCIons.R[1]) );

      Psi.updateCoeff(CHdistance); // update coefficients after ion move
    }
}

/// Constructor.
DMCUpdateAllWithIons::DMCUpdateAllWithIons(MCWalkerConfiguration& w,
					   TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator_t& rg, 
					   ParticleSet& ions, std::vector<int> ionsToMove, bool Restart, bool Interpolate):
  QMCUpdateBase(w,psi,h,rg), DMCIons(ions), Restart(Restart), Interpolate(Interpolate)
{ 

  UpdatePbyP=false;

  // recording distances between ions that have mass
  for (int i=0;i<ionsToMove.size();++i)
    ion_index.push_back(ionsToMove[i]);
  int Ndist=DMCIons.R.size();
  Ndist = Ndist*(Ndist-1)/2;
  RealType etaisyys[Ndist];
  for (int i=0;i<Ndist;++i){
    etaisyys[i]=0.0;
    nuclei_dist_step.push_back(0.0);
  }
  int jj=0;
  for (int j=0;j<DMCIons.R.size();++j) {
    for (int k=j+1;k<DMCIons.R.size();++k) {
      RealType etaisyys1=0.0;  
      for (int i=0;i<3;++i)
	      etaisyys1 += (DMCIons.R[j][i]-DMCIons.R[k][i])*(DMCIons.R[j][i]-DMCIons.R[k][i]);
      
      etaisyys[jj] = std::sqrt(etaisyys1);
      ion_rc.push_back(etaisyys[jj]);
      ++jj;
    }
  }

}

/// destructor
DMCUpdateAllWithIons::~DMCUpdateAllWithIons() { }

void DMCUpdateAllWithIons::dorotateshift(PosType& origin, PosType &second, RealType &getpsi )
{
  //Calculate wavefunction
  RealType logpsi(Psi.evaluateLog(W));
  getpsi = logpsi;  
 return ;
}


void DMCUpdateAllWithIons::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end,
    bool measure)
{


  PosType np1,np2;
  RealType logpsi;
  int Ndist=DMCIons.R.size();
  Ndist = Ndist*(Ndist-1)/2;

  RealType nucleiKineticE=0.0;
  RealType nucleiKineticE2=0.0;
  int kineticCount=0;

  RealType etaisyys[Ndist];
  for (int i=0;i<Ndist;++i)
    etaisyys[i]=0.0;
  int et_count=0;

  for(; it != it_end; ++it)
  {
    Walker_t& thisWalker(**it);

    PosType ionR[thisWalker.ionPos.size()];
    PosType ionRorig[thisWalker.ionPos.size()];
    for (int i=0;i<thisWalker.ionPos.size();++i)
    {
      ionR[i]=thisWalker.ionPos[i];
      ionRorig[i]=ionR[i];
      DMCIons.R[ion_index[i]] = ionR[i]; 	
    }    
    DMCIons.update();
    updateCoeff();

    W.loadWalker(thisWalker,true);

    //create a 3N-Dimensional Gaussian with variance=1
    RealType nodecorr=setScaledDriftPbyPandNodeCorr(Tau,MassInvP,W.G,drift);
    makeGaussRandomWithEngine(deltaR,RandomGen);

    // defined in QMCUpdateBase.h 
    //(can be given as input for ions also, e.g. <parameter name="ionGrid">0.00001</parameter>)
    RealType GridSpacing=UniformGrid_h*1000; 

    // update ion distances
    int jj=0;
    for (int j=0;j<DMCIons.R.size();++j) {
      for (int k=j+1;k<DMCIons.R.size();++k) {
        RealType etaisyys1=0.0;
	      for (int i=0;i<3;++i)
          etaisyys1 += (DMCIons.R[j][i]-DMCIons.R[k][i])*(DMCIons.R[j][i]-DMCIons.R[k][i]);

        etaisyys[jj] += std::sqrt(etaisyys1);
        ++jj;
      }
    }
    ++et_count;    


    // calculate drift
    RealType f2old=nuclei_wfs(DMCIons.R);
    RealType wfs[3][2];
    RealType f=std::exp(thisWalker.Properties(LOGPSI));
    RealType h=GridSpacing;
    
    PosType ionsDrift[thisWalker.ionPos.size()];
    
    for (int j=0;j<thisWalker.ionPos.size();++j)
    { 
      RealType tau_eff=Tau;
      // Norm First update grad and lap
	    for (int i=0;i<3;++i)
	    {
        DMCIons.R[ion_index[j]][i] = ionR[j][i]-h;
        DMCIons.update();
        updateCoeff();
        W.update();	
        np1 = DMCIons.R[ion_index[0]]; np2 = DMCIons.R[ion_index[1]];
        dorotateshift(np1,np2,wfs[i][0]);
        wfs[i][0] = std ::exp(wfs[i][0]);  
        
        DMCIons.R[ion_index[j]][i] = ionR[j][i]+h;
        DMCIons.update();
        updateCoeff();
        W.update();	  

        np1 = DMCIons.R[ion_index[0]]; np2 = DMCIons.R[ion_index[1]];
        dorotateshift(np1,np2,wfs[i][1]);
        wfs[i][1] = std ::exp(wfs[i][1]);  

        DMCIons.R[ion_index[j]][i] = ionR[j][i];
        DMCIons.update();
        updateCoeff();
        W.update();	  	  
	    }
	
	    ionsDrift[j] = ionDrift(wfs,h,f,tau_eff/DMCIons.Mass[ion_index[j]]/2.0)+
         nuclei_wfs_gradient(DMCIons.R,ion_index[j])*tau_eff/DMCIons.Mass[ion_index[j]];
    }

    //get wfs and update laplacians and gradients
    //commented out by norm  //RealType logpsi(Psi.evaluateLog(W));
    np1 = DMCIons.R[ion_index[0]]; np2 = DMCIons.R[ion_index[1]];
    dorotateshift(np1,np2,logpsi);
    
    // move ions
    RealType ilogGf=0.0;
    for (int i=0;i<thisWalker.ionPos.size();++i)
    { if (i==ION0) {
      RealType tau_eff=Tau;
      ilogGf += -0.5*dot(deltaR[i],deltaR[i]);
      ionR[i]=ionR[i]+std::sqrt(Tau/DMCIons.Mass[ion_index[i]])*deltaR[i] + ionsDrift[i];
    } }

    makeGaussRandomWithEngine(deltaR,RandomGen);


    if (!W.makeMoveWithDrift(thisWalker,drift ,deltaR,SqrtTauOverMass))
    {
       
      for (int i=0;i<thisWalker.ionPos.size();++i)
      {
        DMCIons.R[ion_index[i]] = thisWalker.ionPos[i];
      }
      DMCIons.update();
      updateCoeff();
      H.update_source(DMCIons);
      W.update();      
      
      H.rejectedMove(W,thisWalker);
      
      for (int i=0;i<thisWalker.ionPos.size();++i)
      {
        thisWalker.ionPos[i]=ionRorig[i];
        W.ionPos[i]=ionRorig[i];
      }
       
      continue;
    }

    //save old local energy
    RealType eold    = thisWalker.Properties(LOCALENERGY);
    RealType signold = thisWalker.Properties(SIGN);
    RealType enew  = eold;

    for (int i=0;i<thisWalker.ionPos.size();++i)
    {
	    DMCIons.R[ion_index[i]] = ionR[i]; 
    }    
    
    DMCIons.update(); // update ions
    updateCoeff();
    H.update_source(DMCIons); // update ions in Hamiltonian 
    W.update();  // update "wave function" i.e. distances

    RealType f2 = nuclei_wfs(DMCIons.R);
    RealType g2 = f2*f2/f2old/f2old;
    np1 = DMCIons.R[ion_index[0]]; np2 = DMCIons.R[ion_index[1]];
    dorotateshift(np1,np2,logpsi);

    if(UseTMove)
      nonLocalOps.reset();
    bool accepted=false;
    RealType rr_accepted = 0.0;
    nodecorr=0.0;
    //if(branchEngine->phaseChanged(Psi.getPhaseDiff()))
    //else
    if(UseTMove)
      enew=H.evaluate(W,nonLocalOps.Txy);
    else
    {	  
      RealType wfs[3][2];
      RealType f=std::exp(logpsi);
      RealType h=GridSpacing;
      
      RealType ionsKineticE=0.0;
      
      for (int j=0;j<thisWalker.ionPos.size();++j)
      {
        
        for (int i=0;i<3;++i)
        {
          DMCIons.R[ion_index[j]][i] = ionR[j][i]-h;
          DMCIons.update();
          updateCoeff();
          W.update();	  
          np1 = DMCIons.R[ion_index[0]]; np2 = DMCIons.R[ion_index[1]];
          dorotateshift(np1,np2,wfs[i][0]);
          wfs[i][0] = std ::exp(wfs[i][0]);  

          DMCIons.R[ion_index[j]][i] = ionR[j][i]+h;
          DMCIons.update();
          updateCoeff();
          W.update();	  
          np1 = DMCIons.R[ion_index[0]]; np2 = DMCIons.R[ion_index[1]];
          dorotateshift(np1,np2,wfs[i][1]);
          wfs[i][1] = std ::exp(wfs[i][1]);  

          DMCIons.R[ion_index[j]][i] = ionR[j][i];
          DMCIons.update();
          updateCoeff();
          W.update();	  	  
        }
        //ionsKineticE += IonKineticEnergy(wfs,f,h,DMCIons.Mass[j]);
        //ionsKineticE += IonKineticEnergy3(wfs,wfs2,f,f2,h,DMCIons.Mass[j]);

        ionsKineticE += IonKineticEnergy3_2(wfs,f,h,
            DMCIons.Mass[ion_index[j]],ion_index[j],DMCIons.R);
        RealType tau_eff=Tau;

        if (j==ION0) ionsDrift[j] = ionDrift(wfs,h,f,tau_eff/DMCIons.Mass[ion_index[j]]/2.0)
          +nuclei_wfs_gradient(DMCIons.R,ion_index[j])*tau_eff/DMCIons.Mass[ion_index[j]];
        
      }

      //(get wfs and) update laplacians and gradients
      //normcommented out  logpsi = Psi.evaluateLog(W);
      np1 = DMCIons.R[ion_index[0]]; np2 = DMCIons.R[ion_index[1]];
      dorotateshift(np1,np2,logpsi);
      enew=H.evaluate(W)+ionsKineticE;
      nucleiKineticE2 = ionsKineticE;
    }
      
    // forward Green function exponent, (ilogGf if Green function exponent for ions)
    RealType logGf = -0.5*Dot(deltaR,deltaR)+ilogGf;
    nodecorr=setScaledDriftPbyPandNodeCorr(Tau,MassInvP,W.G,drift);
    deltaR = thisWalker.R - W.R - drift;

    // backward Green function exponent
    RealType ilogGb=0.0;
    PosType r_temp[thisWalker.ionPos.size()];
    for (int i=0;i<thisWalker.ionPos.size();++i)
    { if (i==ION0) {
      RealType tau_eff=Tau;
      for (int j=0;j<3;++j)
        r_temp[i][j]=ionRorig[i][j]-ionR[i][j]-ionsDrift[i][j];
      ilogGb += -dot(r_temp[i],r_temp[i])/2.0/tau_eff*DMCIons.Mass[ion_index[i]];
    } }
    RealType logGb=logBackwardGF(deltaR)+ilogGb;
    RealType prob= std::min(g2*std::exp(logGb-logGf +2.0*(logpsi-thisWalker.Properties(LOGPSI))),1.0);
    deltaR = W.R-thisWalker.R;
    RealType rr_proposed = Dot(deltaR,deltaR);
    if(RandomGen() > prob)
    {
      thisWalker.Age++;	

      for (int i=0;i<thisWalker.ionPos.size();++i)
      {
        DMCIons.R[ion_index[i]] = thisWalker.ionPos[i];
      }
      DMCIons.update();
      updateCoeff();
      H.update_source(DMCIons);
      W.update();
    
      enew=eold;
      thisWalker.Properties(R2ACCEPTED)=0.0;
      thisWalker.Properties(R2PROPOSED)=rr_proposed;
      H.rejectedMove(W,thisWalker);
    
      //W.loadWalker(thisWalker,true);
      for (int i=0;i<thisWalker.ionPos.size();++i)
      {
        thisWalker.ionPos[i]=ionRorig[i];
        W.ionPos[i]=ionRorig[i];
      }      
    }
    else
	  {
      nucleiKineticE += nucleiKineticE2;
      ++kineticCount;

      accepted=true;
      for (int i=0;i<thisWalker.ionPos.size();++i)
      {
        W.ionPos[i]=ionR[i];
        thisWalker.ionPos[i]=ionR[i];
      }	 
      PosType ntneworigin, ntnewsec,ntdispl,nttestsec;
      ntneworigin = ionR[0] - ionR[0];
      ntnewsec = ionR[0] - ionR[0];
      ntdispl = ionR[1] - ionR[0];
      RealType ntodist= std::sqrt(ntdispl[0]*ntdispl[0]+ntdispl[1]*ntdispl[1]+ntdispl[2]*ntdispl[2]);
      ntnewsec[2] = ntodist;  //align sec along z direction   
      W.ionPos[0]=ntneworigin;
      W.ionPos[1]=ntnewsec;
      thisWalker.ionPos[0]= ntneworigin;
      thisWalker.ionPos[1]= ntnewsec;
      DMCIons.R[ion_index[0]] = ntneworigin;
      DMCIons.R[ion_index[1]] = ntnewsec;
      DMCIons.update();
      updateCoeff();
      H.update_source(DMCIons);
      W.update();
      W.makeShiftRotate(deltaR,ionR[0],ionR[1]);
      dorotateshift(ntneworigin,ntnewsec,logpsi);

      thisWalker.Age=0;
      W.saveWalker(thisWalker);
      rr_accepted = rr_proposed;
      thisWalker.resetProperty(logpsi,Psi.getPhase(),enew,rr_accepted,rr_proposed,nodecorr);
      H.auxHevaluate(W,thisWalker);
      H.saveProperty(thisWalker.getPropertyBase());
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
        logpsi=Psi.evaluateLog(W);  //This is in tmove, ignore
        setScaledDrift(Tau,W.G,drift);
        thisWalker.resetProperty(logpsi,Psi.getPhase(),eold);
        thisWalker.R[iat] = W.R[iat];
        ++NonLocalMoveAccepted;
      }
    }
    thisWalker.Weight *= branchEngine->branchWeight(enew,eold);
    
    /*// always test that DMC with no weights = VMC
    thisWalker.Weight = 1;
    thisWalker.Age=0;
    branchEngine->accumulate(eold,1);
    */

    if(accepted)
      ++nAccept;
    else
      ++nReject;
  }
/*
  std::ofstream myfile;
  myfile.open("nucleiKineticE.txt", ios::out | ios::app);
  myfile << nucleiKineticE/kineticCount << endl;
  myfile.close();
*/
  for (int i=0;i<Ndist;++i)
    nuclei_dist_step[i] += etaisyys[i]/et_count;

}

}

/***************************************************************************
 * $RCSfile: DMCUpdateAll.cpp,v $   $Author: jnkim $
 * $Revision: 5884 $   $Date: 2013-06-10 08:50:16 -0500 (Mon, 10 Jun 2013) $
 * $Id: DMCUpdateAll.cpp 5884 2013-06-10 13:50:16Z jnkim $
 ***************************************************************************/
