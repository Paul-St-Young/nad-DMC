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
    // update determinant coeffients with ion position
    // !!! this is the first implementation, hard-code to do CH
    if (Interpolate){
      RealType CHdistance=std::sqrt( dot(DMCIons.R[1],DMCIons.R[1]) );
      RealType CHdistanceo=2.116493107; // Bohr

      Psi.updateCoeff(CHdistance-CHdistanceo); // update coefficients after ion move
    }
}

/// Constructor.
DMCUpdateAllWithIons::DMCUpdateAllWithIons(MCWalkerConfiguration& w,
					   TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator_t& rg, 
					   ParticleSet& ions, std::vector<int> ionsToMove, bool Restart, bool Interpolate):
  QMCUpdateBase(w,psi,h,rg), DMCIons(ions), Restart(Restart), Interpolate(Interpolate)
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
      //cout << "DMC rc " << jj << " " << etaisyys[jj] << endl;
      ++jj;
    }
  }



}

/// destructor
DMCUpdateAllWithIons::~DMCUpdateAllWithIons() { }



void DMCUpdateAllWithIons::dorotateshift(PosType& origin, PosType &second, RealType &getpsi )
{
  PosType neworigin, newsec,displ,testsec;
  neworigin = origin-origin;
  newsec = origin - origin;  //zero out newsec
  displ = second - origin;
  RealType odist= std::sqrt(displ[0]*displ[0]+displ[1]*displ[1]+displ[2]*displ[2]);
  newsec[2] = odist;  //align sec along z direction 
  //Now we have two new positions of ions neworigin and new sec
/*
  RealType rot1 = std::atan(displ[0]/displ[2]);
  RealType rot2 = std::atan(displ[1]/(std::sin(rot1)*displ[0]+std::cos(rot1)*displ[2]));
  RealType rangle;
  PosType  newpos = displ;
  rangle = rot1;
 cout << "before newpos " << newpos << endl; 
 testsec[0]  = std::cos(rangle)*newpos[0]-std::sin(rangle)*newpos[2];
 testsec[1]  = newpos[1];
 testsec[2]  = std::sin(rangle)*newpos[0]+std::cos(rangle)*newpos[2];
 
 cout << "testsec " << testsec << endl; 
 
  rangle = rot2;
 newpos[0] = testsec[0];
 newpos[1] = testsec[1];
 newpos[2] = testsec[2];

 testsec[0]  = newpos[0];
 testsec[1]  = std::cos(rangle)*newpos[1]-std::sin(rangle)*newpos[2];
 testsec[2]  = std::sin(rangle)*newpos[1]+std::cos(rangle)*newpos[2];

 cout << "finalsec " << testsec << endl; 
 cout << "and compare " << newsec << endl;
*/
 ///Rotate electrons
   W.makeShiftRotate(deltaR,origin,second);
   DMCIons.R[ion_index[0]] = neworigin;
   DMCIons.R[ion_index[1]] = newsec;
   DMCIons.update();
   updateCoeff();
   H.update_source(DMCIons);
   W.update();
  
 ///Calculate wavefunction
   RealType logpsi(Psi.evaluateLog(W));
 //   cout << "my psi middle " <<logpsi << endl;
   getpsi = logpsi;  
  
 ///Rotate back electrons
   W.invmakeShiftRotate(deltaR,origin,second);
 //Rotate gradient vector, noshift
   W.gradRotate(deltaR,origin,second);
 ///Put back initial ion 
   DMCIons.R[ion_index[0]] = origin;
   DMCIons.R[ion_index[1]] = second;
   DMCIons.update();
   updateCoeff();
   H.update_source(DMCIons);
   W.update();
 ///
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
    //RealType nodecorr = setScaledDriftPbyPandNodeCorr(m_tauovermass,W.G,drift);
    makeGaussRandomWithEngine(deltaR,RandomGen);
    //if(!W.makeMoveWithDrift(thisWalker,drift,deltaR, m_sqrttau))

    // defined in QMCUpdateBase.h 
    //(can be given as input for ions also, e.g. <parameter name="ionGrid">0.00001</parameter>)
    RealType GridSpacing=UniformGrid_h*1000; 
    //RealType GridSpacing=0.0001; 

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
    //RealType f2old=nuclei_wfs(ionR,thisWalker.ionPos.size());
    RealType f2old=nuclei_wfs(DMCIons.R);
    RealType wfs[3][2];
    //RealType wfs2[3][2];
    //RealType f=std::exp(Psi.evaluateLog(W));
    RealType f=std::exp(thisWalker.Properties(LOGPSI));
    RealType h=GridSpacing;
    
    PosType ionsDrift[thisWalker.ionPos.size()];
    
    for (int j=0;j<thisWalker.ionPos.size();++j)
      {

	//if(tauCountFreq[ion_index[j]]>1 && tauCount%tauCountFreq[ion_index[j]]) continue;

	//RealType tau_eff=(RealType)tauCountFreq[ion_index[j]]*Tau;

	RealType tau_eff=Tau;
        //// Norm First update grad and lap
	for (int i=0;i<3;++i)
	  {
	    DMCIons.R[ion_index[j]][i] = ionR[j][i]-h;
	    DMCIons.update();
   updateCoeff();
	    W.update();	
            np1 = DMCIons.R[ion_index[0]]; np2 = DMCIons.R[ion_index[1]];
            dorotateshift(np1,np2,wfs[i][0]);
            wfs[i][0] = std ::exp(wfs[i][0]);  
	//    wfs[i][0]=std::exp(Psi.evaluateLog(W));
	    //wfs2[i][0]=nuclei_wfs(DMCIons.R);
	    
	    DMCIons.R[ion_index[j]][i] = ionR[j][i]+h;
	    DMCIons.update();
   updateCoeff();
	    W.update();	  
            np1 = DMCIons.R[ion_index[0]]; np2 = DMCIons.R[ion_index[1]];
            dorotateshift(np1,np2,wfs[i][1]);
            wfs[i][1] = std ::exp(wfs[i][1]);  
	//    wfs[i][1]=std::exp(Psi.evaluateLog(W));
	    //wfs2[i][1]=nuclei_wfs(DMCIons.R);
	    
	    DMCIons.R[ion_index[j]][i] = ionR[j][i];
	    DMCIons.update();
   updateCoeff();
	    W.update();	  	  
	  }
	
	//ionsDrift[j] = ionDrift(wfs,h,f,Tau/DMCIons.Mass[j]/2.0)+ionDrift(wfs2,h,f2old,Tau/DMCIons.Mass[j]/2.0);
	
	ionsDrift[j] = ionDrift(wfs,h,f,tau_eff/DMCIons.Mass[ion_index[j]]/2.0)+nuclei_wfs_gradient(DMCIons.R,ion_index[j])*tau_eff/DMCIons.Mass[ion_index[j]];

	//cout << "drift 1 " << j << " " << ionDrift(wfs,h,f,Tau/DMCIons.Mass[j]/2.0) << endl;
	//cout << "drift 2 " << j << " " << ionDrift(wfs2,h,f2old,Tau/DMCIons.Mass[j]/2.0) << endl;
      }

    //get wfs and update laplacians and gradients
    //commented out by norm  //RealType logpsi(Psi.evaluateLog(W));
    np1 = DMCIons.R[ion_index[0]]; np2 = DMCIons.R[ion_index[1]];
    dorotateshift(np1,np2,logpsi);
     
    // move ions
    RealType ilogGf=0.0;
    for (int i=0;i<thisWalker.ionPos.size();++i)
      {
	//if(tauCountFreq[ion_index[i]]>1 && tauCount%tauCountFreq[ion_index[i]]) continue;

	//RealType tau_eff=(RealType)tauCountFreq[ion_index[i]]*Tau;

	
	RealType tau_eff=Tau;


	ilogGf += -0.5*dot(deltaR[i],deltaR[i]);

        // Norm moves ions, but nothing needs to be done as long as I put
        // everything back above	
	ionR[i]=ionR[i]+std::sqrt(Tau/DMCIons.Mass[ion_index[i]])*deltaR[i] + ionsDrift[i];
	//DMCIons.R[i] = ionR[i]; 
      }

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
      
      
      //W.loadWalker(thisWalker,true);
      for (int i=0;i<thisWalker.ionPos.size();++i)
	{
	  thisWalker.ionPos[i]=ionRorig[i];
	  W.ionPos[i]=ionRorig[i];
	}
      
    //  W.saveWalker(thisWalker);
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

    //evaluate wave function
    //RealType f2 = nuclei_wfs(ionR,thisWalker.ionPos.size());
    RealType f2 = nuclei_wfs(DMCIons.R);
    RealType g2 = f2*f2/f2old/f2old;
    //RealType logpsi(Psi.evaluateLog(W));//get wfs and update laplacians and gradients
    //commented out by norm //  logpsi=Psi.evaluateLog(W);//get wfs and update laplacians and gradients
    np1 = DMCIons.R[ion_index[0]]; np2 = DMCIons.R[ion_index[1]];
    dorotateshift(np1,np2,logpsi);


    if(UseTMove)
      nonLocalOps.reset();
    bool accepted=false;
    RealType rr_accepted = 0.0;
    nodecorr=0.0;
//    if(branchEngine->phaseChanged(Psi.getPhaseDiff()))
     if(0)
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
      
      H.rejectedMove(W,thisWalker);
      
      //W.loadWalker(thisWalker,true);
      for (int i=0;i<thisWalker.ionPos.size();++i)
	{
	  thisWalker.ionPos[i]=ionRorig[i];
	  W.ionPos[i]=ionRorig[i];
	}
                
   //   W.saveWalker(thisWalker);  
    }
    else
    {
      if(UseTMove)
        enew=H.evaluate(W,nonLocalOps.Txy);
      else
	{	  
	  RealType wfs[3][2];
	  //RealType wfs2[3][2];
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
//		  wfs[i][0]=std::exp(Psi.evaluateLog(W));
		  //wfs2[i][0]=nuclei_wfs(DMCIons.R);
		  
		  DMCIons.R[ion_index[j]][i] = ionR[j][i]+h;
		  DMCIons.update();
   updateCoeff();
		  W.update();	  
                  np1 = DMCIons.R[ion_index[0]]; np2 = DMCIons.R[ion_index[1]];
                  dorotateshift(np1,np2,wfs[i][1]);
                  wfs[i][1] = std ::exp(wfs[i][1]);  
	//	  wfs[i][1]=std::exp(Psi.evaluateLog(W));
		  //wfs2[i][1]=nuclei_wfs(DMCIons.R);
		  
		  DMCIons.R[ion_index[j]][i] = ionR[j][i];
		  DMCIons.update();
   updateCoeff();
		  W.update();	  	  
		}
	      //ionsKineticE += IonKineticEnergy(wfs,f,h,DMCIons.Mass[j]);
	      //ionsKineticE += IonKineticEnergy3(wfs,wfs2,f,f2,h,DMCIons.Mass[j]);

	      ionsKineticE += IonKineticEnergy3_2(wfs,f,h,DMCIons.Mass[ion_index[j]],ion_index[j],DMCIons.R);
	      //should also give the gradient as output

	      //ionsDrift[j] = ionDrift(wfs,h,f,Tau/DMCIons.Mass[j]/2.0)+ionDrift(wfs2,h,f2,Tau/DMCIons.Mass[j]/2.0);

	      //if(tauCountFreq[ion_index[j]]>1 && tauCount%tauCountFreq[ion_index[j]]) continue;

	      //RealType tau_eff=(RealType)tauCountFreq[ion_index[j]]*Tau;
	      

	      RealType tau_eff=Tau;

	      ionsDrift[j] = ionDrift(wfs,h,f,tau_eff/DMCIons.Mass[ion_index[j]]/2.0)+nuclei_wfs_gradient(DMCIons.R,ion_index[j])*tau_eff/DMCIons.Mass[ion_index[j]];
	      
	      /*
		cout << "drift fd " << ionDrift(wfs2,h,f2,Tau/DMCIons.Mass[j]/2.0) << endl;
		cout << "drift an " << nuclei_wfs_gradient(DMCIons.R,j)*Tau/DMCIons.Mass[j] << endl;
		cout << "E fd " << IonKineticEnergy3(wfs,wfs2,f,f2,h,DMCIons.Mass[j]) << endl;
		cout << "E an " << IonKineticEnergy3_2(wfs,f,h,DMCIons.Mass[j],j,DMCIons.R) << endl;
	      */
	    }
	  

	  //(get wfs and) update laplacians and gradients
	//normcommented out  logpsi = Psi.evaluateLog(W);
          np1 = DMCIons.R[ion_index[0]]; np2 = DMCIons.R[ion_index[1]];
          dorotateshift(np1,np2,logpsi);
	  enew=H.evaluate(W)+ionsKineticE;
	  
	  
	  nucleiKineticE2 = ionsKineticE;
	  
	}
      
      RealType logGf = -0.5*Dot(deltaR,deltaR)+ilogGf;
      //RealType nodecorr = setScaledDriftPbyPandNodeCorr(m_tauovermass,W.G,drift);
      RealType nodecorr=setScaledDriftPbyPandNodeCorr(Tau,MassInvP,W.G,drift);
      deltaR = thisWalker.R - W.R - drift;
      PosType r_temp[thisWalker.ionPos.size()];
      RealType ilogGb=0.0;
      for (int i=0;i<thisWalker.ionPos.size();++i)
	{

	  //if(tauCountFreq[ion_index[i]]>1 && tauCount%tauCountFreq[ion_index[i]]) continue;	  
	  //RealType tau_eff=(RealType)tauCountFreq[ion_index[i]]*Tau;	  
	  

	  RealType tau_eff=Tau;

	  for (int j=0;j<3;++j)
	    r_temp[i][j]=ionRorig[i][j]-ionR[i][j]-ionsDrift[i][j];
	  ilogGb += -dot(r_temp[i],r_temp[i])/2.0/tau_eff*DMCIons.Mass[ion_index[i]];
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
	
//	W.saveWalker(thisWalker);
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
//testing   testing testing
 // cout << "new position 1 " <<  ionR[0] << endl;
 // cout << "new position 2 " <<  ionR[1] << endl;
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
 // cout << "rotate position 1 " <<  ntneworigin << endl;
//  cout << "rotate position 2 " <<  ntnewsec << endl;
//  cout << "enew " <<  enew <<" eval  " << H.evaluate(W) << " ik " << nucleiKineticE  <<endl;
//end testing   testing testing
 
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
        logpsi=Psi.evaluateLog(W);  //This is in tmove, ignore
        setScaledDrift(Tau,W.G,drift);
        thisWalker.resetProperty(logpsi,Psi.getPhase(),eold);
        thisWalker.R[iat] = W.R[iat];
        ++NonLocalMoveAccepted;
      }
    }
    thisWalker.Weight *= branchEngine->branchWeight(enew,eold);
   // thisWalker.Weight = 1;
//	  thisWalker.Age=0;
    //branchEngine->accumulate(eold,1);
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
