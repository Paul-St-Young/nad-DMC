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
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCDrivers/VMC/VMCUpdateAll.h"
#include "QMCDrivers/DriftOperators.h"
#include "Message/OpenMP.h"

namespace qmcplusplus
{

VMCUpdateAllWithIons::VMCUpdateAllWithIons(MCWalkerConfiguration& w, TrialWaveFunction& psi,
					   QMCHamiltonian& h, RandomGenerator_t& rg, ParticleSet& ions,
					   std::vector<int> ionsToMove, bool Restart, bool Interpolate)
  : QMCUpdateBase(w,psi,h,rg), VMCIons(ions), Restart(Restart), Interpolate(Interpolate)
// "ions" is initialized with coordinates read from the ptcl.xml
{
  UpdatePbyP=false;

  MCWalkerConfiguration::iterator
    wit(W.begin()), wit_end(W.end());
  app_log() << "initializing QMCDriver: VMCUpdateAllWithIons" << endl; 

  // the ions to be moved (i.e. mass>0) need to be introduced first the ions ParticleSet
  // at the moment
  for(; wit != wit_end; ++wit)
    {
      Walker_t& thisWalker(**wit);
      Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
      W.loadWalker(thisWalker,true);
      for (int i=0;i<VMCIons.R.size();++i)
	{
	  if (VMCIons.Mass[i]>0.0) 
	    {
            if ( Restart ) // this implies W.ionPos is filled by restart
                VMCIons.R[i]=W.ionPos[i];
            else  // initialize ion position in walker using info from inputfile
                W.ionPos[i]=VMCIons.R[i];
            // end if Restart
            thisWalker.ionPos[i]=VMCIons.R[i];
	    }
	}
      W.saveWalker(thisWalker);
    }

  for (int i=0;i<ionsToMove.size();++i)
    ion_index.push_back(ionsToMove[i]);
  
  int Ndist=VMCIons.R.size();
  Ndist = Ndist*(Ndist-1)/2;  
  RealType etaisyys[Ndist];
  for (int i=0;i<Ndist;++i){
    etaisyys[i]=0.0;
    nuclei_dist_step.push_back(0.0);
  }
  int jj=0;
  for (int j=0;j<VMCIons.R.size();++j) {
    for (int k=j+1;k<VMCIons.R.size();++k) {
      RealType etaisyys1=0.0;  
      for (int i=0;i<3;++i)
	etaisyys1 += (VMCIons.R[j][i]-VMCIons.R[k][i])*(VMCIons.R[j][i]-VMCIons.R[k][i]);
      
      etaisyys[jj] = std::sqrt(etaisyys1);
      ion_rc.push_back(etaisyys[jj]);
      //cout << "VMC rc " << jj << " " << etaisyys[jj] << endl;
      ++jj;
    }
  }
}

VMCUpdateAllWithIons::~VMCUpdateAllWithIons()
{
}


void VMCUpdateAllWithIons::updateCoeff(){
    // update determinant coeffients with ion separation
    if (Interpolate){
      RealType CHdistance=std::sqrt( dot(VMCIons.R[1],VMCIons.R[1]) );

      Psi.updateCoeff(CHdistance); // update coefficients after ion move
    }
}


void VMCUpdateAllWithIons::dorotateshift(PosType& origin, PosType &second, RealType &getpsi )
{
  PosType neworigin, newsec,displ,testsec;
  neworigin = origin-origin;
  newsec = origin - origin;  //zero out newsec
  displ = second - origin;
  RealType odist= std::sqrt(displ[0]*displ[0]+displ[1]*displ[1]+displ[2]*displ[2]);
  newsec[2] = odist;  //align sec along z direction 
  //Now we have two new positions of ions neworigin and new sec

  //Rotate electrons
  W.makeShiftRotate(deltaR,origin,second);
  VMCIons.R[ion_index[0]] = neworigin;
  VMCIons.R[ion_index[1]] = newsec;
  VMCIons.update();
  updateCoeff();
  H.update_source(VMCIons);
  W.update();

  ///Calculate wavefunction
  RealType logpsi(Psi.evaluateLog(W));
  getpsi = logpsi;  

  ///Rotate back electrons
  W.invmakeShiftRotate(deltaR,origin,second);
  ///Put back initial ion 
  VMCIons.R[ion_index[0]] = origin;
  VMCIons.R[ion_index[1]] = second;
  VMCIons.update();
  updateCoeff();
  H.update_source(VMCIons);
  W.update();
  return ;
}

void VMCUpdateAllWithIons::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
{

  int Ndist=VMCIons.R.size();
  PosType  np1,np2;
  Ndist = Ndist*(Ndist-1)/2;
  
  RealType etaisyys[Ndist];
  for (int i=0;i<Ndist;++i)
    etaisyys[i]=0.0;
  int et_count=0;


  for (; it!= it_end; ++it)
  {
    MCWalkerConfiguration::Walker_t& thisWalker(**it);

    Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
    W.loadWalker(thisWalker,true);

    // make random numbers for ion move
    makeGaussRandomWithEngine(deltaR,RandomGen);

    PosType ionR[thisWalker.ionPos.size()];
    for (int i=0;i<thisWalker.ionPos.size();++i)
    {
      ionR[i]=thisWalker.ionPos[i];
      VMCIons.R[ion_index[i]] = ionR[i];
      // cout << ionR[i] << endl; // for debugging restart !!!!!!!!!!!
    }    

    // store distances between ions in etaissys
    int jj=0;
    for (int j=0;j<VMCIons.R.size();++j) {
      for (int k=j+1;k<VMCIons.R.size();++k) {
        RealType etaisyys1=0.0;  
        for (int i=0;i<3;++i)
          etaisyys1 += (VMCIons.R[j][i]-VMCIons.R[k][i])*(VMCIons.R[j][i]-VMCIons.R[k][i]);
        
        etaisyys[jj] += std::sqrt(etaisyys1);
        ++jj;
      }
    }
    ++et_count;        

    RealType f2old=nuclei_wfs(VMCIons.R);

    // make ion move
    for (int i=0;i<thisWalker.ionPos.size();++i)
    {	
      if(tauCountFreq[ion_index[i]]>1 && tauCount%tauCountFreq[ion_index[i]]) continue;	  
      RealType tau_eff=(RealType)tauCountFreq[ion_index[i]]*Tau;	
      //ntdeubg delete me (the factor 2)
      ionR[i]=ionR[i]+std::sqrt(tau_eff/VMCIons.Mass[ion_index[i]])*deltaR[i];

      VMCIons.R[ion_index[i]] = ionR[i];
    }    

    // generate new random numbers for electron move
    makeGaussRandomWithEngine(deltaR,RandomGen);

    if (!W.makeMove(thisWalker,deltaR,SqrtTauOverMass))
    {
      for (int i=0;i<thisWalker.ionPos.size();++i)
      {
	VMCIons.R[ion_index[i]] = thisWalker.ionPos[i];
      }
      VMCIons.update();
   updateCoeff();
      H.update_source(VMCIons);
      W.update();

      H.rejectedMove(W,thisWalker);

      continue; // to next walker
    }

    // why are ion coordinates updated after electron move?
    VMCIons.update();
    updateCoeff();
    H.update_source(VMCIons);
    W.update();    
    // why are ion coordinates updated after electron move?


    RealType f2=nuclei_wfs(VMCIons.R);
    RealType g2 = f2*f2/f2old/f2old;
    RealType logpsi;

    np1 = VMCIons.R[ion_index[0]]; np2 = VMCIons.R[ion_index[1]];
    dorotateshift( np1,np2,logpsi);

    RealType g= std::exp(2.0*(logpsi-thisWalker.Properties(LOGPSI)))*g2;
    if (RandomGen() > g)
    {
      thisWalker.Age++;
      ++nReject;
      for (int i=0;i<thisWalker.ionPos.size();++i)
      {
        VMCIons.R[ion_index[i]] = thisWalker.ionPos[i];
      }
      VMCIons.update();
      updateCoeff();
      H.update_source(VMCIons);
      W.update();
      H.rejectedMove(W,thisWalker);
    }
    else
    { // if move accepted

      // calculate kinetic energy using finite difference 
      //  this block should be inside a function
      RealType wfs[3][2];
      wfs[0][0] = 0;
      wfs[0][1] = 0;
      wfs[1][0] = 0; 
      wfs[1][1] = 0; 
      wfs[2][0] = 0;
      wfs[2][1] = 0;

      RealType f=std::exp(logpsi);
      RealType h=UniformGrid_h*1000; // defined in QMCUpdateBase.h
      
      RealType ionsKineticE=0.0;
      for (int j=0;j<thisWalker.ionPos.size();++j)
	    {	// ionsKinetic 
        for (int i=2;i<2;++i)
        { 
          if(j >1) cout << "This cant be happening, check non adabatic rotate" <<endl;
          //find distances between ions
          int ionidx = 0;
          if(j ==0) ionidx = 1;
          if(j ==1) ionidx = 0;
          RealType zdist = ionR[ionidx][2] - ionR[j][2];
          RealType rangle  = std::atan(h/zdist);
          if(zdist > 0) 
          {
            VMCIons.R[ion_index[j]][2] = ionR[ionidx][2] - std::sqrt(h*h+zdist*zdist);
            VMCIons.update();
            updateCoeff();
            W.update();	  
          }
          if(zdist < 0)
          {
            VMCIons.R[ion_index[j]][2] = ionR[ionidx][2] + std::sqrt(h*h+zdist*zdist);
            VMCIons.update();
            updateCoeff();
            W.update();	  
          }

          W.makeRotate(thisWalker,deltaR,ionR[ionidx],rangle,i);   //rotate forward
          wfs[i][0]=std::exp(Psi.evaluateLog(W)); //eval wavefunction
          W.makeRotate(thisWalker,deltaR,ionR[ionidx],-rangle,i);   //rotate  backward

          W.makeRotate(thisWalker,deltaR,ionR[ionidx],-rangle,i);   //rotate backward
          wfs[i][1]=std::exp(Psi.evaluateLog(W)); //eval wavefunction
          W.makeRotate(thisWalker,deltaR,ionR[ionidx],rangle,i);   //rotate backward

          VMCIons.R[ion_index[j]][2] = ionR[j][2];
          VMCIons.update();
          updateCoeff();
          W.update();	  
        }

        for (int i=0;i<3;++i)
        {
          VMCIons.R[ion_index[j]][i] = ionR[j][i]-h;
          VMCIons.update();
          updateCoeff();
          W.update();	  
          np1 = VMCIons.R[ion_index[0]]; np2 = VMCIons.R[ion_index[1]];
          dorotateshift( np1,np2,wfs[i][0]);
          wfs[i][0] = std ::exp(wfs[i][0]);

          VMCIons.R[ion_index[j]][i] = ionR[j][i]+h;
          VMCIons.update();
          updateCoeff();
          W.update();	  
          np1 = VMCIons.R[ion_index[0]]; np2 = VMCIons.R[ion_index[1]];
          dorotateshift( np1,np2,wfs[i][1]);
          wfs[i][1] = std ::exp(wfs[i][1]);

          VMCIons.R[ion_index[j]][i] = ionR[j][i];
          VMCIons.update();
          updateCoeff();
          W.update();	  	  
        }
        //ionsKineticE += IonKineticEnergy(wfs,f,h,VMCIons.Mass[j]);
        //ionsKineticE += IonKineticEnergy3(wfs,wfs2,f,f2,h,VMCIons.Mass[j]);
        
        ionsKineticE += IonKineticEnergy3_2(wfs,f,h,VMCIons.Mass[ion_index[j]],ion_index[j],VMCIons.R);//should also give the gradient as output (however, not at the moment)

	    } // end ionsKinetic 

      np1 = VMCIons.R[ion_index[0]]; np2 = VMCIons.R[ion_index[1]];
      dorotateshift( np1,np2,logpsi);

      RealType eloc=H.evaluate(W)+ionsKineticE;
      thisWalker.R = W.R;

      for (int i=0;i<thisWalker.ionPos.size();++i)
      {
        W.ionPos[i]=ionR[i];
        thisWalker.ionPos[i]=ionR[i];
      }	
      //testing   testing testing
      PosType ntneworigin, ntnewsec,ntdispl,nttestsec;
      //cout << "new position 1 " <<  ionR[0] << endl;
      //cout << "new position 2 " <<  ionR[1] << endl;
      ntneworigin = ionR[0] - ionR[0];
      ntnewsec = ionR[0] - ionR[0];
      ntdispl = ionR[1] - ionR[0];
      RealType ntodist= std::sqrt(ntdispl[0]*ntdispl[0]+ntdispl[1]*ntdispl[1]+ntdispl[2]*ntdispl[2]);
      ntnewsec[2] = ntodist;  //align sec along z direction   
      W.ionPos[0]=ntneworigin;
      W.ionPos[1]=ntnewsec;
      thisWalker.ionPos[0]= ntneworigin;
      thisWalker.ionPos[1]= ntnewsec;
      // cout << "rotate position 1 " <<  ntneworigin << endl; 
      // cout << "rotate position 2 " <<  ntnewsec << endl;
      // cout << "eloc " <<  eloc <<" eval  " << H.evaluate(W) << " ik " << ionsKineticE <<endl;
      VMCIons.R[ion_index[0]] = ntneworigin;
      VMCIons.R[ion_index[1]] = ntnewsec;
      VMCIons.update();
      updateCoeff();
      H.update_source(VMCIons);
      W.update();
      W.makeShiftRotate(deltaR,ionR[0],ionR[1]);
      dorotateshift(ntneworigin,ntnewsec,logpsi);
      //end testing   testing testing  
      W.saveWalker(thisWalker);

      thisWalker.resetProperty(logpsi,Psi.getPhase(),eloc);
      H.auxHevaluate(W,thisWalker);
      H.saveProperty(thisWalker.getPropertyBase());
      ++nAccept;
    } // end if move accepted

  } // end for each walker

  for (int i=0;i<Ndist;++i)
    nuclei_dist_step[i] += etaisyys[i]/et_count;
  

}

}

/***************************************************************************
 * $RCSfile: VMCUpdateAll.cpp,v $   $Author: jnkim $
 * $Revision: 1.25 $   $Date: 2006/10/18 17:03:05 $
 * $Id: VMCUpdateAll.cpp,v 1.25 2006/10/18 17:03:05 jnkim Exp $
 ***************************************************************************/
