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

#include <iostream>

namespace qmcplusplus
{

VMCUpdateAllWithIons::VMCUpdateAllWithIons(MCWalkerConfiguration& w, TrialWaveFunction& psi,
					   QMCHamiltonian& h, RandomGenerator_t& rg, ParticleSet& ions,
					   std::vector<int> ionsToMove)
  : QMCUpdateBase(w,psi,h,rg), VMCIons(ions)
{
  UpdatePbyP=false;

  MCWalkerConfiguration::iterator
    wit(W.begin()), wit_end(W.end());

  // the ions to be moved (i.e. mass>0) need to be introduced first the ions ParticleSet
  // at the moment
  for(; wit != wit_end; ++wit)
    {
      Walker_t& thisWalker(**wit);
      Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
      W.loadWalker(thisWalker,true);
      ionRo.resize(thisWalker.ionPos.size());
      for (int i=0;i<VMCIons.R.size();++i)
	{
	  if (VMCIons.Mass[i]>0.0) 
	    {
	      W.ionPos[i]=VMCIons.R[i];
	      ionRo[i] = VMCIons.R[i];
              thisWalker.ionPos[i]=VMCIons.R[i];
	    }
	}
      W.saveWalker(thisWalker);
    }

  for (int i=0;i<ionsToMove.size();++i)
    ion_index.push_back(ionsToMove[i]);
  
  vector<OrbitalBase*> Z=psi.getOrbitals();
  vector<OrbitalBase*>::iterator it(Z.begin());
  vector<OrbitalBase*>::iterator it_end(Z.end());
  app_log() << " !!!!!!!!!!! updating determinant set " << endl;
  app_log() << " number of orbitals = " << Z.size() << endl;
  app_log() << " Z[0]=SlaterDet " << endl;
  app_log() << " Z[1]=two-body jastrow, Z[2]=one-body jastrow " << endl;
  app_log() << " SlaterDet.Dets[i] are DiracDeterminantBase" << endl;
  app_log() << " LogValue are returned by OrbitalSetTraits::evaluateLogAndPhase" << endl;
  //std::ostream myos; 
  //Z[0]->reportStatus(myos);
  //app_log() << myos << endl;

}

VMCUpdateAllWithIons::~VMCUpdateAllWithIons()
{
}



void VMCUpdateAllWithIons::dorotateshift(PosType& origin, PosType &second, RealType &getpsi )
{
   VMCIons.update();
   H.update_source(VMCIons);
   W.update();
  
 ///Calculate wavefunction
   RealType logpsi(Psi.evaluateLog(W));
   getpsi = logpsi;  
 return ;
}

void VMCUpdateAllWithIons::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
{

  PosType  np1,np2;

  for (; it!= it_end; ++it)
  {
    MCWalkerConfiguration::Walker_t& thisWalker(**it);

    Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
    W.loadWalker(thisWalker,true);

    // evaluate old nuclei_wfs
    PosType ionR[thisWalker.ionPos.size()];
    for (int i=0;i<thisWalker.ionPos.size();++i)
      {
	ionR[i]=thisWalker.ionPos[i];
	VMCIons.R[ion_index[i]] = ionR[i];
      }    
    RealType f2old=nuclei_wfs(VMCIons.R);

    // move the ions
    makeGaussRandomWithEngine(deltaR,RandomGen);
    for (int i=0;i<thisWalker.ionPos.size();++i)
      {	
	if(tauCountFreq[ion_index[i]]>1 && tauCount%tauCountFreq[ion_index[i]]) continue;	  
	RealType tau_eff=(RealType)tauCountFreq[ion_index[i]]*Tau;	
   //ntdeubg delete me (the factor 2)
	ionR[i]=ionR[i]+std::sqrt(tau_eff/VMCIons.Mass[ion_index[i]])*deltaR[i];
        
	VMCIons.R[ion_index[i]] = ionR[i];
      }    

    // move the electrons
    makeGaussRandomWithEngine(deltaR,RandomGen);
    if (!W.makeMove(thisWalker,deltaR,SqrtTauOverMass))
    {

      for (int i=0;i<thisWalker.ionPos.size();++i)
      {
	VMCIons.R[ion_index[i]] = thisWalker.ionPos[i];
      }
      VMCIons.update();
      H.update_source(VMCIons);
      W.update();

      H.rejectedMove(W,thisWalker);

      continue;
    }
    VMCIons.update();
    H.update_source(VMCIons);
    W.update();    

    // find acceptance ratio
    RealType f2=nuclei_wfs(VMCIons.R);
    RealType g2 = f2*f2/f2old/f2old;
    RealType logpsi;

     np1 = VMCIons.R[ion_index[0]]; np2 = VMCIons.R[ion_index[1]];
     dorotateshift( np1,np2,logpsi);

    RealType g= std::exp(2.0*(logpsi-thisWalker.Properties(LOGPSI)))*g2;

    // accept/reject
    if (RandomGen() > g)
    {
      thisWalker.Age++;
      ++nReject;
      for (int i=0;i<thisWalker.ionPos.size();++i)
	{
	  VMCIons.R[ion_index[i]] = thisWalker.ionPos[i];
	}
      VMCIons.update();
      H.update_source(VMCIons);
      W.update();
      H.rejectedMove(W,thisWalker);
    }
    else
    { // if accepted, calculate observables
      
      RealType wfs[3][2];
      wfs[0][0] = 0;
      wfs[0][1] = 0;
      wfs[1][0] = 0; 
      wfs[1][1] = 0; 
      wfs[2][0] =  0;
      wfs[2][1] =  0;
      //RealType wfs2[3][2];
      RealType f=std::exp(logpsi);
      RealType h=UniformGrid_h*1000; // defined in QMCUpdateBase.h
      
      RealType ionsKineticE=0.0;
      
      for (int j=0;j<thisWalker.ionPos.size();++j)
	{	  

	  for (int i=0;i<3;++i)
	    {
	      VMCIons.R[ion_index[j]][i] = ionR[j][i]-h;
	      VMCIons.update();
	      W.update();	  
//	      wfs[i][0]=std::exp(Psi.evaluateLog(W));
              np1 = VMCIons.R[ion_index[0]]; np2 = VMCIons.R[ion_index[1]];
              dorotateshift( np1,np2,wfs[i][0]);
              wfs[i][0] = std ::exp(wfs[i][0]);
	      //wfs2[i][0]=nuclei_wfs(VMCIons.R);
	      
	      VMCIons.R[ion_index[j]][i] = ionR[j][i]+h;
	      VMCIons.update();
	      W.update();	  
//	      wfs[i][1]=std::exp(Psi.evaluateLog(W));
              np1 = VMCIons.R[ion_index[0]]; np2 = VMCIons.R[ion_index[1]];
              dorotateshift( np1,np2,wfs[i][1]);
              wfs[i][1] = std ::exp(wfs[i][1]);
	      //wfs2[i][1]=nuclei_wfs(VMCIons.R);
	      
	      VMCIons.R[ion_index[j]][i] = ionR[j][i];
	      VMCIons.update();
	      W.update();	  	  
	    }
	  ionsKineticE += IonKineticEnergy3_2(wfs,f,h,VMCIons.Mass[ion_index[j]],ion_index[j],VMCIons.R);//should also give the gradient as output (however, not at the moment)

	}           

      np1 = VMCIons.R[ion_index[0]]; np2 = VMCIons.R[ion_index[1]];
      dorotateshift( np1,np2,logpsi);
      RealType eloc=H.evaluate(W)+ionsKineticE;
      thisWalker.R = W.R;

      for (int i=0;i<thisWalker.ionPos.size();++i)
	{
	  W.ionPos[i]=ionR[i];
	  thisWalker.ionPos[i]=ionR[i];
	}	
      W.saveWalker(thisWalker);
      
      thisWalker.resetProperty(logpsi,Psi.getPhase(),eloc);
      H.auxHevaluate(W,thisWalker);
      H.saveProperty(thisWalker.getPropertyBase());
      ++nAccept;
    }

  }

}

}

/***************************************************************************
 * $RCSfile: VMCUpdateAll.cpp,v $   $Author: jnkim $
 * $Revision: 1.25 $   $Date: 2006/10/18 17:03:05 $
 * $Id: VMCUpdateAll.cpp,v 1.25 2006/10/18 17:03:05 jnkim Exp $
 ***************************************************************************/
