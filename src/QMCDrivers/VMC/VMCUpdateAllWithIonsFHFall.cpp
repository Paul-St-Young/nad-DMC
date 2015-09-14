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
      for (int i=0;i<VMCIons.R.size();++i)
	{
	  if (VMCIons.Mass[i]>0.0) 
	    {
	      W.ionPos[i]=VMCIons.R[i];
              thisWalker.ionPos[i]=VMCIons.R[i];
	    }
	}
      W.saveWalker(thisWalker);
    }
  


}

VMCUpdateAllWithIons::~VMCUpdateAllWithIons()
{
}

void VMCUpdateAllWithIons::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
{
  for (; it!= it_end; ++it)
  {
    MCWalkerConfiguration::Walker_t& thisWalker(**it);

    Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
    W.loadWalker(thisWalker,true);

    makeGaussRandomWithEngine(deltaR,RandomGen);

    //FHF all

    PosType ionR[thisWalker.ionPos.size()];
    for (int i=0;i<thisWalker.ionPos.size();++i)
      {
	ionR[i]=thisWalker.ionPos[i];
      }    
    
    PosType ionCenter(0.0);

    RealType coeff_h = 1.0;
    RealType coeff_ff = 20.0;
    RealType rc_h = 0.0;
    RealType rc_ff = 4.3;

    
    for (int i=0;i<3;++i)
      ionCenter[i]=(ionR[1][i]+ionR[2][i])/2;
    
    RealType f2old=testGaussian(ionR[0],ionCenter,coeff_h,rc_h)
      +testGaussian(ionR[1],ionR[2],coeff_ff,rc_ff);

    for (int i=0;i<thisWalker.ionPos.size();++i)
      {
	ionR[i]=ionR[i]+std::sqrt(Tau/VMCIons.Mass[i])*deltaR[i];
	VMCIons.R[i] = ionR[i];
      }    

    makeGaussRandomWithEngine(deltaR,RandomGen);

    //if (!W.makeMove(thisWalker,deltaR, m_sqrttau))
    if (!W.makeMove(thisWalker,deltaR,SqrtTauOverMass))
    {

      for (int i=0;i<thisWalker.ionPos.size();++i)
      {
	VMCIons.R[i] = thisWalker.ionPos[i];
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


    for (int i=0;i<3;++i)
      ionCenter[i]=(ionR[1][i]+ionR[2][i])/2;
    
    RealType f2=testGaussian(ionR[0],ionCenter,coeff_h,rc_h)
      +testGaussian(ionR[1],ionR[2],coeff_ff,rc_ff);

    RealType g2 = f2*f2/f2old/f2old;
    
    //W.R = m_sqrttau*deltaR + thisWalker.R;
    //W.update();
    RealType logpsi(Psi.evaluateLog(W));
    RealType g= std::exp(2.0*(logpsi-thisWalker.Properties(LOGPSI)))*g2;
    if (RandomGen() > g)
    {
      thisWalker.Age++;
      ++nReject;
      for (int i=0;i<thisWalker.ionPos.size();++i)
	{
	  VMCIons.R[i] = thisWalker.ionPos[i];
	}
      VMCIons.update();
      H.update_source(VMCIons);
      W.update();
      H.rejectedMove(W,thisWalker);
    }
    else
    {

      
      RealType wfs[3][2];
      RealType wfs2[3][2];
      RealType f=std::exp(logpsi);
      RealType h=0.000001;
      
      RealType ionsKineticE=0.0;
      
      for (int j=0;j<thisWalker.ionPos.size();++j)
	{
	  for (int i=0;i<3;++i)
	    {
	      VMCIons.R[j][i] = ionR[j][i]-h;
	      VMCIons.update();
	      W.update();	  
	      wfs[i][0]=std::exp(Psi.evaluateLog(W));
	      for (int k=0;k<3;++k)
		ionCenter[k]=(VMCIons.R[1][k]+VMCIons.R[2][k])/2;
	      wfs2[i][0]=testGaussian(VMCIons.R[0],ionCenter,coeff_h,rc_h)
		+testGaussian(VMCIons.R[1],VMCIons.R[2],coeff_ff,rc_ff);
	      
	      VMCIons.R[j][i] = ionR[j][i]+h;
	      VMCIons.update();
	      W.update();	  
	      wfs[i][1]=std::exp(Psi.evaluateLog(W));
	      for (int k=0;k<3;++k)
		ionCenter[k]=(VMCIons.R[1][k]+VMCIons.R[2][k])/2;
	      wfs2[i][1]=testGaussian(VMCIons.R[0],ionCenter,coeff_h,rc_h)
		+testGaussian(VMCIons.R[1],VMCIons.R[2],coeff_ff,rc_ff);
	      
	      VMCIons.R[j][i] = ionR[j][i];
	      VMCIons.update();
	      W.update();	  	  
	    }
	  //ionsKineticE += IonKineticEnergy(wfs,f,h,VMCIons.Mass[j]);
	  ionsKineticE += IonKineticEnergy3(wfs,wfs2,f,f2,h,VMCIons.Mass[j]);
	}           
      
      RealType eloc=H.evaluate(W)+ionsKineticE;
      thisWalker.R = W.R;
      for (int i=0;i<thisWalker.ionPos.size();++i)
	{	  
	  thisWalker.ionPos[i]=ionR[i];
	}
      
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
