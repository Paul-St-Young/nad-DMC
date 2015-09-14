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
                           QMCHamiltonian& h, RandomGenerator_t& rg, ParticleSet& ions)
  : QMCUpdateBase(w,psi,h,rg), VMCIons(ions)
{
  UpdatePbyP=false;

  MCWalkerConfiguration::iterator
    wit(W.begin()), wit_end(W.end());

  for(; wit != wit_end; ++wit)
    {
      Walker_t& thisWalker(**wit);
      Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
      W.loadWalker(thisWalker,true);
      for (int i=0;i<VMCIons.R.size();++i)
	{
	  // Mass==0 is default for fixed nuclei (not used)
	  // i.e. if defined Mass>0 then those ions are moved
	  cout << "size " << VMCIons.R.size() << " Mass " <<  VMCIons.Mass[i] << endl;
	  if (VMCIons.Mass[i]>0.0) 
	    W.initializeIonPos(VMCIons.R[i]);
	  cout << "NumIons " << W.getNumIons() << endl;

	}

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

    cout << "W.getNumIons " << W.getNumIons() << endl;
    PosType ionR[W.getNumIons()];
    for (int i=0;i<W.getNumIons();++i)
      {
	ionR[i]=W.getIonPos(i);
	cout << "ion Pos at " << i << " is " << ionR[i] << endl;
	ionR[i]=ionR[i]+std::sqrt(Tau/VMCIons.Mass[i])*deltaR[i];
	VMCIons.R[i] = ionR[i];
      }    
    /*
    RealType etaisyys1=0.0;
    for (int i=0;i<3;++i)
      etaisyys1 += (ionR[0][i]-ionR[1][i])*(ionR[0][i]-ionR[1][i]);
    
    etaisyys1 = std::sqrt(etaisyys1);
    */

    //if (etaisyys1>1.6 || etaisyys1<1.2)
    //  continue;    

    makeGaussRandomWithEngine(deltaR,RandomGen);

    //if (!W.makeMove(thisWalker,deltaR, m_sqrttau))
    if (!W.makeMove(thisWalker,deltaR,SqrtTauOverMass))
    {
      H.rejectedMove(W,thisWalker);
      for (int i=0;i<W.getNumIons();++i)
      {
	ionR[i]=W.getIonPos(i);
	VMCIons.R[i] = ionR[i];
      }
      //VMCIons.update();
      //H.update_source(VMCIons);
      //W.update();
      continue;
    }
    VMCIons.update();
    H.update_source(VMCIons);
    W.update();    

    //W.R = m_sqrttau*deltaR + thisWalker.R;
    //W.update();
    RealType logpsi(Psi.evaluateLog(W));
    RealType g= std::exp(2.0*(logpsi-thisWalker.Properties(LOGPSI)));
    if (RandomGen() > g)
    {
      thisWalker.Age++;
      ++nReject;
      for (int i=0;i<W.getNumIons();++i)
	{
	  ionR[i]=W.getIonPos(i);
	  VMCIons.R[i] = ionR[i];
	}
      VMCIons.update();
      H.update_source(VMCIons);
      W.update();
      H.rejectedMove(W,thisWalker);
    }
    else
    {

      
      RealType wfs[3][2];
      RealType f=std::exp(logpsi);
      RealType h=0.001;
      
      RealType ionsKineticE=0.0;
      
      for (int j=0;j<W.getNumIons();++j)
	{
	  for (int i=0;i<3;++i)
	    {
	      VMCIons.R[j][i] = ionR[j][i]-h;
	      VMCIons.update();
	      W.update();	  
	      wfs[i][0]=std::exp(Psi.evaluateLog(W));
	      
	      VMCIons.R[j][i] = ionR[j][i]+h;
	      VMCIons.update();
	      W.update();	  
	      wfs[i][1]=std::exp(Psi.evaluateLog(W));
	      
	      VMCIons.R[j][i] = ionR[j][i];
	      VMCIons.update();
	      W.update();	  	  
	    }
	  ionsKineticE += IonKineticEnergy(wfs,f,h,VMCIons.Mass[j]);
	}           
      
      cout << "kinetic E " << ionsKineticE << endl;
      
      RealType eloc=H.evaluate(W)+ionsKineticE;
      thisWalker.R = W.R;
      W.setIonPosAll(ionR);

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
