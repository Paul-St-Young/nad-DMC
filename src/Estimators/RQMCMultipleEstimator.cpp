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
#include "Estimators/RQMCMultipleEstimator.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "Message/CommOperators.h"
#include "QMCDrivers/DriftOperators.h"

namespace qmcplusplus
{

/** constructor
 * @param h QMCHamiltonian to define the components
 * @param hcopy number of copies of QMCHamiltonians
 */
RQMCEstimator::RQMCEstimator(QMCHamiltonian& h, int hcopy) : CurrentWalker(0)
{
  NumCopies=hcopy;
  NumOperators = h.sizeOfObservables();
  FirstHamiltonian=h.startIndex();
  esum.resize(NumCopies,LE_INDEX);
  esum_name.resize(LE_INDEX,NumCopies);
  elocal.resize(NumCopies,NumOperators);
  elocal_name.resize(NumCopies,NumOperators);
  char aname[32];
  //make the name tables
  //(localenergy, variance, weight)* + (each hamiltonian term)*
  for(int i=0; i<NumCopies; i++)
  {
    sprintf(aname,"LE%i",i);
    esum_name(ENERGY_INDEX,i)=aname;
    sprintf(aname,"Var%i",i);
    esum_name(ENERGY_SQ_INDEX,i)=aname;
    sprintf(aname,"WPsi%i",i);
    esum_name(WEIGHT_INDEX,i)=aname;
    for(int j=0; j<NumOperators; j++)
    {
      //sprintf(aname,"%s%i",h.getName(j).c_str(),i);
      sprintf(aname,"%s%i",h.getObservableName(j).c_str(),i);
      elocal_name(i,j)=aname;
    }
  }
  int ipair=0;
  for(int i=0; i<NumCopies-1; i++)
  {
    for(int j=i+1; j<NumCopies; j++)
    {
      sprintf(aname,"DiffS%iS%i",i,j);
      ediff_name.push_back(aname);
    }
  }
}

RQMCEstimator::RQMCEstimator(const RQMCEstimator& rest):
  ScalarEstimatorBase(rest),
  NumCopies(rest.NumCopies), NumOperators(rest.NumOperators),
  FirstHamiltonian(rest.FirstHamiltonian), esum(rest.esum),
  elocal(rest.elocal),elocal_name(rest.elocal_name)
{
}

ScalarEstimatorBase* RQMCEstimator::clone()
{
  return new RQMCEstimator(*this);
}
// new initialize
void RQMCEstimator::initialize(MultiChain* reptileRef, int setDirect, double setTau, int setSteps)
{
  Reptile = reptileRef;
  AveEloc.resize(NumCopies);
  AveWeight.resize(NumCopies);
  for(int ipsi=0; ipsi<NumCopies; ipsi++)
  {
    AveEloc[ipsi]=0.0;
    AveWeight[ipsi]=0.0;
  }
  Tau = setTau;
  Directionless = setDirect;
  nSteps = setSteps;
}

//OLD and BROKEN
//  void
//    RQMCEstimator
//    ::initialize(MCWalkerConfiguration& W,
//        vector<QMCHamiltonian*>& h,
//        vector<TrialWaveFunction*>& psi,
//        RealType tau,vector<RealType>& Norm,
//        bool require_register) {
//
//      NumWalkers = W.getActiveWalkers();
//      //allocate UmbrellaEnergy
//      int numPtcls(W.getTotalNum());
//      RatioIJ.resize(NumWalkers,NumCopies*(NumCopies-1)/2);
//
//      MCWalkerConfiguration::iterator it(W.begin());
//      MCWalkerConfiguration::iterator it_end(W.end());
//
//      vector<RealType> sumratio(NumCopies), logpsi(NumCopies);
//      int iw(0);
//      int DataSetSize((*it)->DataSet.size());
//      while(it != it_end) {
//
//        Walker_t& thisWalker(**it);
//        (*it)->DataSet.rewind();
//
//        if(require_register) {
//          W.registerData(thisWalker,(*it)->DataSet);
//        } else {
//          W.R = thisWalker.R;
//          W.update();
//          if(DataSetSize) W.updateBuffer((*it)->DataSet);
//        }
//
//        //evalaute the wavefunction and hamiltonian
//        for(int ipsi=0; ipsi< NumCopies;ipsi++) {
//          psi[ipsi]->G.resize(numPtcls);
//          psi[ipsi]->L.resize(numPtcls);
//          //Need to modify the return value of OrbitalBase::registerData
//          if(require_register) {
//            logpsi[ipsi]=psi[ipsi]->registerData(W,(*it)->DataSet);
//          } else {
//            if(DataSetSize)logpsi[ipsi]=psi[ipsi]->updateBuffer(W,(*it)->DataSet);
//            else logpsi[ipsi]=psi[ipsi]->evaluateLog(W);
//          }
//          psi[ipsi]->G=W.G;
//          thisWalker.Properties(ipsi,LOGPSI)=logpsi[ipsi];
//          thisWalker.Properties(ipsi,LOCALENERGY)=h[ipsi]->evaluate(W);
//          h[ipsi]->saveProperty(thisWalker.getPropertyBase(ipsi));
//          sumratio[ipsi]=1.0;
//        }
//
//        //Check SIMONE's note
//        //Compute the sum over j of Psi^2[j]/Psi^2[i] for each i
//        int indexij(0);
//        RealType *rPtr=RatioIJ[iw];
//        for(int ipsi=0; ipsi< NumCopies-1; ipsi++) {
//          for(int jpsi=ipsi+1; jpsi< NumCopies; jpsi++){
//            RealType r= std::exp(2.0*(logpsi[jpsi]-logpsi[ipsi]));
//            rPtr[indexij++]=r*Norm[ipsi]/Norm[jpsi];
//            sumratio[ipsi] += r;
//            sumratio[jpsi] += 1.0/r;
//          }
//        }
//
//        //Re-use Multiplicity as the sumratio
//        thisWalker.Multiplicity=sumratio[0];
//        //DON't forget DRIFT!!!
//        thisWalker.Drift=0.0;
//
//        for(int ipsi=0; ipsi< NumCopies; ipsi++) {
//          RealType wgt=1.0/sumratio[ipsi];
//          thisWalker.Properties(ipsi,UMBRELLAWEIGHT)=wgt;
//
//          //thisWalker.Drift += wgt*psi[ipsi]->G;
//          PAOps<RealType,DIM>::axpy(wgt,psi[ipsi]->G,thisWalker.Drift);
//        }
//        thisWalker.Drift *= tau;
//        ++it;++iw;
//      }
//    }
//
//  void
//    RQMCEstimator
//    ::initialize(MCWalkerConfiguration& W, vector<ParticleSet*>& WW,
//        SpaceWarp& Warp,
//        vector<QMCHamiltonian*>& h,
//        vector<TrialWaveFunction*>& psi,
//        RealType tau,vector<RealType>& Norm,
//        bool require_register) {
//
//      NumWalkers = W.getActiveWalkers();
//
//      int numPtcls(W.getTotalNum());
//
//      RatioIJ.resize(NumWalkers,NumCopies*(NumCopies-1)/2);
//
//      vector<RealType> invsumratio(NumCopies);
//      MCWalkerConfiguration::ParticlePos_t drift(numPtcls);
//
//      MCWalkerConfiguration::iterator it(W.begin());
//      MCWalkerConfiguration::iterator it_end(W.end());
//
//      vector<RealType> sumratio(NumCopies), logpsi(NumCopies);
//      vector<RealType> Jacobian(NumCopies);
//
//      int jindex=W.addProperty("Jacobian");
//      int iw(0);
//      int DataSetSize((*it)->DataSet.size());
//      while(it != it_end) {
//        Walker_t& thisWalker(**it);
//        (*it)->DataSet.rewind();
//
//        //NECESSARY SINCE THE DISTANCE TABLE OF W ARE USED TO WARP
//        if(require_register) {
//          W.registerData(thisWalker,(*it)->DataSet);
//        } else {
//          W.R = thisWalker.R;
//          W.update();
//          if(DataSetSize) W.updateBuffer((*it)->DataSet);
//        }
//
//        for(int ipsi=0; ipsi<NumCopies; ipsi++) Jacobian[ipsi]=1.e0;
//        for(int iptcl=0; iptcl< numPtcls; iptcl++){
//          Warp.warp_one(iptcl,0);
//          for(int ipsi=0; ipsi<NumCopies; ipsi++){
//            WW[ipsi]->R[iptcl]=W.R[iptcl]+Warp.get_displacement(iptcl,ipsi);
//            Jacobian[ipsi]*=Warp.get_Jacobian(iptcl,ipsi);
//          }
//          if(require_register || DataSetSize) Warp.update_one_ptcl_Jacob(iptcl);
//        }
//
//        for(int ipsi=0; ipsi<NumCopies; ipsi++){
//          thisWalker.Properties(ipsi,jindex)=Jacobian[ipsi];
//        }
//
//        //update distance table and bufferize it if necessary
//        if(require_register) {
//          for(int ipsi=0; ipsi<NumCopies; ipsi++){
//            WW[ipsi]->registerData((*it)->DataSet);
//          }
//          Warp.registerData(WW,(*it)->DataSet);
//        } else {
//          for(int ipsi=0; ipsi<NumCopies; ipsi++){
//            WW[ipsi]->update();
//            if(DataSetSize) WW[ipsi]->updateBuffer((*it)->DataSet);
//          }
//          if(DataSetSize) Warp.updateBuffer((*it)->DataSet);
//        }
//
//
//
//        //evalaute the wavefunction and hamiltonian
//        for(int ipsi=0; ipsi< NumCopies;ipsi++) {
//          psi[ipsi]->G.resize(numPtcls);
//          psi[ipsi]->L.resize(numPtcls);
//          //Need to modify the return value of OrbitalBase::registerData
//          if(require_register) {
//            logpsi[ipsi]=psi[ipsi]->registerData(*WW[ipsi],(*it)->DataSet);
//          }else{
//            if(DataSetSize)logpsi[ipsi]=psi[ipsi]->updateBuffer(*WW[ipsi],(*it)->DataSet);
//            else logpsi[ipsi]=psi[ipsi]->evaluateLog(*WW[ipsi]);
//          }
//          psi[ipsi]->G=WW[ipsi]->G;
//          thisWalker.Properties(ipsi,LOGPSI)=logpsi[ipsi];
//          thisWalker.Properties(ipsi,LOCALENERGY)=h[ipsi]->evaluate(*WW[ipsi]);
//          h[ipsi]->saveProperty(thisWalker.getPropertyBase(ipsi));
//          sumratio[ipsi]=1.0;
//        }
//
//        //Check SIMONE's note
//        //Compute the sum over j of Psi^2[j]/Psi^2[i] for each i
//        int indexij(0);
//        RealType *rPtr=RatioIJ[iw];
//        for(int ipsi=0; ipsi< NumCopies-1; ipsi++) {
//          for(int jpsi=ipsi+1; jpsi< NumCopies; jpsi++){
//            RealType r= std::exp(2.0*(logpsi[jpsi]-logpsi[ipsi]))*Norm[ipsi]/Norm[jpsi];
//            //BEWARE: RatioIJ DOES NOT INCLUDE THE JACOBIANS!
//            rPtr[indexij++]=r;
//            r*=(Jacobian[jpsi]/Jacobian[ipsi]);
//            sumratio[ipsi] += r;
//            sumratio[jpsi] += 1.0/r;
//          }
//        }
//
//        //Re-use Multiplicity as the sumratio
//        thisWalker.Multiplicity=sumratio[0];
//
//        /*START COMMENT
//          QMCTraits::PosType WarpDrift;
//          RealType denom(0.e0),wgtpsi;
//          thisWalker.Drift=0.e0;
//          for(int ipsi=0; ipsi< NumCopies; ipsi++) {
//          wgtpsi=1.e0/sumratio[ipsi];
//          thisWalker.Properties(ipsi,UMBRELLAWEIGHT)=wgtpsi;
//          denom += wgtpsi;
//          for(int iptcl=0; iptcl< numPtcls; iptcl++){
//          WarpDrift=dot( psi[ipsi]->G[iptcl], Warp.get_Jacob_matrix(iptcl,ipsi)  )
//          +5.0e-1*Warp.get_grad_ln_Jacob(iptcl,ipsi) ;
//          thisWalker.Drift[iptcl] += (wgtpsi*WarpDrift);
//          }
//          }
//        //Drift = denom*Drift;
//        thisWalker.Drift *= (tau/denom);
//        END COMMENT*/
//        for(int ipsi=0; ipsi< NumCopies ;ipsi++){
//          invsumratio[ipsi]=1.0/sumratio[ipsi];
//          thisWalker.Properties(ipsi,UMBRELLAWEIGHT)=invsumratio[ipsi];
//        }
//        setScaledDrift(tau,psi[0]->G,drift);
//        thisWalker.Drift=invsumratio[0]*drift;
//        for(int ipsi=1; ipsi< NumCopies ;ipsi++) {
//          setScaledDrift(tau,psi[ipsi]->G,drift);
//          thisWalker.Drift += (invsumratio[ipsi]*drift);
//        }
//        ++it;++iw;
//      }
//    }
//
/**  add the local energy, variance and all the Hamiltonian components to the scalar record container
 *@param record storage of scalar records (name,value)
 */
void
RQMCEstimator::add2Record(RecordNamedProperty<RealType>& record)
{
  if(ediff_name.size())
  {
    FirstColumnIndex = record.add(ediff_name[0].c_str());
    for(int i=1; i<ediff_name.size(); i++)
      record.add(ediff_name[i].c_str());
    for(int i=0; i<esum_name.size(); i++)
      record.add(esum_name(i).c_str());
  }
  else
  {
    FirstColumnIndex = record.add(esum_name(0).c_str());
    for(int i=1; i<esum_name.size(); i++)
      record.add(esum_name(i).c_str());
  }
  for(int i=0; i<elocal_name.size(); i++)
    record.add(elocal_name(i).c_str());
  //FirstColumnIndex = record.add(esum_name(0).c_str());
  //for(int i=1; i<esum_name.size(); i++) record.add(esum_name(i).c_str());
  //for(int i=0; i<ediff_name.size(); i++) record.add(ediff_name[i].c_str());
  //for(int i=0; i<elocal_name.size(); i++) record.add(elocal_name(i).c_str());
}

void
RQMCEstimator::accumulate(const Walker_t& awalker, RealType wgt)
{
  // taken from RQMCMultiple.cpp
  for(int ipsi=0; ipsi<NumCopies; ipsi++)
  {
    //cerr << "  " << ipsi+1 << " of " << NumCopies;
    double WeightedEloc=Reptile->UmbrellaWeight[ipsi]*
                        ( Reptile->front()->Action(ipsi,Directionless)
                          +Reptile->back()->Action(ipsi,Directionless) );
    AveEloc[ipsi]+=WeightedEloc;
    AveWeight[ipsi]+=Reptile->UmbrellaWeight[ipsi];
    //cerr << "  Estim8 " << ipsi << " " << AveEloc[ipsi] << " " << AveWeight[ipsi] << endl;
  }
  ////const RealType* restrict etot=UmbrellaEnergy[CurrentWalker];
  ////const RealType* restrict wsum=UmbrellaWeight[CurrentWalker];
  //int iloc=0;
  //RealType *restrict esumPtr = esum.data();
  //for(int i=0; i<NumCopies; i++) {
  //  //get the pointer to the i-th row
  //  const RealType* restrict prop=awalker.getPropertyBase(i);
  //  RealType invr = prop[UMBRELLAWEIGHT];
  //  RealType e = prop[LOCALENERGY];
  //  *esumPtr++ += invr*e;   //esum(i,ENERGY_INDEX)    += invr*e;
  //  *esumPtr++ += invr*e*e; //esum(i,ENERGY_SQ_INDEX) += invr*e*e; //how to variance
  //  *esumPtr++ += invr;     //esum(i,WEIGHT_INDEX)    += invr;
  //  //accumulate elocal(i,j) with the weight,
  //  //The content the Walker::DynProperties other than LOCALENERGY are NOT weighted.
  //  //const RealType *t= awalker.getEnergyBase(i);
  //  //for(int j=0; j<NumOperators; j++) {
  //  //  *elocPtr++ += invr*(*t++);
  //  //}
  //
  //  //const RealType *t= awalker.getPropertyBase(i)+FirstHamiltonian;
  //  for(int j=0,h=FirstHamiltonian; j<NumOperators; j++,h++) {
  //    elocal(iloc++) += invr*prop[h];
  //  }
  //}
  ++CurrentWalker;
  //reset to zero
  if(CurrentWalker == NumWalkers)
    CurrentWalker=0;
  //cerr << "Leaving accumulate" << endl;
}

void RQMCEstimator::registerObservables(vector<observable_helper*>& h5dec, hid_t gid)
{
  //IMPLEMENT for hdf5
}

//  ///Set CurrentWalker to zero so that accumulation is done in a vectorized way
//  void RQMCEstimator::reset() {
//    CurrentWalker=0;
//    elocal=0.0; esum=0.0;
//  }

//  void RQMCEstimator::copy2Buffer(BufferType& msg)
//  {
//    msg.put(esum.begin(),esum.end());
//  }
//
//  /** calculate the averages and reset to zero
//   * @param record a container class for storing scalar records (name,value)
//   * @param wgtinv the inverse weight
//   *
//   * Disable collection. RQMCEstimator does not need to communiate at all.
//   */
//  void RQMCEstimator::report(RecordNamedProperty<RealType>& record, RealType wgtinv)
//  {
//    // taken from RQMCMultiple.cpp
//    int ir(FirstColumnIndex);
//    for(int ipsi=0; ipsi<NumCopies; ipsi++)
//      AveEloc[ipsi]/=(AveWeight[ipsi]*Tau+numeric_limits<RealType>::epsilon());
//
//    // not sure if this is the correct format for more than 2 correlated samples
//    for(int ipsi=0; ipsi<NumCopies-1; ipsi++)
//      for(int jpsi=ipsi+1; jpsi<NumCopies; jpsi++)
//        record[ir++] = AveEloc[jpsi] - AveEloc[ipsi];
//
//    for(int ipsi=0; ipsi<NumCopies; ipsi++)
//      record[ir++] = AveEloc[ipsi];
//
//    for(int ipsi=0; ipsi<NumCopies; ipsi++)
//      record[ir++] = 0.0;
//
//    for(int ipsi=0; ipsi<NumCopies; ipsi++)
//      record[ir++] = AveWeight[ipsi]/nSteps;
//
//    // reset containers
//    for(int ipsi=0; ipsi<NumCopies; ipsi++){
//      AveEloc[ipsi]=0.0;
//      AveWeight[ipsi]=0.0;
//    }
//
//    //for(int i=0; i<NumCopies; i++) {
//    //  RealType r = 1.0/esum(i,WEIGHT_INDEX);
//    //  RealType e = esum(i,ENERGY_INDEX)*r;
//    //  esum(i,ENERGY_INDEX)=e;
//    //  esum(i,ENERGY_SQ_INDEX)=esum(i,ENERGY_SQ_INDEX)*r-e*e;
//    //  esum(i,WEIGHT_INDEX)*=wgtinv;
//    //}
//
//    //ediff
//    //int ir(FirstColumnIndex);
//    //for(int i=0; i<NumCopies-1; i++) {
//    //  for(int j=i+1; j<NumCopies; j++) {
//    /// this was changed...
//    //    record[ir++]=AveEloc[i];
//    //  }
//    //}
//    //+(localenergy, variance, weight)*
//    //swap the row/column indices for (localenergy*, variance*, weight*)
//    //for(int l=0; l<LE_INDEX;l++) {
//    //  for(int i=0; i<NumCopies; i++)
//    //    record[ir++]=esum(i,l);
//    //}
//
//    //(each hamiltonian term)*
//    //int iloc=0;
//    //for(int i=0; i<NumCopies; i++) {
//    //  for(int j=0; j<NumOperators; j++, ir++) {
//    //    record[ir]=elocal(iloc++)*wgtinv/esum(i,WEIGHT_INDEX);
//    //  }
//    //}
//    /*for(int i=0; i<elocal.size(); i++, ir++) {
//      record[ir]=wgtinv*elocal(i);
//      }*/
//    //int n(elocal.size());
//    //const RealType* restrict eit(elocal.data());
//    //while(n) {record[ir++]=wgtinv*(*eit++);--n;}
//
//    //set the ScalarEstimator<T>::b_average and b_variace
//    //b_average = esum(0,ENERGY_INDEX);
//    //b_variance = esum(0,ENERGY_SQ_INDEX);
//
//    reset();
//  }
//  void RQMCEstimator::report(RecordNamedProperty<RealType>& record, RealType wgtinv,
//      BufferType& msg) {
//    msg.get(esum.begin(),esum.end());
//    report(record,wgtinv);
//  }
}
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 19:14:53 -0500 (Thu, 25 Apr 2013) $
 * $Id: RQMCMultipleEstimator.cpp 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/
