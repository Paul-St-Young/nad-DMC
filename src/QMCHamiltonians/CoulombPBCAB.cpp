//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim and Kris Delaney
//////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCHamiltonians/CoulombPBCAB.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "Message/Communicate.h"
#include "Utilities/ProgressReportEngine.h"

namespace qmcplusplus
{

CoulombPBCAB::CoulombPBCAB(ParticleSet& ions, ParticleSet& elns,
                                   bool computeForces):
  PtclA(ions), myConst(0.0), myGrid(0),V0(0),ComputeForces(computeForces),
  ForceBase (ions, elns), MaxGridPoints(10000)
{
  // if (ComputeForces)
  // 	InitVarReduction (0.5, 0, 3);
  ReportEngine PRE("CoulombPBCAB","CoulombPBCAB");
  //Use singleton pattern
  //AB = new LRHandlerType(ions);
  myTableIndex=elns.addTable(ions);
  initBreakup(elns);
  prefix="Flocal";
  app_log() << "  Maximum K shell " << AB->MaxKshell << endl;
  app_log() << "  Number of k vectors " << AB->Fk.size() << endl;
}

QMCHamiltonianBase* CoulombPBCAB::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  CoulombPBCAB* myclone=new CoulombPBCAB(PtclA,qp,ComputeForces);
  myclone->FirstForceIndex = FirstForceIndex;
  if(myGrid)
    myclone->myGrid=new GridType(*myGrid);
  for(int ig=0; ig<Vspec.size(); ++ig)
  {
    if(Vspec[ig])
    {
      RadFunctorType* apot=Vspec[ig]->makeClone();
      myclone->Vspec[ig]=apot;
      for(int iat=0; iat<PtclA.getTotalNum(); ++iat)
      {
        if(PtclA.GroupID[iat]==ig)
          myclone->Vat[iat]=apot;
      }
    }
  }
  return myclone;
}

CoulombPBCAB:: ~CoulombPBCAB()
{
  //probably need to clean up
}

void CoulombPBCAB::resetTargetParticleSet(ParticleSet& P)
{
  int tid=P.addTable(PtclA);
  if(tid != myTableIndex)
  {
    APP_ABORT("CoulombPBCAB::resetTargetParticleSet found inconsistent table index");
  }
  AB->resetTargetParticleSet(P);
}

void CoulombPBCAB::addObservables(PropertySetType& plist, BufferType& collectables)
{
  myIndex=plist.add(myName.c_str());
  if (ComputeForces)
    addObservablesF(plist);
}

CoulombPBCAB::Return_t
CoulombPBCAB::evaluate(ParticleSet& P)
{
  if (ComputeForces)
  {
    forces = 0.0;
    Value = evalLRwithForces(P) + evalSRwithForces(P) +myConst;
  }
  else
    Value = evalLR(P) + evalSR(P) +myConst;
  return Value;
}

CoulombPBCAB::Return_t
CoulombPBCAB::evalSR(ParticleSet& P)
{
  const DistanceTableData &d_ab(*P.DistTables[myTableIndex]);
  RealType res=0.0;
  //Loop over distinct eln-ion pairs
  for(int iat=0; iat<NptclA; iat++)
  {
    RealType esum = 0.0;
    RadFunctorType* rVs=Vat[iat];
    for(int nn=d_ab.M[iat], jat=0; nn<d_ab.M[iat+1]; ++nn,++jat)
    {
      // if(d_ab.r(nn)>=(myRcut-0.1)) continue;
      esum += Qat[jat]*d_ab.rinv(nn)*rVs->splint(d_ab.r(nn));;
    }
    //Accumulate pair sums...species charge for atom i.
    res += Zat[iat]*esum;
  }
  return res;
}


CoulombPBCAB::Return_t
CoulombPBCAB::evalLR(ParticleSet& P)
{
  const int slab_dir=OHMMS_DIM-1;
  RealType res=0.0;
  const StructFact& RhoKA(*(PtclA.SK));
  const StructFact& RhoKB(*(P.SK));
  if(RhoKA.SuperCellEnum==SUPERCELL_SLAB)
  {
    const DistanceTableData &d_ab(*P.DistTables[myTableIndex]);
    for(int iat=0; iat<NptclA; ++iat)
    {
      RealType u=0;
#if !defined(USE_REAL_STRUCT_FACTOR)
      for(int nn=d_ab.M[iat], jat=0; nn<d_ab.M[iat+1]; ++nn,++jat)
        u += Qat[jat]*AB->evaluate_slab(d_ab.dr(nn)[slab_dir], RhoKA.KLists.kshell, RhoKA.eikr[iat], RhoKB.eikr[jat]);
#endif
      res += Zat[iat]*u;
    }
  }
  else
  {
    for(int i=0; i<NumSpeciesA; i++)
    {
      RealType esum=0.0;
      for(int j=0; j<NumSpeciesB; j++)
      {
#if defined(USE_REAL_STRUCT_FACTOR)
        esum += Qspec[j]*AB->evaluate(RhoKA.KLists.kshell
                                      , RhoKA.rhok_r[i],RhoKA.rhok_i[i] , RhoKB.rhok_r[j],RhoKB.rhok_i[j]);
#else
        esum += Qspec[j]*AB->evaluate(RhoKA.KLists.kshell, RhoKA.rhok[i],RhoKB.rhok[j]);
#endif
      } //speceln
      res += Zspec[i]*esum;
    }
  }//specion
  return res;
}

/** Evaluate the background term. Other constants are handled by AA potentials.
 *
 * \f$V_{bg}^{AB}=-\sum_{\alpha}\sum_{\beta} N^{\alpha} N^{\beta} q^{\alpha} q^{\beta} v_s(k=0) \f$
 * @todo Here is where the charge system has to be handled.
 */
CoulombPBCAB::Return_t
CoulombPBCAB::evalConsts()
{
  RealType Consts=0.0;
  RealType vs_k0 = AB->evaluateSR_k0();
  for(int i=0; i<NumSpeciesA; i++)
  {
    RealType q=Zspec[i]*NofSpeciesA[i];
    for(int j=0; j<NumSpeciesB; j++)
    {
      Consts -= vs_k0*Qspec[j]*NofSpeciesB[j]*q;
    }
  }
  app_log() << "   Constant of PBCAB " << Consts << endl;
  return Consts;
}

void CoulombPBCAB::initBreakup(ParticleSet& P)
{
  SpeciesSet& tspeciesA(PtclA.getSpeciesSet());
  SpeciesSet& tspeciesB(P.getSpeciesSet());
  int ChargeAttribIndxA = tspeciesA.addAttribute("charge");
  int MemberAttribIndxA = tspeciesA.addAttribute("membersize");
  int ChargeAttribIndxB = tspeciesB.addAttribute("charge");
  int MemberAttribIndxB = tspeciesB.addAttribute("membersize");
  NptclA = PtclA.getTotalNum();
  NptclB = P.getTotalNum();
  NumSpeciesA = tspeciesA.TotalNum;
  NumSpeciesB = tspeciesB.TotalNum;
  //Store information about charges and number of each species
  Zat.resize(NptclA);
  Zspec.resize(NumSpeciesA);
  Qat.resize(NptclB);
  Qspec.resize(NumSpeciesB);
  NofSpeciesA.resize(NumSpeciesA);
  NofSpeciesB.resize(NumSpeciesB);
  for(int spec=0; spec<NumSpeciesA; spec++)
  {
    Zspec[spec] = tspeciesA(ChargeAttribIndxA,spec);
    NofSpeciesA[spec] = static_cast<int>(tspeciesA(MemberAttribIndxA,spec));
  }
  for(int spec=0; spec<NumSpeciesB; spec++)
  {
    Qspec[spec] = tspeciesB(ChargeAttribIndxB,spec);
    NofSpeciesB[spec] = static_cast<int>(tspeciesB(MemberAttribIndxB,spec));
  }
  RealType totQ=0.0;
  for(int iat=0; iat<NptclA; iat++)
    totQ+=Zat[iat] = Zspec[PtclA.GroupID[iat]];
  for(int iat=0; iat<NptclB; iat++)
    totQ+=Qat[iat] = Qspec[P.GroupID[iat]];
//    if(totQ>numeric_limits<RealType>::epsilon())
//    {
//      LOGMSG("PBCs not yet finished for non-neutral cells");
//      OHMMS::Controller->abort();
//    }
  ////Test if the box sizes are same (=> kcut same for fixed dimcut)
  kcdifferent = (std::abs(PtclA.Lattice.LR_kc - P.Lattice.LR_kc) > numeric_limits<RealType>::epsilon());
  minkc = std::min(PtclA.Lattice.LR_kc,P.Lattice.LR_kc);
  //AB->initBreakup(*PtclB);
  //initBreakup is called only once
  //AB = LRCoulombSingleton::getHandler(*PtclB);
  AB = LRCoulombSingleton::getHandler(P);
  myConst=evalConsts();
  myRcut=AB->get_rc();//Basis.get_rc();
  // create the spline function for the short-range part assuming pure potential
  if(V0==0)
  {
    V0 = LRCoulombSingleton::createSpline4RbyVs(AB,myRcut,myGrid);
    if(Vat.size())
    {
      app_log() << "  Vat is not empty. Something is wrong" << endl;
      OHMMS::Controller->abort();
    }
    Vat.resize(NptclA,V0);
    Vspec.resize(NumSpeciesA,0);//prepare for PP to overwrite it
  }
}

/** add a local pseudo potential
 * @param groupID species index
 * @param ppot radial functor for \f$rV_{loc}\f$ on a grid
 */
void CoulombPBCAB::add(int groupID, RadFunctorType* ppot)
{
  if(myGrid ==0)
  {
    myGrid = new LinearGrid<RealType>;
    int ng = min(MaxGridPoints, static_cast<int>(myRcut/1e-3)+1);
    app_log() << "    CoulombPBCAB::add \n Setting a linear grid=[0,"
              << myRcut << ") number of grid =" << ng << endl;
    myGrid->set(0,myRcut,ng);
  }
  if(Vspec[groupID]==0)
  {
    app_log() << "    Creating the short-range pseudopotential for species " << groupID << endl;
    int ng=myGrid->size();
    vector<RealType> v(ng);
    for(int ig=1; ig<ng-2; ig++)
    {
      RealType r=(*myGrid)[ig];
      //need to multiply r for the LR
      v[ig]=r*AB->evaluateLR(r)+ppot->splint(r);
    }
    v[0] = 2.0*v[1] - v[2];
    //by construction, v has to go to zero at the boundary
    v[ng-2]=0.0;
    v[ng-1]=0.0;
    RadFunctorType* rfunc=new RadFunctorType(myGrid,v);
    RealType deriv=(v[1]-v[0])/((*myGrid)[1]-(*myGrid)[0]);
    rfunc->spline(0,deriv,ng-1,0.0);
    Vspec[groupID]=rfunc;
    for(int iat=0; iat<NptclA; iat++)
    {
      if(PtclA.GroupID[iat]==groupID)
        Vat[iat]=rfunc;
    }
  }
  if (ComputeForces)
  {
    if(OHMMS::Controller->rank()==0)
    {
      FILE *fout = fopen ("Vlocal.dat", "w");
      for (double r=1.0e-8; r<myRcut; r+=1.0e-2)
      {
        double d_rV_dr, d2_rV_dr2;
        double Vr = Vat[0]->splint(r, d_rV_dr, d2_rV_dr2);
        Vr = Vat[0]->splint(r);
        fprintf (fout, "%1.8e %1.12e %1.12e %1.12e\n", r, Vr, d_rV_dr, d2_rV_dr2);
      }
      fclose(fout);
    }
  }
}

CoulombPBCAB::Return_t
CoulombPBCAB::registerData(ParticleSet& P, BufferType& buffer)
{
  P.SK->DoUpdate=true;
  SRpart.resize(NptclB);
  LRpart.resize(NptclB);
  Value=evaluateForPyP(P);
  buffer.add(SRpart.begin(),SRpart.end());
  buffer.add(LRpart.begin(),LRpart.end());
  buffer.add(Value);
  return Value;
}

/** The functions for PbyP move for reptation */
CoulombPBCAB::Return_t
CoulombPBCAB::updateBuffer(ParticleSet& P, BufferType& buffer)
{
  Value=evaluateForPyP(P);
  buffer.put(SRpart.begin(),SRpart.end());
  buffer.put(LRpart.begin(),LRpart.end());
  buffer.put(Value);
  return Value;
}

void CoulombPBCAB::copyFromBuffer(ParticleSet& P, BufferType& buffer)
{
  buffer.get(SRpart.begin(),SRpart.end());
  buffer.get(LRpart.begin(),LRpart.end());
  buffer.get(Value);
}

void CoulombPBCAB::copyToBuffer(ParticleSet& P, BufferType& buffer)
{
  buffer.put(SRpart.begin(),SRpart.end());
  buffer.put(LRpart.begin(),LRpart.end());
  buffer.put(Value);
}

CoulombPBCAB::Return_t
CoulombPBCAB::evaluateForPyP(ParticleSet& P)
{
  Return_t res=myConst;
#if defined(USE_REAL_STRUCT_FACTOR)
  APP_ABORT("CoulombPBCAB::evaluateForPyP(ParticleSet& P)");
#else
  SRpart=0.0;
  const DistanceTableData* d_ab=P.DistTables[myTableIndex];
  for(int iat=0; iat<NptclA; ++iat)
  {
    RealType z=Zat[iat];
    RadFunctorType* rVs=Vat[iat];
    for(int nn=d_ab->M[iat], jat=0; nn<d_ab->M[iat+1]; nn++,jat++)
    {
      RealType e=z*Qat[jat]*d_ab->rinv(nn)*rVs->splint(d_ab->r(nn));
      SRpart[jat]+=e;
      res+=e;
    }
  }
  LRpart=0.0;
  const StructFact& RhoKA(*(PtclA.SK));
  const StructFact& RhoKB(*(P.SK));
  // const StructFact& RhoKB(*(PtclB->SK));
  for(int i=0; i<NumSpeciesA; i++)
  {
    RealType z=Zspec[i];
    for(int jat=0; jat<P.getTotalNum(); ++jat)
    {
      RealType e=z*Qat[jat]*AB->evaluate(RhoKA.KLists.kshell, RhoKA.rhok[i],RhoKB.eikr[jat]);
      LRpart[jat]+=e;
      res+=e;
    }
  }
#endif
  return res;
}


CoulombPBCAB::Return_t
CoulombPBCAB::evaluatePbyP(ParticleSet& P, int active)
{
#if defined(USE_REAL_STRUCT_FACTOR)
  APP_ABORT("CoulombPBCAB::evaluatePbyP(ParticleSet& P, int active)");
#else
  const std::vector<DistanceTableData::TempDistType> &temp(P.DistTables[myTableIndex]->Temp);
  RealType q=Qat[active];
  SRtmp=0.0;
  for(int iat=0; iat<NptclA; ++iat)
  {
    SRtmp+=Zat[iat]*q*temp[iat].rinv1*Vat[iat]->splint(temp[iat].r1);
  }
  LRtmp=0.0;
  const StructFact& RhoKA(*(PtclA.SK));
  //const StructFact& RhoKB(*(PtclB->SK));
  const StructFact& RhoKB(*(P.SK));
  for(int i=0; i<NumSpeciesA; i++)
    LRtmp+=Zspec[i]*q*AB->evaluate(RhoKA.KLists.kshell, RhoKA.rhok[i],RhoKB.eikr_temp.data());
#endif
  return NewValue=Value+(SRtmp-SRpart[active])+(LRtmp-LRpart[active]);
  //return NewValue=Value+(SRtmp-SRpart[active]);
}

void CoulombPBCAB::acceptMove(int active)
{
  SRpart[active]=SRtmp;
  LRpart[active]=LRtmp;
  Value=NewValue;
}


CoulombPBCAB::Return_t
CoulombPBCAB::evalLRwithForces(ParticleSet& P)
{
  const StructFact& RhoKA(*(PtclA.SK));
  const StructFact& RhoKB(*(P.SK));
  vector<TinyVector<RealType,DIM> > grad(PtclA.getTotalNum());
  for(int j=0; j<NumSpeciesB; j++)
  {
    for (int iat=0; iat<grad.size(); iat++)
      grad[iat] = TinyVector<RealType,DIM>(0.0, 0.0, 0.0);
    AB->evaluateGrad(PtclA, P, j, Zat, grad);
    for (int iat=0; iat<grad.size(); iat++)
      forces[iat] += Qspec[j]*grad[iat];
  } // electron species
  return evalLR(P);
}

CoulombPBCAB::Return_t
CoulombPBCAB::evalSRwithForces(ParticleSet& P)
{
  const DistanceTableData &d_ab(*P.DistTables[myTableIndex]);
  RealType res=0.0;
  //Loop over distinct eln-ion pairs
  for(int iat=0; iat<NptclA; iat++)
  {
    RealType esum = 0.0;
    RadFunctorType* rVs=Vat[iat];
    for(int nn=d_ab.M[iat], jat=0; nn<d_ab.M[iat+1]; ++nn,++jat)
    {
      RealType rV, d_rV_dr, d2_rV_dr2, V;
      rV = rVs->splint(d_ab.r(nn), d_rV_dr, d2_rV_dr2);
      V = rV *d_ab.rinv(nn);
      PosType drhat = d_ab.rinv(nn) * d_ab.dr(nn);
      esum += Qat[jat]*d_ab.rinv(nn)*rV;
      forces[iat] += Zat[iat]*Qat[jat] * //g(d_ab.r(nn)) *
                     (d_rV_dr - V)*d_ab.rinv(nn) *drhat;
    }
    //Accumulate pair sums...species charge for atom i.
    res += Zat[iat]*esum;
  }
  return res;
}
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 5832 $   $Date: 2013-05-10 07:46:50 -0500 (Fri, 10 May 2013) $
 * $Id: CoulombPBCAB.cpp 5832 2013-05-10 12:46:50Z jnkim $
 ***************************************************************************/

