//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim and Kris Delaney
//////////////////////////////////////////////////////////////////
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
#include "QMCHamiltonians/CoulombPBCAA.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "Utilities/ProgressReportEngine.h"
#include <numeric>

namespace qmcplusplus
{

CoulombPBCAA::CoulombPBCAA(ParticleSet& ref, bool active,
                                   bool computeForces) :
  AA(0), myGrid(0), rVs(0),
  is_active(active), FirstTime(true), myConst(0.0),
  ComputeForces(computeForces), ForceBase(ref,ref)
{
  ReportEngine PRE("CoulombPBCAA","CoulombPBCAA");

  //save source tag
  SourceID=ref.tag();

  //create a distance table: just to get the table name
  DistanceTableData *d_aa = DistanceTable::add(ref);
  PtclRefName=d_aa->Name;
  initBreakup(ref);
  prefix="F_AA";
  app_log() << "  Maximum K shell " << AA->MaxKshell << endl;
  app_log() << "  Number of k vectors " << AA->Fk.size() << endl;
  if(!is_active)
  {
    d_aa->evaluate(ref);
    update_source(ref);
    //RealType eL(0.0), eS(0.0);
    //if (computeForces)
    //{
    //  forces = 0.0;
    //  eS=evalSRwithForces(ref);
    //  // 1.3978248322
    //  eL=evalLRwithForces(ref);
    //  // 2.130267378
    //}
    //else
    //{
    //  eL=evalLR(ref);
    //  eS=evalSR(ref);
    //}
    //NewValue=Value = eL+eS+myConst;
    //app_log() << "  Fixed Coulomb potential for " << ref.getName();
    //app_log() << "\n    e-e Madelung Const. =" << MC0
    //          << "\n    Vtot     =" << Value << endl;
  }
  app_log() << "  Fixed Coulomb potential for " << ref.getName();
  app_log() << "\n    e-e Madelung Const. =" << MC0
            << "\n    Vtot     =" << Value << endl;
}

CoulombPBCAA:: ~CoulombPBCAA() { }

void CoulombPBCAA::addObservables(PropertySetType& plist, BufferType& collectables)
{
  addValue(plist);
  if (ComputeForces)
    addObservablesF(plist);
}

void CoulombPBCAA::update_source(ParticleSet& s)
{
  if(s.tag() == SourceID || s.parent() == SourceID)
  {
    RealType eL(0.0), eS(0.0);
    if (ComputeForces)
    {
      forces = 0.0;
      eS=evalSRwithForces(s);
      eL=evalLRwithForces(s);
    }
    else
    {
      eL=evalLR(s);
      eS=evalSR(s);
    }
    NewValue=Value = eL+eS+myConst;
  }
}

void CoulombPBCAA::resetTargetParticleSet(ParticleSet& P)
{
  if(is_active)
  {
    PtclRefName=P.DistTables[0]->Name;
    AA->resetTargetParticleSet(P);
  }
}

CoulombPBCAA::Return_t
CoulombPBCAA::evaluate(ParticleSet& P)
{
  if(is_active)
    Value = evalLR(P)+evalSR(P)+myConst;
  return Value;
}

CoulombPBCAA::Return_t
CoulombPBCAA::registerData(ParticleSet& P, BufferType& buffer)
{
  if(is_active)
  {
    P.SK->DoUpdate=true;
    SR2.resize(NumCenters,NumCenters);
    dSR.resize(NumCenters);
    del_eikr.resize(P.SK->KLists.numk);
    Value=evaluateForPbyP(P);
    buffer.add(SR2.begin(),SR2.end());
    buffer.add(Value);
  }
  return Value;
}

CoulombPBCAA::Return_t
CoulombPBCAA::updateBuffer(ParticleSet& P, BufferType& buffer)
{
  if(is_active)
  {
    Value=evaluateForPbyP(P);
    buffer.put(SR2.begin(),SR2.end());
    buffer.put(Value);
  }
  return Value;
}

void CoulombPBCAA::copyFromBuffer(ParticleSet& P, BufferType& buffer)
{
  if(is_active)
  {
    buffer.get(SR2.begin(),SR2.end());
    buffer.get(Value);
  }
}

void CoulombPBCAA::copyToBuffer(ParticleSet& P, BufferType& buffer)
{
  if(is_active)
  {
    buffer.put(SR2.begin(),SR2.end());
    buffer.put(Value);
  }
}

CoulombPBCAA::Return_t
CoulombPBCAA::evaluateForPbyP(ParticleSet& P)
{
  if(is_active)
  {
    SR2=0.0;
    Return_t res=myConst;
    const DistanceTableData *d_aa = P.DistTables[0];
    for(int iat=0; iat<NumCenters; iat++)
    {
      Return_t z=0.5*Zat[iat];
      for(int nn=d_aa->M[iat],jat=iat+1; nn<d_aa->M[iat+1]; ++nn,++jat)
      {
        Return_t e=z*Zat[jat]*d_aa->rinv(nn)*rVs->splint(d_aa->r(nn));
        SR2(iat,jat)=e;
        SR2(jat,iat)=e;
        res+=e+e;
      }
    }
    return res+evalLR(P);
  }
  else
    return Value;
}

CoulombPBCAA::Return_t
CoulombPBCAA::evaluatePbyP(ParticleSet& P, int active)
{
  if(is_active)
  {
    const std::vector<DistanceTableData::TempDistType> &temp(P.DistTables[0]->Temp);
    Return_t z=0.5*Zat[active];
    Return_t sr=0;
    const Return_t* restrict sr_ptr=SR2[active];
    for(int iat=0; iat<NumCenters; ++iat,++sr_ptr)
    {
      if(iat==active)
        dSR[active]=0.0;
      else
        sr+=dSR[iat]=(z*Zat[iat]*temp[iat].rinv1*rVs->splint(temp[iat].r1)- (*sr_ptr));
    }
#if defined(USE_REAL_STRUCT_FACTOR)
    APP_ABORT("CoulombPBCAA::evaluatePbyP");
#else
    const StructFact& PtclRhoK(*(P.SK));
    const ComplexType* restrict eikr_new=PtclRhoK.eikr_temp.data();
    const ComplexType* restrict eikr_old=PtclRhoK.eikr[active];
    ComplexType* restrict d_ptr=del_eikr.data();
    for(int k=0; k<del_eikr.size(); ++k)
      *d_ptr++ = (*eikr_new++ - *eikr_old++);
    int spec2=SpeciesID[active];
    for(int spec1=0; spec1<NumSpecies; ++spec1)
    {
      Return_t zz=z*Zspec[spec1];
      sr += zz*AA->evaluate(PtclRhoK.KLists.kshell,PtclRhoK.rhok[spec1],del_eikr.data());
    }
    sr+= z*z*AA->evaluate(PtclRhoK.KLists.kshell,del_eikr.data(),del_eikr.data());
    //// const StructFact& PtclRhoK(*(PtclRef->SK));
    //const StructFact& PtclRhoK(*(P.SK));
    //const ComplexType* restrict eikr_new=PtclRhoK.eikr_temp.data();
    //const ComplexType* restrict eikr_old=PtclRhoK.eikr[active];
    //ComplexType* restrict d_ptr=del_eikr.data();
    //for(int k=0; k<del_eikr.size(); ++k) *d_ptr++ = (*eikr_new++ - *eikr_old++);
    //for(int iat=0;iat<NumCenters; ++iat)
    //{
    //  if(iat!=active)
    //    sr += z*Zat[iat]*AA->evaluate(PtclRhoK.KLists.kshell, PtclRhoK.eikr[iat],del_eikr.data());
    //}
#endif
    return NewValue=Value+2.0*sr;
  }
  else
    return Value;
}

void CoulombPBCAA::acceptMove(int active)
{
  if(is_active)
  {
    Return_t* restrict sr_ptr=SR2[active];
    Return_t* restrict pr_ptr=SR2.data()+active;
    for(int iat=0; iat<NumCenters; ++iat, ++sr_ptr,pr_ptr+=NumCenters)
      *pr_ptr = *sr_ptr += dSR[iat];
    Value=NewValue;
  }
}

void CoulombPBCAA::initBreakup(ParticleSet& P)
{
  //SpeciesSet& tspecies(PtclRef->getSpeciesSet());
  SpeciesSet& tspecies(P.getSpeciesSet());
  //Things that don't change with lattice are done here instead of InitBreakup()
  ChargeAttribIndx = tspecies.addAttribute("charge");
  MemberAttribIndx = tspecies.addAttribute("membersize");
  NumCenters = P.getTotalNum();
  NumSpecies = tspecies.TotalNum;
  Zat.resize(NumCenters);
  Zspec.resize(NumSpecies);
  NofSpecies.resize(NumSpecies);
  for(int spec=0; spec<NumSpecies; spec++)
  {
    Zspec[spec] = tspecies(ChargeAttribIndx,spec);
    NofSpecies[spec] = static_cast<int>(tspecies(MemberAttribIndx,spec));
  }
  SpeciesID.resize(NumCenters);
  for(int iat=0; iat<NumCenters; iat++)
  {
    SpeciesID[iat]=P.GroupID[iat];
    Zat[iat] = Zspec[P.GroupID[iat]];
  }
  AA = LRCoulombSingleton::getHandler(P);
  //AA->initBreakup(*PtclRef);
  myConst=evalConsts();
  myRcut=AA->get_rc();//Basis.get_rc();
  if(rVs==0)
  {
    rVs = LRCoulombSingleton::createSpline4RbyVs(AA,myRcut,myGrid);
  }
}

CoulombPBCAA::Return_t
CoulombPBCAA::evalSR(ParticleSet& P)
{
  const DistanceTableData &d_aa(*P.DistTables[0]);
  RealType SR=0.0;
  for(int ipart=0; ipart<NumCenters; ipart++)
  {
    RealType esum = 0.0;
    for(int nn=d_aa.M[ipart],jpart=ipart+1; nn<d_aa.M[ipart+1]; nn++,jpart++)
    {
      //if(d_aa->r(nn)>=myRcut) continue;
      //esum += Zat[jpart]*AA->evaluate(d_aa->r(nn),d_aa->rinv(nn));
      esum += Zat[jpart]*d_aa.rinv(nn)*rVs->splint(d_aa.r(nn));
    }
    //Accumulate pair sums...species charge for atom i.
    SR += Zat[ipart]*esum;
  }
  return SR;
}

CoulombPBCAA::Return_t
CoulombPBCAA::evalLR(ParticleSet& P)
{
  const int slab_dir=OHMMS_DIM-1;
  RealType res=0.0;
  const StructFact& PtclRhoK(*(P.SK));
  if(PtclRhoK.SuperCellEnum==SUPERCELL_SLAB)
  {
    const DistanceTableData &d_aa(*P.DistTables[0]);
    //distance table handles jat>iat
    for(int iat=0; iat<NumCenters; ++iat)
    {
      RealType u=0;
#if !defined(USE_REAL_STRUCT_FACTOR)
      for(int nn=d_aa.M[iat], jat=iat+1; nn<d_aa.M[iat+1]; ++nn,++jat)
        u += Zat[jat]*AA->evaluate_slab(d_aa.dr(nn)[slab_dir]
                                        , PtclRhoK.KLists.kshell , PtclRhoK.eikr[iat], PtclRhoK.eikr[jat]);
#endif
      res += Zat[iat]*u;
    }
  }
  else
  {
    for(int spec1=0; spec1<NumSpecies; spec1++)
    {
      RealType Z1 = Zspec[spec1];
      for(int spec2=spec1; spec2<NumSpecies; spec2++)
      {
#if defined(USE_REAL_STRUCT_FACTOR)
        RealType temp=AA->evaluate(PtclRhoK.KLists.kshell
                                   , PtclRhoK.rhok_r[spec1], PtclRhoK.rhok_i[spec1]
                                   , PtclRhoK.rhok_r[spec2], PtclRhoK.rhok_i[spec2]);
#else
        RealType temp=AA->evaluate(PtclRhoK.KLists.kshell, PtclRhoK.rhok[spec1], PtclRhoK.rhok[spec2]);
#endif
        if(spec2==spec1)
          temp*=0.5;
        res += Z1*Zspec[spec2]*temp;
      } //spec2
    }//spec1
  }
  return res;
}


CoulombPBCAA::Return_t
CoulombPBCAA::evalLRwithForces(ParticleSet& P)
{
  RealType LR=0.0;
  const StructFact& PtclRhoK(*(P.SK));
  vector<TinyVector<RealType,DIM> > grad(P.getTotalNum());
  for(int spec2=0; spec2<NumSpecies; spec2++)
  {
    RealType Z2 = Zspec[spec2];
    for (int iat=0; iat<grad.size(); iat++)
      grad[iat] = TinyVector<RealType,DIM>(0.0);
    AA->evaluateGrad(P, P, spec2, Zat, grad);
    for (int iat=0; iat<grad.size(); iat++)
      forces[iat] += Z2*grad[iat];
  } //spec2
  return evalLR(P);
}



CoulombPBCAA::Return_t
CoulombPBCAA::evalSRwithForces(ParticleSet& P)
{
  const DistanceTableData *d_aa = P.DistTables[0];
  RealType SR=0.0;
  for(int ipart=0; ipart<NumCenters; ipart++)
  {
    RealType esum = 0.0;
    for(int nn=d_aa->M[ipart],jpart=ipart+1; nn<d_aa->M[ipart+1]; nn++,jpart++)
    {
      RealType rV, d_rV_dr, d2_rV_dr2;
      rV = rVs->splint(d_aa->r(nn), d_rV_dr, d2_rV_dr2);
      RealType V = rV *d_aa->rinv(nn);
      esum += Zat[jpart]*d_aa->rinv(nn)*rV;
      PosType grad = Zat[jpart]*Zat[ipart]*
                     (d_rV_dr - V)*d_aa->rinv(nn)*d_aa->rinv(nn)*d_aa->dr(nn);
      forces[ipart] += grad;
      forces[jpart] -= grad;
    }
    //Accumulate pair sums...species charge for atom i.
    SR += Zat[ipart]*esum;
  }
  return SR;
}


/** evaluate the constant term that does not depend on the position
 *
 * \htmlonly
 * <ul>
 * <li> self-energy: \f$ -\frac{1}{2}\sum_{i} v_l(r=0) q_i^2 = -\frac{1}{2}v_l(r=0) \sum_{alpha} N^{\alpha} q^{\alpha}^2\f$
 * <li> background term \f$ V_{bg} = -\frac{1}{2}\sum_{\alpha}\sum_{\beta} N^{\alpha}q^{\alpha}N^{\beta}q^{\beta} v_s(k=0)\f$
 * </ul>
 * \endhtmlonly
 * CoulombPBCABTemp contributes additional background term which completes the background term
 */
CoulombPBCAA::Return_t
CoulombPBCAA::evalConsts()
{
  //LRHandlerType::BreakupBasisType &Basis(AA->Basis);
  //const Vector<RealType> &coefs(AA->coefs);
  RealType Consts=0.0; // constant term
  //v_l(r=0) including correction due to the non-periodic direction
  RealType vl_r0 = AA->evaluateLR_r0();
  for(int spec=0; spec<NumSpecies; spec++)
  {
    RealType z = Zspec[spec];
    RealType n = NofSpecies[spec];
    Consts -= 0.5*vl_r0*z*z*n;
  }
  app_log() << "   PBCAA self-interaction term " << Consts << endl;
  //Compute Madelung constant: this is not correct for general cases
  MC0 = 0.0;
  for(int i=0; i<AA->Fk.size(); i++)
    MC0 += AA->Fk[i];
  MC0 = 0.5*(MC0 - vl_r0);
  //Neutraling background term
  RealType vs_k0=AA->evaluateSR_k0(); //v_s(k=0)
  for(int speca=0; speca<NumSpecies; speca++)
  {
    RealType za = Zspec[speca];
    RealType na = NofSpecies[speca];
    Consts -= 0.5*vs_k0*za*na*za*na;
    for(int specb=speca+1; specb<NumSpecies; specb++)
    {
      RealType zb = Zspec[specb];
      int nb = NofSpecies[specb];
      Consts -= vs_k0*za*zb*na*nb;
    }
  }
  app_log() << "   PBCAA total constant " << Consts << endl;
  //app_log() << "   MC0 of PBCAA " << MC0 << endl;
  return Consts;
}

QMCHamiltonianBase* CoulombPBCAA::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  if(is_active)
    return new CoulombPBCAA(qp,is_active,ComputeForces);
  else
    return new CoulombPBCAA(*this);//nothing needs to be re-evaluated
}
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 5915 $   $Date: 2013-08-05 08:11:42 -0500 (Mon, 05 Aug 2013) $
 * $Id: CoulombPBCAA.cpp 5915 2013-08-05 13:11:42Z jnkim $
 ***************************************************************************/

