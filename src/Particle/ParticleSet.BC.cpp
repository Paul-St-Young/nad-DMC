//////////////////////////////////////////////////////////////////
// (c) Copyright 2004- by Jeongnim Kim
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
/**@file ParticleSet.BC.cpp
 * @brief definition of functions controlling Boundary Conditions
 */
#include "Particle/ParticleSet.h"
#include "Particle/FastParticleOperators.h"
#include "Utilities/OhmmsInfo.h"
#include "Message/OpenMP.h"
#include "LongRange/StructFact.h"

namespace qmcplusplus
{

/** Creating StructureFactor
 *
 * Currently testing only 1 component for PBCs.
 */
void ParticleSet::createSK()
{
  //if(!sorted_ids && !reordered_ids)
  //{
  //  //save ID and GroupID
  //  orgID=ID;
  //  orgGroupID=GroupID;
  //  if(groups()<1)
  //  {
  //    int nspecies=mySpecies.getTotalNum();
  //    vector<int> ppg(nspecies,0);
  //    for(int iat=0; iat<GroupID.size(); ++iat) ppg[GroupID[iat]]+=1;
  //    SubPtcl.resize(nspecies+1);
  //    SubPtcl[0]=0;
  //    for(int i=0; i<nspecies; ++i) SubPtcl[i+1]=SubPtcl[i]+ppg[i];
  //    int new_id=0;
  //    for(int i=0; i<nspecies; ++i)
  //      for(int iat=0; iat<GroupID.size(); ++iat) if(GroupID[iat]==i) orgID[new_id++]=ID[iat];
  //    bool grouped=true;
  //    for(int iat=0; iat<ID.size(); ++iat) grouped &= (orgID[iat]==ID[iat]);
  //    if(grouped)
  //    {
  //      app_log() << "  ParticleSet is grouped. No need to reorder." << endl;
  //    }
  //    else
  //    {
  //      app_log() << "  Need to reorder. Only R is swapped." << endl;
  //      ParticlePos_t oldR(R);
  //      for(int iat=0; iat<R.size(); ++iat) R[iat]=oldR[orgID[iat]];
  //      for(int i=0; i<groups(); ++i)
  //        for(int iat=first(i); iat<last(i); ++iat) GroupID[iat]=i;
  //      reordered_ids=true;
  //    }
  //  }//once group is set, nothing to be done
  //  sorted_ids=true;
  //}
  //int membersize= mySpecies.addAttribute("membersize");
  //for(int ig=0; ig<mySpecies.size(); ++ig)
  //  SubPtcl[ig+1]=SubPtcl[ig]+mySpecies(membersize,ig);
  convert2Cart(R); //make sure that R is in Cartesian coordinates
  //if(Lattice.BoxBConds[0] && SK == 0)
  if(Lattice.SuperCellEnum != SUPERCELL_OPEN)
  {
    Lattice.SetLRCutoffs();
    if(SK)
    {
      app_log() << "\n  Structure Factor is reset by " << Lattice.LR_kc << endl;
      SK->UpdateNewCell(Lattice.LR_kc);
    }
    else
    {
      app_log() << "\n  Creating Structure Factor for periodic systems." <<endl;
      SK = new StructFact(*this,Lattice.LR_kc);
    }
    //Lattice.print(app_log());
    //This uses the copy constructor to avoid recomputing the data.
    //SKOld = new StructFact(*SK);
  }
  //set the mass array
  int beforemass=mySpecies.numAttributes();
  int massind= mySpecies.addAttribute("mass");
  if(beforemass == massind)
  {
    app_log() << "  ParticleSet::createSK setting mass of  " << getName() << " to 1.0" << endl;
    for(int ig=0; ig<mySpecies.getTotalNum(); ++ig)
      mySpecies(massind,ig)=1.0;
  }
  for(int iat=0; iat<GroupID.size(); iat++)
    Mass[iat]=mySpecies(massind,GroupID[iat]);
}

void ParticleSet::convert(const ParticlePos_t& pin, ParticlePos_t& pout)
{
  if(pin.getUnit() == pout.getUnit())
  {
    pout = pin;
    return;
  }
  if(pin.getUnit() == PosUnit::LatticeUnit)
    //convert to CartesianUnit
  {
    ConvertPosUnit<ParticlePos_t,Tensor_t,DIM,OHMMS_ORTHO>::apply(pin,Lattice.R,pout,0,pin.size());
  }
  else
    //convert to LatticeUnit
  {
    ConvertPosUnit<ParticlePos_t,Tensor_t,DIM,OHMMS_ORTHO>::apply(pin,Lattice.G,pout,0,pin.size());
  }
}

void ParticleSet::convert2Unit(const ParticlePos_t& pin, ParticlePos_t& pout)
{
  pout.setUnit(PosUnit::LatticeUnit);
  if(pin.getUnit() == PosUnit::LatticeUnit)
    pout = pin;
  else
    ConvertPosUnit<ParticlePos_t,Tensor_t,DIM,OHMMS_ORTHO>::apply(pin,Lattice.G,pout,0,pin.size());
}

void ParticleSet::convert2Cart(const ParticlePos_t& pin, ParticlePos_t& pout)
{
  pout.setUnit(PosUnit::CartesianUnit);
  if(pin.getUnit() == PosUnit::CartesianUnit)
    pout = pin;
  else
    ConvertPosUnit<ParticlePos_t,Tensor_t,DIM,OHMMS_ORTHO>::apply(pin,Lattice.R,pout,0,pin.size());
}

void ParticleSet::convert2Unit(ParticlePos_t& pinout)
{
  if(pinout.getUnit() == PosUnit::LatticeUnit)
    return;
  else
  {
    pinout.setUnit(PosUnit::LatticeUnit);
    ConvertPosUnit<ParticlePos_t,Tensor_t,DIM,OHMMS_ORTHO>::apply(pinout,Lattice.G,0,pinout.size());
  }
}

void ParticleSet::convert2Cart(ParticlePos_t& pinout)
{
  if(pinout.getUnit() == PosUnit::CartesianUnit)
    return;
  else
  {
    pinout.setUnit(PosUnit::CartesianUnit);
    ConvertPosUnit<ParticlePos_t,Tensor_t,DIM,OHMMS_ORTHO>::apply(pinout,Lattice.R,0,pinout.size());
  }
}

void ParticleSet::applyBC(const ParticlePos_t& pin, ParticlePos_t& pout)
{
  applyBC(pin,pout,0,pin.size());
}

void ParticleSet::applyBC(const ParticlePos_t& pin, ParticlePos_t& pout, int first, int last)
{
  const bool orthogonal = ParticleLayout_t::IsOrthogonal;
  int mode = pin.getUnit()*2+pout.getUnit();
  switch(mode)
  {
  case(0):
    ApplyBConds<ParticlePos_t,Tensor_t,DIM,orthogonal>::Cart2Cart(pin,Lattice.G,Lattice.R,pout,first,last);
    break;
  case(1):
    ApplyBConds<ParticlePos_t,Tensor_t,DIM,orthogonal>::Cart2Unit(pin,Lattice.G,pout,first,last);
    break;
  case(2):
    ApplyBConds<ParticlePos_t,Tensor_t,DIM,orthogonal>::Unit2Cart(pin,Lattice.R,pout,first,last);
    break;
  case(3):
    ApplyBConds<ParticlePos_t,Tensor_t,DIM,orthogonal>::Unit2Unit(pin,pout,first,last);
    break;
  }
}

void ParticleSet::applyBC(ParticlePos_t& pos)
{
  const bool orthogonal = ParticleLayout_t::IsOrthogonal;
  if(pos.getUnit()==PosUnit::LatticeUnit)
  {
    ApplyBConds<ParticlePos_t,Tensor_t,DIM,orthogonal>::Unit2Unit(pos,0,LocalNum);
  }
  else
  {
    ApplyBConds<ParticlePos_t,Tensor_t,DIM,orthogonal>::Cart2Cart(pos,Lattice.G,Lattice.R,0,LocalNum);
  }
}

void ParticleSet::applyMinimumImage(ParticlePos_t& pinout)
{
  if(Lattice.SuperCellEnum==SUPERCELL_OPEN)
    return;
  for(int i=0; i<pinout.size(); ++i)
    MinimumImageBConds<RealType,DIM>::apply(Lattice.R,Lattice.G,pinout[i]);
}

void ParticleSet::convert2UnitInBox(const ParticlePos_t& pin, ParticlePos_t& pout)
{
  pout.setUnit(PosUnit::LatticeUnit);
  convert2Unit(pin,pout); // convert to crystalline unit
  put2box(pout);
}

void ParticleSet::convert2CartInBox(const ParticlePos_t& pin, ParticlePos_t& pout)
{
  convert2UnitInBox(pin,pout); // convert to crystalline unit
  convert2Cart(pout);
}
}
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 19:14:53 -0500 (Thu, 25 Apr 2013) $
 * $Id: ParticleSet.BC.cpp 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/

