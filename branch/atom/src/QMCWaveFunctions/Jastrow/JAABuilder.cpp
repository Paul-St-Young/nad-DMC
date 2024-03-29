//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
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
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCWaveFunctions/Jastrow/JAABuilder.h"
#include "QMCWaveFunctions/Jastrow/ModPadeFunctor.h"
#include "QMCWaveFunctions/Jastrow/PadeFunctors.h"
#include "QMCWaveFunctions/Jastrow/McMillanJ2Functor.h"
#include "QMCWaveFunctions/Jastrow/McMillanJ2GFunctor.h"
#include "QMCWaveFunctions/Jastrow/GaussianFunctor.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/DiffTwoBodyJastrowOrbital.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

JAABuilder::JAABuilder(ParticleSet& p, TrialWaveFunction& psi):OrbitalBuilderBase(p,psi),
  IgnoreSpin(true)
{
}
/** Create a two-body Jatrow function with a template
 *@param cur the current xmlNode
 *@param dummy null pointer used to identify FN
 *
 *The template class JeeType is a functor which handles the
 *evaluation of the function value, gradients and laplacians using
 *distance tables. This is a specialized builder function for
 *spin-dependent Jastrow function,e.g., for electrons, two functions
 *are created for uu(dd) and ud(du).
 */
template <class FN> TwoBodyJastrowOrbital<FN>* JAABuilder::createJAA(xmlNodePtr cur, const string& jname)
{
  string corr_tag("correlation");
  int ng = targetPtcl.groups();
  int ia=0, ib=0, iab=0;
  xmlNodePtr gridPtr=NULL;
  cur = cur->children;
  const SpeciesSet& species(targetPtcl.getSpeciesSet());
  typedef TwoBodyJastrowOrbital<FN> JeeType;
  int taskid=(targetPsi.is_manager())?targetPsi.getGroupID():-1;
  JeeType *J2 = new JeeType(targetPtcl,taskid);
  typedef DiffTwoBodyJastrowOrbital<FN> dJ2Type;
  dJ2Type *dJ2 = new dJ2Type(targetPtcl);
  RealType rc=targetPtcl.Lattice.WignerSeitzRadius;
  int pairs=0;
  while (cur != NULL)
  {
    string cname((const char*)(cur->name));
    if (cname == corr_tag)
    {
      string spA("u");
      string spB("u");
      OhmmsAttributeSet rAttrib;
      rAttrib.add(spA, "speciesA");
      rAttrib.add(spA, "species1");
      rAttrib.add(spB, "speciesB");
      rAttrib.add(spB, "species2");
      rAttrib.put(cur);
      if (spA==targetPsi.getName()) //could have used the particle name
      {
        spA=species.speciesName[0];
        spB=species.speciesName[0];
      }
      int ia = species.findSpecies(spA);
      int ib = species.findSpecies(spB);
      if (ia==species.size() || ia == species.size())
      {
        APP_ABORT("JAABuilder::createJAA is trying to use invalid species");
      }
      string pairID=spA+spB;
      FN *j= new FN;
      j->cutoff_radius=rc;
      j->put(cur);
      J2->addFunc( ia,ib,j);
      dJ2->addFunc( ia,ib,j);
      ++pairs;
    }
    cur = cur->next;
  } // while cur
  if (pairs)
  {
    J2->dPsi=dJ2;
    string j2name="J2_"+jname;
    targetPsi.addOrbital(J2,j2name);
    return J2;
  }
  else
  {
    //clean up and delete the twobody orbitals
    APP_ABORT("JAABuilder::put Failed to create Two-Body with "+jname);
    return 0;
  }
}

bool JAABuilder::put(xmlNodePtr cur)
{
  string spinOpt("no");
  string typeOpt("Two-Body");
  string jastfunction("pade");
  OhmmsAttributeSet aAttrib;
  aAttrib.add(spinOpt,"spin");
  aAttrib.add(typeOpt,"type");
  aAttrib.add(jastfunction,"function");
  aAttrib.put(cur);
  IgnoreSpin=(spinOpt=="no");
  OrbitalBase* newJ = 0;
  //if(jastfunction == "pade") {
  //  app_log() << "  Two-Body Jastrow Function = " << jastfunction << endl;
  //  PadeJastrow<RealType> *dummy = 0;
  //  success = createJAA(cur,dummy);
  //} else
  if (jastfunction == "short")
  {
    app_log() << "  Modified Jastrow function Two-Body Jastrow Function = " << jastfunction << endl;
    IgnoreSpin=true;
    //ModPadeFunctor<RealType> *dummy = 0;
    newJ = createJAA<ModPadeFunctor<RealType> >(cur,jastfunction);
  }
  else if (jastfunction == "modmcmillan")
  {
    app_log() << "  Modified McMillan Jastrow function Two-Body Jastrow Function = " << jastfunction << endl;
    IgnoreSpin=true;
    newJ = createJAA<ModMcMillanJ2Functor<RealType> >(cur,jastfunction);
  }
  else if (jastfunction == "combomcmillan")
  {
    app_log() << "  Combo McMillan Jastrow function Two-Body Jastrow Function = " << jastfunction << endl;
    IgnoreSpin=true;
    newJ = createJAA<comboMcMillanJ2Functor<RealType> >(cur,jastfunction);
  }
  else if (jastfunction == "mcmillan")
  {
    app_log() << "  McMillan (LONG RANGE!) Two-Body Jastrow Function = " << jastfunction << endl;
    IgnoreSpin=true;
    SpeciesSet& species(targetPtcl.getSpeciesSet());
    TwoBodyJastrowOrbital<McMillanJ2Functor<RealType> >* a = createJAA<McMillanJ2Functor<RealType> >(cur,jastfunction);
    species(species.addAttribute("J2_A"),species.addSpecies(species.speciesName[targetPtcl.GroupID[0]])) = (a->F[a->F.size()-1])->A;
    species(species.addAttribute("J2_B"),species.addSpecies(species.speciesName[targetPtcl.GroupID[0]])) = (a->F[a->F.size()-1])->B;
    newJ = a;
  }
  else if (jastfunction == "mcmillanj2g")
  {
    app_log() << "  McMillan Two-Body Jastrow Function (Gaussian for r < 2.5) = " << jastfunction << endl;
    IgnoreSpin=true;
    SpeciesSet& species(targetPtcl.getSpeciesSet());
    TwoBodyJastrowOrbital<McMillanJ2GFunctor<RealType> >* a = createJAA<McMillanJ2GFunctor<RealType> >(cur,jastfunction);
    species(species.addAttribute("J2_A"),species.addSpecies(species.speciesName[targetPtcl.GroupID[0]])) = (a->F[a->F.size()-1])->A;
    species(species.addAttribute("J2_B"),species.addSpecies(species.speciesName[targetPtcl.GroupID[0]])) = (a->F[a->F.size()-1])->B;
    newJ = a;
  }
  else if (jastfunction == "gaussian")
  {
    app_log() << "  Gaussian function Two-Body Jastrow Function = " << jastfunction << endl;
    IgnoreSpin=true;
    newJ = createJAA<GaussianFunctor<RealType> >(cur,jastfunction);
  }
  else if (jastfunction == "shiftedgaussian")
  {
    app_log() << "  Gaussian function Two-Body Jastrow Function = " << jastfunction << endl;
    IgnoreSpin=true;
    newJ = createJAA<TruncatedShiftedGaussianFunctor<RealType> >(cur,jastfunction);
  }
  else if (jastfunction == "padetwo2ndorderfunctor")
  {
    app_log() << "  PadeTwo2ndOrderFunctor Jastrow function Two-Body Jastrow Function = " << jastfunction << endl;
    //IgnoreSpin=true;
    newJ = createJAA<PadeTwo2ndOrderFunctor<RealType> >(cur,jastfunction);
  }
//} else if(jastfunction == "rpa") {
  //  app_log() << "  Two-Body Jastrow Function = " << jastfunction << endl;
  //  RPAJastrow<RealType> *dummy = 0;
  //  success = createJAA(cur,dummy);
  //}
  return (newJ != 0);
}
}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 5875 $   $Date: 2013-05-30 13:10:06 -0500 (Thu, 30 May 2013) $
 * $Id: JAABuilder.cpp 5875 2013-05-30 18:10:06Z jnkim $
 ***************************************************************************/
