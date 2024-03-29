//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
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
#include "QMCWaveFunctions/AtomicOrbitals/HFAtomicSTOSetBuilder.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCWaveFunctions/SlaterDeterminant.h"
#include "QMCWaveFunctions/MultiSlaterDeterminant.h"
#include "Utilities/OhmmsInfo.h"
#include <set>
#include <sstream>
using std::set;

namespace qmcplusplus
{

template<class T, class POS>
std::ostream&
operator<<(std::ostream& out, const ComboSTO<T,POS>& asto)
{
  out << asto.Name << " LM = " << asto.LM << endl;
  for(int i=0; i<asto.C.size(); i++)
    out << asto.Rnl[i]->ID << " " << asto.C[i] << endl;
  return out;
}

HFAtomicSTOSetBuilder::HFAtomicSTOSetBuilder(ParticleSet& els, TrialWaveFunction& wfs, ParticleSet& ions):
  OrbitalBuilderBase(els,wfs), Lmax(-1)
{
  d_table = DistanceTable::add(ions,els);
}

bool
HFAtomicSTOSetBuilder::getBasis(xmlNodePtr cur)
{
  typedef TrialWaveFunction::RealType RealType;
  //STONorm<RealType> anorm(4);
  cur = cur->xmlChildrenNode;
  while(cur!=NULL)
  {
    string cname((const char*)(cur->name));
    if(cname == basisfunc_tag)
    {
      XMLReport("Found a radial basis function.")
      int n=atoi((const char*)(xmlGetProp(cur, (const xmlChar *)"n")));
      int l=atoi((const char*)(xmlGetProp(cur, (const xmlChar *)"l")));
      Lmax = max(Lmax,l);
      string newRnl;
      const xmlChar* aptr(xmlGetProp(cur, (const xmlChar *)"id"));
      if(aptr)
      {
        newRnl = (const char*)aptr;
      }
      else
      {
        std::ostringstream idstream(newRnl.c_str());
        idstream << 'R' << RnlID.size();
        XMLReport("WARNING!!!! Rnl does not have id. Assign a default" << newRnl)
      }
      double screen = 1.0;
      xmlNodePtr s = cur->xmlChildrenNode;
      while(s != NULL)
      {
        string vname((const char*)(s->name));
        if(vname == param_tag)
        {
          putContent(screen,s);
        }
        s=s->next;
      }
      map<string,int>::iterator it = RnlID.find(newRnl);
      if(it == RnlID.end())
      {
        //the ID of the next Radial Function
        int id = Rnl.size();
        RnlID[newRnl] = id;
        Rnl.push_back(new RadialOrbital_t(n,l,screen));
        //Rnl.push_back(new RadialOrbital_t(n-l-1,screen,anorm(n-1,screen)));
        XMLReport("STO function (n,l,screen) " << n << " "  << l << " " << screen)
        Rnl[id]->ID = id;
      }
    }
    cur=cur->next;
  }
  return true;
}

HFAtomicSTOSet*
HFAtomicSTOSetBuilder::getOrbital(xmlNodePtr cur)
{
  typedef TrialWaveFunction::RealType RealType;
  HFAtomicSTOSet* psi = new HFAtomicSTOSet(Lmax);
  XMLReport("The maximum angular momentum " << Lmax)
  //set of  RadialFunctions(integers)
  set<int> RnlSet;
  int spin = atoi((const char*)(xmlGetProp(cur, (const xmlChar *)"spin")));
  xmlNodePtr orb=cur->xmlChildrenNode;
  while(orb!=NULL)
  {
    string cname((const char*)(orb->name));
    //a single-particle orbtial consiting of a number of Orbitals is added
    if(cname == spo_tag)
    {
      XMLReport("Found a SingleParticleOrbital.")
      //using an existing Orbtial
      xmlChar *att=xmlGetProp(orb,(const xmlChar*)"ref");
      if(att)
      {
        string aname((const char*)att);
        XMLReport("Using an existing Orbital id =" << aname)
        map<string,SPO_t*>::iterator it = OrbSet.find(aname);
        if(it == OrbSet.end())
        {
          ERRORMSG("Orbital " << aname << " does not exisit. Failed.")
        }
        else
        {
          XMLReport("Using an existing " << aname << " Orbtial")
          //create a ComboSTO using a copy Constructor
          int lm = (*it).second->LM;
          if(xmlHasProp(orb, (const xmlChar *)"l"))
          {
            int l=atoi((const char*)(xmlGetProp(orb, (const xmlChar *)"l")));
            int m=atoi((const char*)(xmlGetProp(orb, (const xmlChar *)"m")));
            lm = psi->Ylm.index(l,m);
          }
          SPO_t* aSTO = new SPO_t(lm, psi->Ylm,(*it).second->Rnl, &((*it).second->C[0]));
          psi->Orbital.push_back(aSTO);
          for(int nl=0; nl<aSTO->Rnl.size(); nl++)
            RnlSet.insert(aSTO->Rnl[nl]->ID);
        }
      }
      else
        // add a new orbital
      {
        //get the properties first
        int l=atoi((const char*)(xmlGetProp(orb, (const xmlChar *)"l")));
        int m=atoi((const char*)(xmlGetProp(orb, (const xmlChar *)"m")));
        string aname((const char*)(xmlGetProp(orb, (const xmlChar *)"id")));
        XMLReport("Adding a new Orbital id =" << aname << " l= " << l << " m=" << m)
        //temporary storage for coefficients and radial orbtials
        vector<RealType> c;
        vector<RadialOrbital_t*> sto;
        RealType ctmp;
        xmlNodePtr s = orb->xmlChildrenNode;
        while(s != NULL)
        {
          string sname((const char*)(s->name));
          if(sname == basisfunc_tag)
          {
            //getting the index of the radial basis functions
            string  rfn((const char*)(xmlGetProp(s, (const xmlChar *)"ref")));
            //c.push_back(atof((const char*)(xmlGetProp(s, (const xmlChar *)"C"))));
            xmlNodePtr z=s->xmlChildrenNode;
            while(z != NULL)
            {
              string vname((const char*)(z->name));
              if(vname == param_tag)
              {
                putContent(ctmp,z);
                c.push_back(ctmp);
              }
              z=z->next;
            }
            XMLReport("Adding STO " << rfn << " C = " << c.back())
            map<string,int>::iterator it = RnlID.find(rfn);
            sto.push_back(Rnl[(*it).second]);
            RnlSet.insert((*it).second);
          }
          s = s->next;
        }
        ///create a ComboSTO
        SPO_t* aSTO = new SPO_t(psi->Ylm.index(l,m), psi->Ylm, sto, &c[0]);
        aSTO->Name = aname;
        ///add the ComboSTO to HFAtomicSTOSet::Orbital
        psi->Orbital.push_back(aSTO);
        ///register name for re-use
        OrbSet[aname] = aSTO;
      }///
    }///found Orbtial
    orb = orb->next;
  }
  ///time to add a set of Rnl to HFAtomicSTOSet for pre-calculations
  set<int>::iterator irnl = RnlSet.begin();
  while(irnl != RnlSet.end())
  {
    psi->RnlPool.push_back(Rnl[*irnl]);
    irnl++;
  }
  return psi;
}

bool
HFAtomicSTOSetBuilder::put(xmlNodePtr cur)
{
  typedef TrialWaveFunction::RealType RealType;
  typedef DiracDeterminant<HFAtomicSTOSet> Det_t;
  typedef SlaterDeterminant<HFAtomicSTOSet> SlaterDeterminant_t;
  SlaterDeterminant_t *asymmpsi;
  vector<SlaterDeterminant_t*> slaterdets;
  vector<RealType> C;
  int is=0;
  cur = cur->xmlChildrenNode;
  while(cur != NULL)
  {
    string cname((const char*)(cur->name));
    if(cname == basisset_tag)
    {
      getBasis(cur);
    }
    else
      if(cname == sd_tag)
      {
        int first = 0;
        xmlNodePtr tcur = cur->xmlChildrenNode;
        slaterdets.push_back(new SlaterDeterminant_t);
        C.push_back(1.0);
        //if(xmlHasProp(tcur,(const xmlChar*)"C")) {
        //C[is] = atof((const char*)(xmlGetProp(tcur, (const xmlChar *)"C")));
        //}
        while(tcur != NULL)
        {
          string vname((const char*)(tcur->name));
          if(vname == param_tag)
          {
            putContent(C[is],tcur);
          }
          else
            if(vname == det_tag)
            {
              HFAtomicSTOSet* psi = getOrbital(tcur);
              psi->setTable(d_table);
              Det_t * adet = new Det_t(*psi,first);
              adet->set(first,psi->size());
              XMLReport("Adding a determinant to the SlaterDeterminant " << first<< " " << psi->size())
              slaterdets[is]->add(adet);
              first += psi->size();
            }
          tcur = tcur->next;
        }
        is++;
      }
    cur = cur->next;
  }
  XMLReport("Done with he initialization HFAtomicSTSet")
  if(slaterdets.size() > 1)
  {
    MultiSlaterDeterminant<HFAtomicSTOSet> *multidet=
      new MultiSlaterDeterminant<HFAtomicSTOSet>;
    HFAtomicSTOSet::BasisSet_t *bs=new HFAtomicSTOSet::BasisSet_t;
    for(int i=0; i<slaterdets.size(); i++)
    {
      cout << "Multi determinant " << C[i] << endl;
      slaterdets[i]->setBasisSet(bs);
      multidet->add(slaterdets[i],C[i]);
    }
    targetPsi.addOrbital(multidet);
  }
  else
  {
    slaterdets[0]->setBasisSet(new HFAtomicSTOSet::BasisSet_t);
    targetPsi.addOrbital(slaterdets[0]);
  }
  return true;
}
}
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 19:14:53 -0500 (Thu, 25 Apr 2013) $
 * $Id: HFAtomicSTOSetBuilder.cpp 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/
