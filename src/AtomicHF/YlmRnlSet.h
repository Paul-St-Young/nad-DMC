//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim and Jordan Vincent
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
/*@author Jeongnim Kim and Jordan Vincent
 */
#ifndef OHMMS_YLMRNLSET_ONGRID_H
#define OHMMS_YLMRNLSET_ONGRID_H
#include <map>
#include <set>
#include <vector>
#include <stdlib.h>
#include "OhmmsData/OhmmsElementBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"
#include "Numerics/OneDimIntegration.h"
#include "OhmmsPETE/TinyVector.h"


struct ltnlm
{
  bool operator()(const TinyVector<int,3>& a,
                  const TinyVector<int,3>& b) const
  {
    if(a[0] > b[0])
      return false;
    else
      if(a[0] == b[0])
      {
        //   return a[1] < b[1];
        if (a[1] < b[1])
          return true;
        else
          if(a[1] == b[1])
            return a[2] < b[2];
          else
            return true;
      }
      else
      {
        return true;
      }
  }
};

struct ltnl
{
  bool operator()(const TinyVector<int,2>& a,
                  const TinyVector<int,2>& b) const
  {
    if(a[0] > b[0])
      return false;
    else
      if(a[0] == b[0])
      {
        return a[1] < b[1];
      }
      else
      {
        return true;
      }
  }
};

template<unsigned D>
struct equal_nlms
{

  bool operator()(const TinyVector<int,D>& a,
                  const TinyVector<int,D>& b) const
  {
    int i=0;
    while(i<D)
    {
      if(a[i]!= b[i])
        return true;
      i++;
    }
    return false;
  }
};


template<class T>
struct YlmRnlSet
{

  typedef T                                        value_type;
  typedef OneDimGridBase<value_type>               RadialGrid_t;
  typedef OneDimCubicSpline<value_type,value_type> RadialOrbital_t;
  typedef vector<RadialOrbital_t>                  RadialOrbitalSet_t;

  typedef TinyVector<int,2> NLIndex;
  typedef TinyVector<int,3> NLMIndex;
  typedef TinyVector<int,4> NLMSIndex;

  typedef std::map<NLIndex,int,ltnl >  NL_Map_t;
  typedef std::map<NLMIndex,int,ltnlm > NLM_Map_t;
  typedef std::map<NLMSIndex,int,equal_nlms<4> >  NLMS_Map_t;

  ///constructor
  YlmRnlSet(): m_grid(NULL), NumOrb(0),
    NumUniqueOrb(0), Restriction("none"), CuspParam(0.0),
    EigenBoundParam(0.0) { }

  bool add(int n, int l, int m, int s, value_type occ);

  void applyRestriction(int norb);

  ///normailize all the orbitals
  void normalize(int norb)
  {
    for(int i=0; i < norb; i++)
      normalize_RK2(psi[i]);
  }

  ///assigns orbital \f$ \psi_{iorb}(r) \f$
  inline RadialOrbital_t& operator()(int iorb)
  {
    return psi[iorb];
  }

  ///returns a reference to orbital \f$ \psi_{iorb}(r) \f$
  inline const RadialOrbital_t& operator()(int iorb) const
  {
    return psi[iorb];
  }

  ///return the number of orbitals
  inline int size() const
  {
    return psi.size();
  }

  ///return \f$ \psi_{iorb}(r_j) \f$
  inline value_type operator()(int iorb, int j) const
  {
    return psi[iorb](j);
  }

  ///assigns orbital \f$ \psi_{iorb}(r_j) \f$
  inline value_type& operator()(int iorb, int j)
  {
    return psi[iorb](j);
  }

  inline const value_type dR0(int iorb)
  {
    return psi[iorb].yprime0;
  }

  inline const value_type dRN(int iorb)
  {
    return psi[iorb].yprimeN;
  }

  ///reset the values of orbitals
  void reset()
  {
    for(int i=0; i < psi.size(); i++)
      psi[i].m_Y = 0.0;
  }

  ///restriction type;
  string Restriction;

  ///number of orbitals
  int NumOrb;

  ///number of unique orbitals
  int NumUniqueOrb;

  ///a common grid
  RadialGrid_t* m_grid;

  /**@defgroup QuantumNumber
   * Quantum numbers of an orbital.
   * The i-th orbital is represented by (N[i],L[i],M[i],S[i])
   *@{n = Principal quantum number
   */
  vector<int> N;

  /*@ l = angular momentum */
  vector<int> L;

  /*@ M = z-angular momentum */
  vector<int> M;

  /*@ S = spin in unit of 1/2 */
  vector<int> S;

  /**@}*///end of group QuantumNumber

  ///coefficient for each orbtial (1 if occupied, 0 if not)
  vector<value_type> Occ;

  ///parameter for calculating orbital cusp conditions
  value_type CuspParam;

  ///parameter for calculating orbital cusp conditions
  value_type EigenBoundParam;

  ///map for nl quantum numbers
  NL_Map_t NL;
  ///map for nlm quantum numbers
  NLM_Map_t NLM;
  ///map for occupation number
  NLMS_Map_t OccNo;

  ///assign an id for each orbital
  vector<int> ID;
  ///keeps track of the number of orbitals with the same id
  vector<int> IDcount;
  // vector<int> IDmap;
  ///the radial grid orbitals
  vector<RadialOrbital_t> psi;

  bool put(xmlNodePtr cur);
  bool get(std::ostream& os);
//  bool getqmc(std::ostream& os);
//  bool getsiesta(std::ostream& os, int orb);
  bool print(const std::string&);
};
#include "AtomicHF/YlmRnlSet.cpp"
#include "AtomicHF/YlmRnlSet.IO.h"
#endif

/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 19:14:53 -0500 (Thu, 25 Apr 2013) $
 * $Id: YlmRnlSet.h 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/

