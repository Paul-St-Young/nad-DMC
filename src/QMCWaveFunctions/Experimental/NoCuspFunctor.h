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
/** @file NoCuspFunctor.h
 * @brief Define NoCupsFunctor \f$ u(r) = \frac{A}{1+Br^2} \f$
 */
#ifndef QMCPLUSPLUS_NOCUSPJASTROW_H
#define QMCPLUSPLUS_NOCUSPJASTROW_H
#include "Numerics/OptimizableFunctorBase.h"

namespace qmcplusplus
{
/** NoCusp functor \f$ u(r) = \frac{A}{1+Br^2} \f$
 *
 * Prototype of the template parameter of TwoBodyJastrow and OneBodyJastrow
 */
template<class T>
struct NoCuspFunctor: public OptimizableFunctorBase<T>
{

  typedef typename OptimizableFunctorBase<T>::real_type real_type;
  typedef typename OptimizableFunctorBase<T>::OptimizableSetType OptimizableSetType;
  ///coefficients
  real_type A, B, AB2;
  ///id of A
  string ID_A;
  ///id of B
  string ID_B;
  ///constructor
  NoCuspFunctor(real_type a=1.0, real_type b=1.0)
  {
    reset(a,b);
  }

  OptimizableFunctorBase<T>* makeClone() const
  {
    return new NoCuspFunctor<T>(*this);
  }

  /**
   *@param a NoCusp Jastrow parameter a
   *@param b NoCusp Jastrow parameter b
   *@brief reset the internal variables.
   */
  void reset(real_type a, real_type b)
  {
    A=a;
    B=b;
    AB2=2.0*a*b;
  }

  /**@param r the distance
    @return \f$ u(r) = a/(1+br^2) \f$
    */
  inline real_type evaluate(real_type r)
  {
    return A/(1.0+B*r*r);
  }

  /**@param r the distance
    @param dudr return value  \f$ du/dr = -2abr/(1+br^2)^2 \f$
    @param d2udr2 return value  \f$ d^2u/dr^2 =
    -2ab(1-3br^2)/(1+br^2)^3 \f$
    @return \f$ u(r) = a/(1+br^2) \f$
    */
  inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2)
  {
    real_type rsq = r*r;
    real_type u = 1.0/(1.0+B*rsq);
    real_type usq = u*u;
    dudr = -AB2*r*usq;
    d2udr2 = -AB2*usq*u*(1.0-3.0*B*rsq);
    return A*u;
  }

  real_type f(real_type r)
  {
    return evaluate(r);
  }

  real_type df(real_type r)
  {
    real_type dudr,d2udr2;
    real_type res=evaluate(r,dudr,d2udr2);
    return dudr;
  }

  /** implements virtual function
   * @param cur xml node
   */
  bool put(xmlNodePtr cur)
  {
    real_type Atemp,Btemp;
    //jastrow[iab]->put(cur->xmlChildrenNode,wfs_ref.RealVars);
    xmlNodePtr tcur = cur->xmlChildrenNode;
    while(tcur != NULL)
    {
      string cname((const char*)(tcur->name));
      if(cname == "parameter" || cname == "Var")
      {
        string aname((const char*)(xmlGetProp(tcur,(const xmlChar *)"name")));
        string idname((const char*)(xmlGetProp(tcur,(const xmlChar *)"id")));
        if(aname == "A")
        {
          ID_A = idname;
          putContent(Atemp,tcur);
        }
        else
          if(aname == "B")
          {
            ID_B = idname;
            putContent(Btemp,tcur);
          }
      }
      tcur = tcur->next;
    }
    reset(Atemp,Btemp);
    LOGMSG("  NoCuspFunctor Parameters ")
    LOGMSG("    A (" << ID_A << ") = " << A  << "  B (" << ID_B << ") =  " << B)
    return true;
  }

  void addOptimizables(OptimizableSetType& vlist)
  {
    vlist[ID_A]=A;
    vlist[ID_B]=B;
  }

  /** reset the internal variables.
   *
   * USE_resetParameters
   */
  void resetParameters(OptimizableSetType& optVariables)
  {
    typename OptimizableSetType::iterator it_a(optVariables.find(ID_A));
    if(it_a != optVariables.end())
      A=(*it_a).second;
    typename OptimizableSetType::iterator it_b(optVariables.find(ID_B));
    if(it_b != optVariables.end())
      B=(*it_b).second;
    AB2 = 2.0*A*B;
  }
};

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 19:14:53 -0500 (Thu, 25 Apr 2013) $
 * $Id: NoCuspFunctor.h 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/

