//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
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
#ifndef QMCPLUSPLUS_MODIFIED_PADEFUNCTION_H
#define QMCPLUSPLUS_MODIFIED_PADEFUNCTION_H

#include "Numerics/OptimizableFunctorBase.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{
/** ModPade Jastrow functional
 *
 * \f[ u(r) = \frac{1}{2*A}\left[1-\exp(-A*r)\right], \f]
 */
template<class T>
struct ModPadeFunctor: public OptimizableFunctorBase
{

  ///coefficients
  real_type A;
  real_type B;
  real_type Zeff;

  real_type Coeff;
  real_type mAB;

  string ID_B;

  /** constructor
   * @param a A coefficient
   * @param samespin boolean to indicate if this function is for parallel spins
   */
  ModPadeFunctor(real_type a=-0.5, real_type b=1): A(a),B(b),Zeff(1.0)
  {
    reset();
  }

  OptimizableFunctorBase* makeClone() const
  {
    return new ModPadeFunctor<T>(*this);
  }

  /** reset the internal variables.
   *@param a New Jastrow parameter a
   */
  inline void reset()
  {
    Coeff=A/B;
    mAB = -A*B;
  }

  /** evaluate the value at r
   * @param r the distance
   * @return \f$ u(r) = \frac{A}{r}\left[1-\exp(-\frac{r}{F})\right]\f$
   */
  inline real_type evaluate(real_type r)
  {
    return Coeff*(1.0-std::exp(-B*r));
  }

  /**@param r the distance
    @param dudr first derivative
    @param d2udr second derivative
    @return the value
    */
  inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2)
  {
    real_type expmar=std::exp(-B*r);
    dudr=A*expmar;
    d2udr2=mAB*expmar;
    return Coeff*(1.0-expmar);
  }

  /**@param r the distance
    @param dudr first derivative
    @param d2udr second derivative
    @return the value
    */
  inline real_type evaluate(real_type r, real_type& dudr,
                            real_type& d2udr2, real_type d3udr3)
  {
    std::cerr << "Third derivative not implemented for ModPadeFunctor.\n";
    real_type expmar=std::exp(-B*r);
    dudr=A*expmar;
    d2udr2=mAB*expmar;
    return Coeff*(1.0-expmar);
  }


  /** return a value at r
  */
  real_type f(real_type r)
  {
    return evaluate(r);
  }

  /** return a derivative at r
  */
  real_type df(real_type r)
  {
    return A*std::exp(-B*r);
  }

  /** Read in the parameter from the xml input file.
   * @param cur current xmlNode from which the data members are reset
   */
  bool put(xmlNodePtr cur)
  {
    Zeff=1.0;
    //jastrow[iab]->put(cur->xmlChildrenNode,wfs_ref.RealVars);
    cur = cur->xmlChildrenNode;
    while(cur != NULL)
    {
      //@todo Var -> <param(eter) role="opt"/>
      string cname((const char*)(cur->name));
      if(cname == "parameter" || cname == "Var")
      {
        string aname("B"),idname("0");
        OhmmsAttributeSet p;
        p.add(aname,"name");
        p.add(idname,"id");
        p.put(cur);
        if(aname == "A")
        {
          putContent(A,cur);
        }
        else
          if(aname == "B")
          {
            ID_B = idname;
            putContent(B,cur);
          }
          else
            if(aname == "Z")
            {
              putContent(Zeff,cur);
            }
      }
      cur = cur->next;
    }
    reset();
    myVars.insert(ID_B,B);
    return true;
  }

  void checkInVariables(opt_variables_type& active)
  {
    active.insertFrom(myVars);
  }

  void checkOutVariables(const opt_variables_type& active)
  {
    myVars.getIndex(active);
  }

  inline void resetParameters(const opt_variables_type& active)
  {
    int loc=myVars.where(0);
    if(loc>-1)
      B = myVars[0] = active[loc];
    Coeff=A/B;
    mAB = -A*B;
  }

};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 19:14:53 -0500 (Thu, 25 Apr 2013) $
 * $Id: ModPadeFunctor.h 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/

