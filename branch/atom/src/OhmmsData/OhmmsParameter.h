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
/** @file OhmmsParameter.h
 * @brief Declaration of OhmmsParameter class.
 */
#ifndef OHMMS_OHMMSPARAMETER_H
#define OHMMS_OHMMSPARAMETER_H

#include "OhmmsData/OhmmsElementBase.h"

/** generic class for parameter xmlNode
 *
 * <parameter/> node is used to generically add a named parameter whose
 *value is the content of an xmlNode. The definition confirms docbook::parameter
 *\htmlonly
 Usage is:
 &lt;parameter name="<b>aname</b>" condition="<b>unit</b>"&gt; <b>value</b> &lt;/parameter&gt;
 <br>
  Example is:
 &lt;parameter name="temperature" condition="K"&gt; 200 &lt;/parameter&gt;
 <br>
 <ul>
 <li> <b>name</b> is the name of the parameter.
 <li> <b>condition</b> is the unit of the parameter: default = none for dimensionless values
 <li> <b>value</b> is the content of the parameter of type T.
 </ul>
 <ul> Two kinds of template parameter T are valid:
 <li> intrinsic C/C++ variables, such as int,double
 <li> Ohmms basic data types in OhmmsPETE, such as  TinyVector&lt;T,D&gt; and  Tensor&lt;T,D&gt;
 </ul>
 Using other types are valid, as far as the operators << and >> can handle the xmlChar*.

 \endhtmlonly
 */
template<class T>
class OhmmsParameter: public OhmmsElementBase
{

  //@{
  ///reference to a value of type T
  T& ref_;
  ///the unit of this object
  std::string unit_;
  ///pointer to the correponding xmlNode
  xmlNodePtr node_;
  //@}

public:

  /*!\fn OhmmsParameter(T& a, const char* aname, const char* uname)
   *\param a the value to be referenced
   *\param aname the name of this object
   *\param uname the unit
   */
  OhmmsParameter(T& a, const char* aname, const char* uname="none"):
    OhmmsElementBase(aname), ref_(a), unit_(uname), node_(NULL)
  {
  }

  ///print to an ostream
  inline bool get(std::ostream& os) const
  {
    os << "<parameter name=\""<< myName << "\" condition=\""
       << unit_ << "\">"
       << ref_ << "</parameter>" << std::endl;
    return true;
  }

  /*!inline bool put(xmlNodePtr cur)
   *\param cur the current xmlNode whose content is assigned to ref_
   */
  inline bool put(xmlNodePtr cur)
  {
    node_ = cur;
    putContent(ref_,cur);
    return true;
  }


  ///read from istream
  inline bool put(std::istream& is)
  {
    is >> ref_;
    return true;
  }

  /*!\fn bool add(xmlNodePtr parent)
   *\param parent the parent node to which a xmlNode for this object is appended.
   *\brief This function is used by the users to add a xmlNode, when the
   *input file does not contain the corresponding <parameter/>. The content
   *of the new xmlNode is the current value of ref_.
   */
  bool add(xmlNodePtr parent)
  {
    if(!node_)
    {
      node_ =xmlNewChild(parent,parent->ns,(const xmlChar*)"parameter",NULL);
      xmlNewProp(node_,(const xmlChar*)"name",(const xmlChar*)(myName.c_str()));
      xmlNewProp(node_,(const xmlChar*)"condition",(const xmlChar*)(unit_.c_str()));
      getContent(ref_,node_);
    }
    return true;
  }

  inline void setValue(T x)
  {
    ref_=x;
  }

  ///reset member data
  inline void reset()
  {
    getContent(ref_,node_);
  }

};

/*!\class OhmmsParameter<bool>
 *\brief A specialization of OhmmsParameter<T> for T = boolean.
 */
template<>
class OhmmsParameter<bool>: public OhmmsElementBase
{

  //@{
  ///reference to a value of type T
  bool& ref_;
  ///the unit of this object
  std::string unit_;
  ///pointer to the correponding xmlNode
  xmlNodePtr node_;
  //@}

public:

  /*!\fn OhmmsParameter(bool& a, const char* aname, const char* uname)
   *\param a the boolean to be referenced.
   *\param aname the name of this object
   *\param uname the unit
   */
  OhmmsParameter(bool& a, const char* aname, const char* uname="none"):
    OhmmsElementBase(aname), ref_(a), unit_(uname), node_(NULL)
  {
  }

  ///print to an ostream
  inline bool get(std::ostream& os) const
  {
    os << "<parameter name=\""<< myName << "\">" << ref_ << "</parameter>" << std::endl;
    return true;
  }

  /*!inline bool put(xmlNodePtr cur)
   *\param cur the current xmlNode whose content is assigned to ref_
   *\brief If the content is empty, the negation of the current value is taken.
   *Example is <parameter name="force"/> to turn on the force-evaluation flag
   *of NoPropagator.
   */
  inline bool put(xmlNodePtr cur)
  {
    node_ = cur;
    const char* ac = (const char*)(xmlNodeListGetString(cur->doc, cur->xmlChildrenNode, 1));
    if(ac)
    {
      std::istringstream stream(ac);
      return stream >> ref_;
    }
    else
    {
      ref_ = !ref_;//flip the bit
      return true;
    }
  }

  inline void setValue(bool x)
  {
    ref_=x;
  }

  ///read from istream
  inline bool put(std::istream& is)
  {
    std::string yes;
    is >> yes;
    if(yes == "yes" || yes == "true" || yes == "1")
      ref_=true;
    else
      ref_ = false;
    return true;
  }

  /*!\fn bool add(xmlNodePtr parent)
   *\param parent the parent node to which a xmlNode for this object is appended.
   *\brief This function is used by the users to add a xmlNode, when the
   *input file does not contain the corresponding <parameter/>. The content
   *of the new xmlNode is the current value of ref_.
   */
  bool add(xmlNodePtr parent)
  {
    if(!node_)
    {
      node_ =xmlNewChild(parent,parent->ns,(const xmlChar*)"parameter",NULL);
      xmlNewProp(node_,(const xmlChar*)"name",(const xmlChar*)(myName.c_str()));
      xmlNewProp(node_,(const xmlChar*)"condition",(const xmlChar*)(unit_.c_str()));
      getContent(ref_,node_);
    }
    return true;
  }

  ///reset member data
  inline void reset()
  {
    getContent(ref_,node_);
  }

};
#endif /*OHMMS_OHMMSPARAMETER_H*/
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 19:14:53 -0500 (Thu, 25 Apr 2013) $
 * $Id: OhmmsParameter.h 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/
