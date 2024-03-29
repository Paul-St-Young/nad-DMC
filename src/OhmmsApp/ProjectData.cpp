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

#include "OhmmsApp/ProjectData.h"
#include "Message/Communicate.h"
#include "Platforms/sysutil.h"

namespace qmcplusplus
{

//----------------------------------------------------------------------------
// ProjectData
//----------------------------------------------------------------------------
// constructors and destructors
ProjectData::ProjectData(const char* aname):
  m_title("asample"),
  m_user("none"),
  m_host("none"),
  m_date("none"),
  m_series(0),
  m_cur(0)
{
  myComm=OHMMS::Controller;
  if(aname==0)
  {
    m_title=getDateAndTime("%Y%m%dT%H%M");
    setName(m_title);
  }
  else
    setName(aname);
}

void  ProjectData::setCommunicator(Communicate* c)
{
  myComm=c;
}

bool ProjectData::get(ostream& os) const
{
  os << "  Project = " << m_title << "\n";
  os << "  date    = " << getDateAndTime("%Y-%m-%d %H:%M:%S %Z\n");
  os << "  host    = " << m_host << "\n";
  os << "  user    = " << m_user << "\n";
  return true;
}

bool ProjectData::put(istream& is)
{
#if defined(ENABLE_GUI)
  // get the data from window
  wxTextCtrl* temp = (wxTextCtrl*)(wxWindow::FindWindowById(ID_PROJECT_TITLE));
  m_title = temp->GetValue();
  temp = (wxTextCtrl*)(wxWindow::FindWindowById(ID_PROJECT_SID));
  wxString t = temp->GetValue();
  m_series = atoi(t.c_str());
  temp = (wxTextCtrl*)(wxWindow::FindWindowById(ID_DATE));
  wxDateTime now = wxDateTime::Now();
  m_date = now.FormatISODate();
  temp->SetValue(m_date.c_str());
  temp = (wxTextCtrl*)(wxWindow::FindWindowById(ID_HOST));
  m_host = wxGetFullHostName();
  temp->SetValue(m_host.c_str());
  temp = (wxTextCtrl*)(wxWindow::FindWindowById(ID_USER_ID));
  m_user = wxGetUserId();
  temp->SetValue(m_user.c_str());
#else
  string t1;
  while(!is.eof())
  {
    if(isdigit(t1[0]))
      m_series = atoi(t1.c_str());
    is >> t1;
    if(t1 == "series")
      is >> m_series;
    else
      if(t1 == "user")
        is >> m_user;
      else
        if(t1 == "host")
          is >> m_host;
        else
          if(t1 == "date")
            is >> m_date;
          else
            m_title = t1;
  }
#endif
  reset();
  return true;
}

void ProjectData::advance()
{
  m_series++;
  reset();
}

void ProjectData::rewind()
{
  if(m_series>0)
    m_series--;
  reset();
}

/**\fn void ProjectData::reset()
 *\brief Construct the root name with m_title and m_series.
 */
void ProjectData::reset()
{
  int nproc_g = OHMMS::Controller->size();
  int nproc = myComm->size();
  int nodeid = myComm->rank();
  int groupid=myComm->getGroupID();
  char fileroot[256], nextroot[256];
  if(nproc_g == nproc)
    sprintf(fileroot,"%s.s%03d",m_title.c_str(),m_series);
  else
    sprintf(fileroot,"%s.g%03d.s%03d",m_title.c_str(),groupid,m_series);
  m_projectmain=fileroot;
  //set the communicator name
  myComm->setName(fileroot);
  if(nproc_g == nproc)
  {
    if(nproc > 1)
    {
      sprintf(fileroot,".s%03d.p%03d", m_series,nodeid);
      sprintf(nextroot,".s%03d.p%03d", m_series+1,nodeid);
    }
    else
    {
      sprintf(fileroot,".s%03d", m_series);
      sprintf(nextroot,".s%03d", m_series+1);
    }
  }
  else
  {
    if(nproc > 1)
    {
      sprintf(fileroot,".g%03d.s%03d.p%03d", groupid,m_series,nodeid);
      sprintf(nextroot,".g%03d.s%03d.p%03d", groupid,m_series+1,nodeid);
    }
    else
    {
      sprintf(fileroot,".g%03d.s%03d", groupid,m_series);
      sprintf(nextroot,".g%03d.s%03d", groupid,m_series+1);
    }
  }
  m_projectroot = m_title;
  m_projectroot.append(fileroot);
  m_nextroot = m_title;
  m_nextroot.append(nextroot);
  std::stringstream s;
  s << m_series+1;
  if(m_cur)
    xmlSetProp(m_cur, (const xmlChar *) "series", (const xmlChar *)(s.str().c_str()));
}

bool ProjectData::PreviousRoot(string& oldroot) const
{
  oldroot.erase(oldroot.begin(), oldroot.end());
  if(m_series)
  {
    char fileroot[128];
    int nproc_g = OHMMS::Controller->size();
    int nproc = myComm->size();
    int nodeid = myComm->rank();
    int groupid=myComm->getGroupID();
    if(nproc_g == nproc)
    {
      if(nproc > 1)
        sprintf(fileroot,".s%03d.p%03d", m_series-1,nodeid);
      else
        sprintf(fileroot,".s%03d", m_series-1);
    }
    else
    {
      if(nproc > 1)
        sprintf(fileroot,".g%03d.s%03d.p%03d", groupid,m_series-1,nodeid);
      else
        sprintf(fileroot,".g%03d.s%03d",groupid,m_series-1);
    }
    oldroot = m_title;
    oldroot.append(fileroot);
    return true;
  }
  else
  {
    return false;
  }
}

#if defined(HAVE_LIBXML2)

bool ProjectData::put(xmlNodePtr cur)
{
  m_cur = cur;
  xmlDocPtr doc = cur->doc;
  m_title = (const char*)(xmlGetProp(cur, (const xmlChar *) "id"));
  m_series = atoi((const char*)(xmlGetProp(cur, (const xmlChar *) "series")));
  ///first, overwrite the existing xml nodes
  cur = cur->xmlChildrenNode;
  while (cur != NULL)
  {
    string cname((const char*)(cur->name));
    if(cname == "user")
    {
      m_user = getUserName();
      xmlNodeSetContent(cur,(const xmlChar*)(m_user.c_str()));
    }
    if (cname == "host")
    {
      m_host = getHostName();
      xmlNodeSetContent(cur,(const xmlChar*)(m_host.c_str()));
    }
    if (cname == "date")
    {
      m_date = getDateAndTime();
      xmlNodeSetContent(cur,(const xmlChar*)(m_date.c_str()));
    }
    cur = cur->next;
  }
  ///second, add xml nodes, if missing
  if(m_host == "none")
  {
    m_host =  getHostName();
    xmlNewChild(m_cur,m_cur->ns,(const xmlChar*)"host",
                (const xmlChar*)(m_host.c_str()));
  }
  if(m_date == "none")
  {
    m_date = getDateAndTime();
    xmlNewChild(m_cur,m_cur->ns,
                (const xmlChar*)"date",(const xmlChar*)(m_date.c_str()));
  }
  if(m_user == "none")
  {
    m_user = getUserName();
    xmlNewChild(m_cur,m_cur->ns,
                (const xmlChar*)"user",(const xmlChar*)(m_user.c_str()));
  }
  reset();
  return true;
}

#endif
}
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 19:14:53 -0500 (Thu, 25 Apr 2013) $
 * $Id: ProjectData.cpp 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/
