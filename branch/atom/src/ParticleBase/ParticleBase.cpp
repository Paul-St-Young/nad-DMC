//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//
//   Jeongnim Kim
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-

namespace qmcplusplus
{
/** @file ParticleBase.cpp
 * @brief Declarations of the template functions of ParticleBase
 */
/** Default constructor */
template<class PL>
ParticleBase<PL>::ParticleBase():Counter(0), LocalNum(0), GlobalNum(0)
{
  initBase();
}

template<class PL>
void ParticleBase<PL>::initBase()
{
  R.setTypeName(ParticleTags::postype_tag);
  R.setObjName(ParticleTags::position_tag);
  ID.setTypeName(ParticleTags::indextype_tag);
  ID.setObjName(ParticleTags::id_tag);
  GroupID.setTypeName(ParticleTags::indextype_tag);
  GroupID.setObjName(ParticleTags::ionid_tag);
  //map[type-name] to enum
  AttribTypeMap[ParticleTags::indextype_tag]  = PA_IndexType;
  AttribTypeMap[ParticleTags::scalartype_tag] = PA_ScalarType;
  AttribTypeMap[ParticleTags::stringtype_tag] = PA_StringType;
  AttribTypeMap[ParticleTags::postype_tag]    = PA_PositionType;
  AttribTypeMap[ParticleTags::tensortype_tag] = PA_TensorType;
  //add basic attributes
  addAttribute(R);
  addAttribute(ID);
  addAttribute(GroupID);
  //curR is always in unit
  curR.setUnit(PosUnit::LatticeUnit);
}

template<class PL>
ParticleBase<PL>::~ParticleBase()
{
  //LOGMSG("Calling ParticleBase<PL>::~ParticleBase")
  //remove only those created by addAttribute(tname,oname)
  for(int i=0; i<myAttribList.size(); i++)
    delete myAttribList[i];
}

/** check if named attribute exists
 *  \param aname String for the attribute name
 *  \return true if the name is found.
 */
template<class PL>
bool ParticleBase<PL>::hasAttrib(const std::string& aname)
{
  return (Name2Index.find(aname) != Name2Index.end());
}

/** Search the name in VarMap and add the attribute if not found.
 *  \param tname String for the attribute type in ParticleTags
 *  \param oname String for the name of the new attribute
 *  \return the iterator of a new attribute
 */
template<class PL>
int ParticleBase<PL>::addAttribute(const std::string& tname,
                                   const std::string& oname)
{
  ///Name2Index connects the name of an object to a location in the vector of a specific type
  ///An object with oname exists. Return the value of oname
  map<string,int>::iterator it= Name2Index.find(oname);
  if(it != Name2Index.end())
    return  (*it).second;
  ///a new object to be added either INDEX, VAL, POS or TENZOR
  ///update Name2Index and AttribList
  int o_id=-1;
  int t_id=AttribTypeMap[tname];
  OhmmsObject* obj = NULL;
  if(t_id == PA_IndexType)
  {
    o_id = INDEX.size();
    INDEX.push_back(new  ParticleIndex_t(ParticleTags::indextype_tag,oname,getLocalNum()));
    obj = INDEX.back();
  }
  else
    if(t_id ==  PA_ScalarType)
    {
      o_id = VAL.size();
      VAL.push_back(new  ParticleScalar_t(ParticleTags::scalartype_tag,oname,getLocalNum()));
      obj = VAL.back();
    }
    else
      if(t_id == PA_PositionType)
      {
        o_id = POS.size();
        POS.push_back(new  ParticlePos_t(ParticleTags::postype_tag,oname,getLocalNum()));
        obj = POS.back();
      }
      else
        if(t_id == PA_TensorType)
        {
          o_id = TENZOR.size();
          TENZOR.push_back(new  ParticleTensor_t(ParticleTags::tensortype_tag,oname,getLocalNum()));
          obj = TENZOR.back();
        }
  if(obj)
  {
    Name2Index[oname]=o_id;
    obj->setID(o_id);
    myAttribList.push_back(obj);
  }
  return o_id;
}

/** Add ParticleIndex_t attribute, if not found
 * \param pa  ParticleIndex_t to be added
 * \return true the locator (iterator) of the pa in the vector<ParticleIndex_t*>
 */
template<class PL>
int
ParticleBase<PL>::addAttribute(typename ParticleBase<PL>::ParticleIndex_t& pa)
{
  map<string,int>::iterator it= Name2Index.find(pa.objName());
  if(it != Name2Index.end())
    return  (*it).second;
  if(pa.size() < getLocalNum())
    pa.resize(getLocalNum());
  int oid=Name2Index[pa.objName()]=INDEX.size();
  AttribList[pa.objName()]=&pa;
  INDEX.push_back(&pa);
  pa.setID(oid);
  return oid;
}

/** Add ParticleScalar_t attribute, if not found
 * \param pa  ParticleScalar_t to be added
 * \return true the locator (iterator) of the pa in the vector<ParticleScalar_t*>
 */
template<class PL>
int
ParticleBase<PL>::addAttribute(typename ParticleBase<PL>::ParticleScalar_t& pa)
{
  map<string,int>::iterator it= Name2Index.find(pa.objName());
  if(it != Name2Index.end())
    return  (*it).second;
  if(pa.size() < getLocalNum())
    pa.resize(getLocalNum());
  int oid=Name2Index[pa.objName()]=VAL.size();
  AttribList[pa.objName()]=&pa;
  VAL.push_back(&pa);
  pa.setID(oid);
  return oid;
}

/** Add ParticlePos_t attribute, if not found
 * \param pa  ParticlePos_t to be added
 * \return true the locator (iterator) of the pa in the vector<ParticlePos_t*>
 */
template<class PL>
int
ParticleBase<PL>::addAttribute(typename ParticleBase<PL>::ParticlePos_t& pa)
{
  map<string,int>::iterator it= Name2Index.find(pa.objName());
  if(it != Name2Index.end())
    return  (*it).second;
  if(pa.size() < getLocalNum())
    pa.resize(getLocalNum());
  int oid=Name2Index[pa.objName()]=POS.size();
  POS.push_back(&pa);
  pa.setID(oid);
  AttribList[pa.objName()]=&pa;
  return oid;
}

/** Add ParticleTensor_t attribute, if not found
 * \param pa  ParticleTensor_t to be added
 * \return true the locator (iterator) of the pa in the vector<ParticleTensor_t*>
 */
template<class PL>
int
ParticleBase<PL>::addAttribute(typename ParticleBase<PL>::ParticleTensor_t& pa)
{
  map<string,int>::iterator it= Name2Index.find(pa.objName());
  if(it != Name2Index.end())
    return  (*it).second;
  if(pa.size() < getLocalNum())
    pa.resize(getLocalNum());
  int oid=Name2Index[pa.objName()]=TENZOR.size();
  AttribList[pa.objName()]=&pa;
  TENZOR.push_back(&pa);
  pa.setID(oid);
  return oid;
}

#if defined(QMC_COMPLEX)
/** Add ParticleLaplacian_t  attribute, if not found
 * \param pa  ParticleLaplacian_t to be added
 * \return true the locator (iterator) of the pa in the vector<ParticlePos_t*>
 *
 * This function is only requred when QMC_COMPLEX is defined
 */
template<class PL>
int
ParticleBase<PL>::addAttribute(typename ParticleBase<PL>::ParticleLaplacian_t& pa)
{
  map<string,int>::iterator it= Name2Index.find(pa.objName());
  if(it != Name2Index.end())
    return  (*it).second;
  if(pa.size() < getLocalNum())
    pa.resize(getLocalNum());
  int oid=Name2Index[pa.objName()]=LAPS.size();
  LAPS.push_back(&pa);
  pa.setID(oid);
  AttribList[pa.objName()]=&pa;
  return oid;
}
/** Add ParticleGradient_t  attribute, if not found
 * \param pa  ParticleGradient_t to be added
 * \return true the locator (iterator) of the pa in the vector<ParticlePos_t*>
 *
 * This function is only requred when QMC_COMPLEX is defined
 */
template<class PL>
int
ParticleBase<PL>::addAttribute(typename ParticleBase<PL>::ParticleGradient_t& pa)
{
  map<string,int>::iterator it= Name2Index.find(pa.objName());
  if(it != Name2Index.end())
    return  (*it).second;
  if(pa.size() < getLocalNum())
    pa.resize(getLocalNum());
  int oid=Name2Index[pa.objName()]=GRADS.size();
  GRADS.push_back(&pa);
  pa.setID(oid);
  AttribList[pa.objName()]=&pa;
  return oid;
}
#endif

/*@{
 * \fn ParticlePos_t* ParticleBase<PL>::getVectorAttrib(const std::string& aname)
 * \param aname String for the attribute name
 * \return A pointer to the vector attribute with the name.
 * \brief If the named attribute does not exist, create an attribute and add to the list.
 */
template<class PL>
typename ParticleBase<PL>::ParticleIndex_t*
ParticleBase<PL>::getIndexAttrib(const std::string& aname)
{
  map<string,OhmmsObject*>::iterator it= AttribList.find(aname);
  if(it != AttribList.end())
  {
    return  dynamic_cast<ParticleIndex_t*>((*it).second);
  }
  ParticleIndex_t *pa
  = new  ParticleIndex_t(ParticleTags::indextype_tag,aname,getLocalNum());
  int oid=Name2Index[aname]=INDEX.size();
  INDEX.push_back(pa);
  AttribList[aname]=pa;
  pa->setID(oid);
  return pa;
}

template<class PL>
typename ParticleBase<PL>::ParticleScalar_t*
ParticleBase<PL>::getScalarAttrib(const std::string& aname)
{
  map<string,OhmmsObject*>::iterator it= AttribList.find(aname);
  if(it != AttribList.end())
  {
    return  dynamic_cast<ParticleScalar_t*>((*it).second);
  }
  ParticleScalar_t *pa
  = new  ParticleScalar_t(ParticleTags::scalartype_tag,aname,getLocalNum());
  int oid=Name2Index[aname]=VAL.size();
  VAL.push_back(pa);
  AttribList[aname]=pa;
  pa->setID(oid);
  return pa;
}

template<class PL>
typename ParticleBase<PL>::ParticlePos_t*
ParticleBase<PL>::getVectorAttrib(const std::string&  aname)
{
  map<string,OhmmsObject*>::iterator it= AttribList.find(aname);
  if(it != AttribList.end())
  {
    return  dynamic_cast<ParticlePos_t*>((*it).second);
  }
  ParticlePos_t *pa
  = new  ParticlePos_t(ParticleTags::postype_tag,aname,getLocalNum());
  int oid=Name2Index[aname]=POS.size();
  POS.push_back(pa);
  AttribList[aname]=pa;
  pa->setID(oid);
  return pa;
}

template<class PL>
typename ParticleBase<PL>::ParticleTensor_t*
ParticleBase<PL>::getTensorAttrib(const std::string& aname)
{
  map<string,OhmmsObject*>::iterator it= AttribList.find(aname);
  if(it != AttribList.end())
    return  dynamic_cast<ParticleTensor_t*>((*it).second);
  ParticleTensor_t *pa
  = new  ParticleTensor_t(ParticleTags::tensortype_tag,aname,getLocalNum());
  int oid=Name2Index[aname]=TENZOR.size();
  TENZOR.push_back(pa);
  AttribList[aname]=pa;
  pa->setID(oid);
  return pa;
}
/*@}*/

template<class PL>
void ParticleBase<PL>::create(unsigned m)
{
  for(int i=0; i< INDEX.size(); i++)
    INDEX[i]->create(m);
  for(int i=0; i< VAL.size(); i++)
    VAL[i]->create(m);
  for(int i=0; i< POS.size(); i++)
    POS[i]->create(m);
  for(int i=0; i< TENZOR.size(); i++)
    TENZOR[i]->create(m);
#if defined(QMC_COMPLEX)
  for(int i=0; i< GRADS.size(); i++)
    GRADS[i]->create(m);
  for(int i=0; i< LAPS.size(); i++)
    LAPS[i]->create(m);
#endif
  curR.create(m);
  LocalNum += m;
  GlobalNum += m;
}

template<class PL>
void ParticleBase<PL>::resize(unsigned m)
{
  for(int i=0; i< INDEX.size(); i++)
    INDEX[i]->resize(m);
  for(int i=0; i< VAL.size(); i++)
    VAL[i]->resize(m);
  for(int i=0; i< POS.size(); i++)
    POS[i]->resize(m);
  for(int i=0; i< TENZOR.size(); i++)
    TENZOR[i]->resize(m);
#if defined(QMC_COMPLEX)
  for(int i=0; i< GRADS.size(); i++)
    GRADS[i]->resize(m);
  for(int i=0; i< LAPS.size(); i++)
    LAPS[i]->resize(m);
#endif
  curR.resize(m);
  LocalNum = m;
  GlobalNum = m;
}

template<class PL>
void ParticleBase<PL>::clear()
{
  for(int i=0; i< INDEX.size(); i++)
    INDEX[i]->clear();
  for(int i=0; i< VAL.size(); i++)
    VAL[i]->clear();
  for(int i=0; i< POS.size(); i++)
    POS[i]->clear();
  for(int i=0; i< TENZOR.size(); i++)
    TENZOR[i]->clear();
#if defined(QMC_COMPLEX)
  for(int i=0; i< GRADS.size(); i++)
    GRADS[i]->clear();
  for(int i=0; i< LAPS.size(); i++)
    LAPS[i]->clear();
#endif
  curR.clear();
  LocalNum = 0;
  GlobalNum = 0;
}

/**function to create N-particle system
 *@param agroup an integer array containing the number of particles belonging
 *to a subgroup.
 *@brief allocate the ParticleAttributes, such as R, G, L.
 *The size of n is the number of distinct subgroups. SubPtcl class
 *is used to efficient evaluate the subgroup properties.
 */
template<class PL>
void ParticleBase<PL>::create(const vector<int>& agroup)
{
  SubPtcl.resize(agroup.size()+1);
  SubPtcl[0] = 0;
  for(int is=0; is<agroup.size(); is++)
    SubPtcl[is+1] = SubPtcl[is]+agroup[is];
  int nsum = SubPtcl[agroup.size()];
  resize(nsum);
  int loc=0;
  for(int i=0; i<agroup.size(); i++)
  {
    for(int j=0; j<agroup[i]; j++,loc++)
      GroupID[loc] = i;
  }
}

}

/*
template<class PL>
void ParticleBase<PL>::update(int imode) {

  if(imode == 1) { // a new ordered configuration

    Lattice.update(this,imode);
    ////////////////////////////////////////////////
    //\todo MPI distribution, global index should be assigned
    ////////////////////////////////////////////////
    for(int i=0; i<getLocalNum(); i++) ID[i] = i;

  } else if(imode == 2) {

    ////////////////////////////////////////////////////////
    // initial cell assignmenet
    ////////////////////////////////////////////////////////
    if(Lattice.dGrid[2]->getTotalNum() == 1) { //1x1x1 grid, do nothing

      // simple partition of particles accoring to the order
      Lattice.dGrid[0]->distribute(getLocalNum());
      Lattice.dGrid[1]->distribute(getLocalNum());
      Lattice.dGrid[2]->distribute(getLocalNum());
      return;
    }

    if(Lattice.dGrid[1]->PtclDist.empty()) {
      // simple partition of particles accoring to the order
      Lattice.dGrid[0]->distribute(getLocalNum());
      Lattice.dGrid[1]->distribute(getLocalNum());
    }

    vector<int> nat(Lattice.dGrid[2]->getTotalNum(),0);
    ParticleIndex_t cellid(getLocalNum());
    ParticleIndex_t idtmp(getLocalNum());
    ParticlePos_t   ptmp(getLocalNum());
    idtmp = GroupID;
    ptmp = R;

    for(int iG=0; iG<getLocalNum(); iG++){
      //int iloc = Lattice.dGrid[2]->loc(Lattice.toUnit(R[iG]));
      int iloc = Lattice.dGrid[2]->loc(R[iG]);
      cellid[iG] = iloc;
      nat[iloc]++;
    }

    Lattice.dGrid[2]->PtclDist.resize(nat.size()+1);
    Lattice.dGrid[2]->PtclDist[0] = 0;
    for(int ic=0; ic<nat.size(); ic++) {
      Lattice.dGrid[2]->PtclDist[ic+1] =
	Lattice.dGrid[2]->PtclDist[ic] + nat[ic];
    }

    int ntot = 0;
    for(int ic=0; ic<Lattice.dGrid[2]->getTotalNum(); ic++) {
      for(int iG=0; iG<getLocalNum(); iG++) {
	if(cellid[iG] == ic) {
	  R[ntot] = ptmp[iG];
	  GroupID[ntot] = idtmp[iG];
	  ID[iG] = ntot; // record the backward map
	  ntot++;
	}
      }
    }
  }
}

template<class PL>
void ParticleBase<PL>::assign(const This_t& ptclin){
  create(ptclin.getLocalNum());
  Lattice = ptclin.Lattice;
  R = ptclin.R;
  ID = ptclin.ID;
  GroupID = ptclin.GroupID;
}

template<class PL>
//void ParticleBase<PL>::update(const ParticleBase<PL>::UpdateMode_t& update) {
void ParticleBase<PL>::update(const UpdateMode_t& ptclupdate) {
  //update[0] -> MPI redistribute
  //update[1] -> OMP redistribute
  //update[2] -> Local reorder
  //update[3] -> reassign ID
  //update[4] -> other field
  if(ptclupdate[2]) {

    ////////////////////////////////////////////////////////
    // initial cell assignmenet
    ////////////////////////////////////////////////////////
    vector<int> nat(Lattice.dGrid[2]->getTotalNum(),0);
#pragma omp parallel
    {
      int myID = omp_get_thread_num();
      typename PL::PtclGrid_t& basegrid = *(Lattice.dGrid[2]);
      typename PL::PtclGrid_t& ompgrid = *(Lattice.dGrid[1]);

      SingleParticlePos_t pos;
      int ni = ompgrid.PtclDist[myID];
      int nf = ompgrid.PtclDist[myID+1];
      vector<int> cellid(nf-ni);
      int iL = 0;
      for(int iG=ni; iG<nf; iG++,iL++){
	pos = Lattice.toUnit(R[iG]);
	int cloc = basegrid.loc(pos);
	cellid[iL] = cloc;
	nat[cloc]++;
      }

      vector<int> idtmp(nf-ni);
      vector<SingleParticlePos_t> ptmp(nf-ni);
      int il=0;
      for(int iG=ni; iG<nf; il++,iG++) idtmp[il] = GroupID[iG];
      il=0;
      for(int iG=ni; iG<nf; il++,iG++) ptmp[il] = R[iG];

      int ntot = ni;
      for(int ic=basegrid.NodeDist[myID]; ic<basegrid.NodeDist[myID+1]; ic++) {
	int iL = 0;
	for(int iG=ni; iG<nf; iL++, iG++) {
	  if(cellid[iL] == ic) {
            R[ntot] = ptmp[iL];
            GroupID[ntot] = idtmp[iL];
	    ID[ntot] = iG; // record the backward map
            ntot++;
	  }
	}
      }
      /////////////////
      //remapping of other particle attributes
      /////////////////
    }//end-of-omp

    Lattice.dGrid[2]->PtclDist.resize(nat.size()+1);
    Lattice.dGrid[2]->PtclDist[0] = 0;
    for(int ic=0; ic<nat.size(); ic++) {
      Lattice.dGrid[2]->PtclDist[ic+1] =
	Lattice.dGrid[2]->PtclDist[ic] + nat[ic];
    }
  }
  //#endif
  if(ptclupdate[3])
    for(int i=0; i<getLocalNum(); i++) ID[i] = i;
}
*/

// template<class PL>
// bool ParticleBase<PL>::write_data(ostream&) {
//   return true;
// }
// template<class PL>
// bool ParticleBase<PL>::read_data(istream&) {
//   return true;
// }

/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 19:14:53 -0500 (Thu, 25 Apr 2013) $
 * $Id: ParticleBase.cpp 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/

