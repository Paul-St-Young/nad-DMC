//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
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
/**@file QMCDriverFactory.cpp
 * @brief Implments QMCMain operators.
 */
#include "QMCApp/QMCDriverFactory.h"
#include "QMCApp/WaveFunctionPool.h"
#include "QMCApp/HamiltonianPool.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCDrivers/VMC/VMCFactory.h"
#include "QMCDrivers/DMC/DMCFactory.h"
#include "QMCDrivers/RMC/RMCFactory.h"
#include "QMCDrivers/QMCOptimize.h"
#include "QMCDrivers/QMCFixedSampleLinearOptimize.h"
#include "QMCDrivers/QMCCorrelatedSamplingLinearOptimize.h"
#include "QMCDrivers/WaveFunctionTester.h"
#include "Utilities/OhmmsInfo.h"
#include <queue>
using namespace std;
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"

namespace qmcplusplus
{

///initialize the static data member
//ParticleSetPool* QMCDriverFactory::ptclPool = new ParticleSetPool;
QMCDriverFactory::QMCDriverFactory(Communicate* c): MPIObjectBase(c),
  qmcSystem(0), qmcDriver(0) , curRunType(DUMMY_RUN)
{
  ////create ParticleSetPool
  ptclPool = new ParticleSetPool(myComm);
  //create WaveFunctionPool
  psiPool = new WaveFunctionPool(myComm);
  psiPool->setParticleSetPool(ptclPool);
  //create HamiltonianPool
  hamPool = new HamiltonianPool(myComm);
  hamPool->setParticleSetPool(ptclPool);
  hamPool->setWaveFunctionPool(psiPool);
}

QMCDriverFactory::~QMCDriverFactory()
{
  delete hamPool;
  delete psiPool;
  delete ptclPool;
}

void QMCDriverFactory::putCommunicator(xmlNodePtr cur)
{
  //BROKEN: myComm is ALWAYS initialized by the constructor
  if(myComm)
    return;
  ParameterSet params;
  int nparts=1;
  params.add(nparts,"groups","int");
  params.add(nparts,"twistAngles","int");
  params.put(cur);
  if(nparts>1)
  {
    app_log() << "  Communiator groups = " << nparts << endl;
    myComm = new Communicate(*OHMMS::Controller,nparts);
  }
}

bool QMCDriverFactory::setQMCDriver(int curSeries, xmlNodePtr cur)
{
  string curName((const char*)cur->name);
  string update_mode("pbyp");
  string qmc_mode("invalid");
  string multi_tag("no");
  string warp_tag("no");
  string append_tag("no");
#if defined(QMC_CUDA)
  string gpu_tag("yes");
#else
  string gpu_tag("no");
#endif
  OhmmsAttributeSet aAttrib;
  aAttrib.add(qmc_mode,"method");
  aAttrib.add(update_mode,"move");
  aAttrib.add(multi_tag,"multiple");
  aAttrib.add(warp_tag,"warp");
  aAttrib.add(append_tag,"append");
  aAttrib.add(gpu_tag,"gpu");
  aAttrib.put(cur);
  bitset<QMC_MODE_MAX>  WhatToDo;
  bool append_run =(append_tag == "yes");
  WhatToDo[SPACEWARP_MODE]= (warp_tag == "yes");
  WhatToDo[MULTIPLE_MODE]= (multi_tag == "yes");
  WhatToDo[UPDATE_MODE]= (update_mode == "pbyp");
#if defined(QMC_CUDA)
  WhatToDo[GPU_MODE      ] = (gpu_tag     == "yes");
#endif
  OhmmsInfo::flush();
  QMCRunType newRunType = DUMMY_RUN;
  if(curName != "qmc")
    qmc_mode=curName;
  int nchars=qmc_mode.size();
  if(qmc_mode.find("linear") < nchars)
  {
    if (qmc_mode.find("cslinear") < nchars)
      newRunType=CS_LINEAR_OPTIMIZE_RUN;
    else
      newRunType=LINEAR_OPTIMIZE_RUN;
  }
  else if(qmc_mode.find("opt") < nchars)
  {
    newRunType=OPTIMIZE_RUN;
  }
  else
  {
    if(qmc_mode.find("ptcl")<nchars)
      WhatToDo[UPDATE_MODE]=1;
    if(qmc_mode.find("mul")<nchars)
      WhatToDo[MULTIPLE_MODE]=1;
    if(qmc_mode.find("warp")<nchars)
      WhatToDo[SPACEWARP_MODE]=1;
//       if (qmc_mode.find("rmcPbyP")<nchars)
//       {
//         newRunType=RMC_PBYP_RUN;
//       }
//       else
    if(qmc_mode.find("rmc")<nchars)
    {
      newRunType=RMC_RUN;
    }
    else if(qmc_mode.find("vmc")<nchars)
    {
      newRunType=VMC_RUN;
    }
    else if(qmc_mode.find("dmc")<nchars)
    {
      newRunType=DMC_RUN;
    }
  }
  unsigned long newQmcMode=WhatToDo.to_ulong();
  //initialize to 0
  QMCDriver::BranchEngineType* branchEngine=0;
  if(qmcDriver)
  {
    if( newRunType != curRunType || newQmcMode != curQmcMode)
    {
      if(curRunType == DUMMY_RUN)
      {
        APP_ABORT("QMCDriverFactory::setQMCDriver\n Other qmc sections cannot come after <qmc method=\"test\">.\n");
      }
      //pass to the new driver
      branchEngine=qmcDriver->getBranchEngine();
      delete qmcDriver;
      //set to 0 so that a new driver is created
      qmcDriver = 0;
      //if the current qmc method is different from the previous one, append_run is set to false
      append_run = false;
    }
    else
    {
      app_log() << "  Reusing " << qmcDriver->getEngineName() << endl;
      //         if(curRunType == DMC_RUN)
      qmcDriver->resetComponents(cur);
    }
  }
  if(curSeries == 0)
    append_run = false;
  //continue with the existing qmcDriver
  if(qmcDriver)
    return append_run;
  //need to create a qmcDriver
  curRunType = newRunType;
  curQmcMode = newQmcMode;
  curQmcModeBits = WhatToDo;
  //create a driver
  createQMCDriver(cur);
  //initialize QMCDriver::myComm
  qmcDriver->initCommunicator(myComm);
  //branchEngine has to be transferred to a new QMCDriver
  if(branchEngine)
    qmcDriver->setBranchEngine(branchEngine);
  OhmmsInfo::flush();
  return append_run;
}

void QMCDriverFactory::createQMCDriver(xmlNodePtr cur)
{
  ///////////////////////////////////////////////
  // get primaryPsi and primaryH
  ///////////////////////////////////////////////
  TrialWaveFunction* primaryPsi= 0;
  QMCHamiltonian* primaryH=0;
  queue<TrialWaveFunction*> targetPsi;//FIFO
  queue<QMCHamiltonian*> targetH;//FIFO
  xmlNodePtr tcur=cur->children;
  while(tcur != NULL)
  {
    if(xmlStrEqual(tcur->name,(const xmlChar*)"qmcsystem"))
    {
      const xmlChar* t= xmlGetProp(tcur,(const xmlChar*)"wavefunction");
      if(t != NULL)
      {
        targetPsi.push(psiPool->getWaveFunction((const char*)t));
      }
      else
      {
        app_warning() << " qmcsystem does not have wavefunction. Assign 0" << endl;
        targetPsi.push(0);
      }
      t= xmlGetProp(tcur,(const xmlChar*)"hamiltonian");
      if(t != NULL)
      {
        targetH.push(hamPool->getHamiltonian((const char*)t));
      }
      else
      {
        app_warning() << " qmcsystem does not have hamiltonian. Assign 0" << endl;
        targetH.push(0);
      }
    }
    tcur=tcur->next;
  }
  //mark the first targetPsi and targetH as the primaries
  if(targetH.empty())
  {
    primaryPsi=psiPool->getPrimary();
    primaryH=hamPool->getPrimary();
  }
  else
  {
    primaryPsi=targetPsi.front();
    targetPsi.pop();
    primaryH=targetH.front();
    targetH.pop();
  }
  //set primaryH->Primary
  primaryH->setPrimary(true);
  ////flux is evaluated only with single-configuration VMC
  //if(curRunType == VMC_RUN && !curQmcModeBits[MULTIPLE_MODE])
  //{
  //  QMCHamiltonianBase* flux=primaryH->getHamiltonian("Flux");
  //  if(flux == 0) primaryH->addOperator(new ConservedEnergy,"Flux");
  //}
  //else
  //{
  //  primaryH->remove("Flux");
  //}
  //(SPACEWARP_MODE,MULTIPE_MODE,UPDATE_MODE)
  if(curRunType == VMC_RUN)
  {
    //VMCFactory fac(curQmcModeBits[UPDATE_MODE],cur);
    VMCFactory fac(curQmcModeBits.to_ulong(),cur);
    qmcDriver = fac.create(*qmcSystem,*primaryPsi,*primaryH,*ptclPool,*hamPool,*psiPool);
    //TESTING CLONE
    //TrialWaveFunction* psiclone=primaryPsi->makeClone(*qmcSystem);
    //qmcDriver = fac.create(*qmcSystem,*psiclone,*primaryH,*ptclPool,*hamPool);
  }
  else if(curRunType == DMC_RUN)
  {
    DMCFactory fac(curQmcModeBits[UPDATE_MODE],
                   curQmcModeBits[GPU_MODE], cur);
    qmcDriver = fac.create(*qmcSystem,*primaryPsi,*primaryH,*hamPool,*psiPool);
  }
#if QMC_BUILD_LEVEL>1
  else if(curRunType == RMC_RUN)
  {
    RMCFactory fac(curQmcModeBits[UPDATE_MODE], cur);
    qmcDriver = fac.create(*qmcSystem,*primaryPsi,*primaryH,*ptclPool,*hamPool,*psiPool);
  }
#endif
//#if QMC_BUILD_LEVEL>1
//    else if(curRunType == RMC_RUN)
//    {
//      app_log() << "Using RQMCMultiple: no warping, no pbyp" << endl;
//      qmcDriver = new RQMCMultiple(*qmcSystem,*primaryPsi,*primaryH,*psiPool);
//    }
//#endif
//#if defined(QMC_BUILD_COMPLETE)
//    else if(curRunType == RMC_RUN)
//    {
//#if defined(QMC_COMPLEX)
//      qmcDriver = new RQMCMultiple(*qmcSystem,*primaryPsi,*primaryH);
//#else
//      if(curQmcModeBits[SPACEWARP_MODE])
//        qmcDriver = new RQMCMultiWarp(*qmcSystem,*primaryPsi,*primaryH, *ptclPool);
//      else
//        qmcDriver = new RQMCMultiple(*qmcSystem,*primaryPsi,*primaryH);
//    }
//    else if (curRunType==RMC_PBYP_RUN)
//    {
//      qmcDriver = new RQMCMultiplePbyP(*qmcSystem,*primaryPsi,*primaryH);
//#endif
//    }
//#endif
  else if(curRunType == OPTIMIZE_RUN)
  {
    QMCOptimize *opt = new QMCOptimize(*qmcSystem,*primaryPsi,*primaryH,*hamPool,*psiPool);
    //ZeroVarianceOptimize *opt = new ZeroVarianceOptimize(*qmcSystem,*primaryPsi,*primaryH );
    opt->setWaveFunctionNode(psiPool->getWaveFunctionNode("psi0"));
    qmcDriver=opt;
  }
  else if(curRunType == LINEAR_OPTIMIZE_RUN)
  {
    QMCFixedSampleLinearOptimize *opt = new QMCFixedSampleLinearOptimize(*qmcSystem,*primaryPsi,*primaryH,*hamPool,*psiPool);
    //ZeroVarianceOptimize *opt = new ZeroVarianceOptimize(*qmcSystem,*primaryPsi,*primaryH );
    opt->setWaveFunctionNode(psiPool->getWaveFunctionNode("psi0"));
    qmcDriver=opt;
  }
  else if(curRunType == CS_LINEAR_OPTIMIZE_RUN)
  {
#if defined(QMC_CUDA)
    app_log() << "cslinear is not supported. Switch to linear method. " << endl;
    QMCFixedSampleLinearOptimize *opt = new QMCFixedSampleLinearOptimize(*qmcSystem,*primaryPsi,*primaryH,*hamPool,*psiPool);
#else
    QMCCorrelatedSamplingLinearOptimize *opt = new QMCCorrelatedSamplingLinearOptimize(*qmcSystem,*primaryPsi,*primaryH,*hamPool,*psiPool);
#endif
    opt->setWaveFunctionNode(psiPool->getWaveFunctionNode("psi0"));
    qmcDriver=opt;
  }
  else
  {
    WARNMSG("Testing wavefunctions. Creating WaveFunctionTester for testing");
    qmcDriver = new WaveFunctionTester(*qmcSystem,*primaryPsi,*primaryH,
                                       *ptclPool,*psiPool);
  }
  if(curQmcModeBits[MULTIPLE_MODE])
  {
    while(targetH.size())
    {
      qmcDriver->add_H_and_Psi(targetH.front(),targetPsi.front());
      targetH.pop();
      targetPsi.pop();
    }
  }
}
}
/***********************************************************************
 * $RCSfilMCMain.cpp,v $   $Author: jnkim $
 * $Revision: 5825 $   $Date: 2013-05-06 20:01:12 -0500 (Mon, 06 May 2013) $
 * $Id: QMCDriverFactory.cpp 5825 2013-05-07 01:01:12Z jnkim $
 ***************************************************************************/
