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
#include "QMCWaveFunctions/Fermion/MultiSlaterDeterminantFast.h"
#include "QMCWaveFunctions/Fermion/MultiDiracDeterminantBase.h"
#include "ParticleBase/ParticleAttribOps.h"

namespace qmcplusplus
{

MultiSlaterDeterminantFast::MultiSlaterDeterminantFast(ParticleSet& targetPtcl, MultiDiracDeterminantBase* up, MultiDiracDeterminantBase* dn):
  RatioTimer("MultiSlaterDeterminantFast::ratio"),
  RatioGradTimer("MultiSlaterDeterminantFast::ratioGrad"),
  RatioAllTimer("MultiSlaterDeterminantFast::ratio(all)"),
  Ratio1Timer("MultiSlaterDeterminantFast::detEval_ratio"),
  Ratio1GradTimer("MultiSlaterDeterminantFast::detEval_ratioGrad"),
  Ratio1AllTimer("MultiSlaterDeterminantFast::detEval_ratio(all)"),
  UpdateTimer("MultiSlaterDeterminantFast::updateBuffer"),
  EvaluateTimer("MultiSlaterDeterminantFast::evaluate"),
  AccRejTimer("MultiSlaterDeterminantFast::Accept_Reject")
{
  registerTimers();
  //Optimizable=true;
  Optimizable=true;
  OrbitalName="MultiSlaterDeterminantFast";
  usingCSF=false;
  NP = targetPtcl.getTotalNum();
  nels_up = targetPtcl.last(0)-targetPtcl.first(0);
  nels_dn = targetPtcl.last(1)-targetPtcl.first(1);
  FirstIndex_up=targetPtcl.first(0);
  FirstIndex_dn=targetPtcl.first(1);
  Dets.resize(2);
  Dets[0]=up;
  Dets[1]=dn;
  myG.resize(NP);
  myL.resize(NP);
  myG_temp.resize(NP);
  myL_temp.resize(NP);
  DetID.resize(NP);
  for(int i=0; i<targetPtcl.groups(); ++i)
    for(int j=targetPtcl.first(i); j<targetPtcl.last(i); ++j)
      DetID[j]=i;
  usingBF=false;
  BFTrans=0;
}

OrbitalBasePtr MultiSlaterDeterminantFast::makeClone(ParticleSet& tqp) const
{
  MultiDiracDeterminantBase* up_clone = new MultiDiracDeterminantBase(*Dets[0]);
  MultiDiracDeterminantBase* dn_clone = new MultiDiracDeterminantBase(*Dets[1]);
  MultiSlaterDeterminantFast* clone = new MultiSlaterDeterminantFast(tqp,up_clone,dn_clone);
  if(usingBF)
  {
    BackflowTransformation *tr = BFTrans->makeClone(tqp);
    clone->setBF(tr);
  }
  clone->resetTargetParticleSet(tqp);
  clone->C2node_up=C2node_up;
  clone->C2node_dn=C2node_dn;
  clone->Optimizable=Optimizable;
  clone->C=C;
  clone->myVars=myVars;
  clone->usingCSF=usingCSF;
  clone->usingBF=usingBF;
  if (usingCSF)
  {
    clone->CSFcoeff=CSFcoeff;
    clone->CSFexpansion=CSFexpansion;
    clone->DetsPerCSF=DetsPerCSF;
  }
  return clone;
}

MultiSlaterDeterminantFast::~MultiSlaterDeterminantFast() { }

void MultiSlaterDeterminantFast::resetTargetParticleSet(ParticleSet& P)
{
  if(usingBF)
  {
    BFTrans->resetTargetParticleSet(P);
    for(int i=0; i<Dets.size(); i++)
      Dets[i]->resetTargetParticleSet(BFTrans->QP);
  }
  else
  {
    for(int i=0; i<Dets.size(); i++)
      Dets[i]->resetTargetParticleSet(P);
  }
}

//  void MultiSlaterDeterminantFast::resize(int n1, int n2)
//  {
//  }

void MultiSlaterDeterminantFast::testMSD(ParticleSet& P, int iat)
{
//     APP_ABORT("Testing disabled for safety");
  app_log() <<"Testing MSDFast. \n";
  int n = nels_up+nels_dn;
  ParticleSet::ParticleGradient_t G(n),G0(n);
  ParticleSet::ParticleLaplacian_t L(n),L0(n);
  ValueType log, log0;
//     log = msd->evaluate(P,G,L);
  log0 = evaluate(P,G0,L0);
  /*
       app_log() <<"Testing evaluate(P,G,L). \n";
       cout<<endl <<endl;
       cout<<"Psi: " <<log <<"   " <<log0 <<"   " <<log/log0 <<endl;

       for(int i=0; i<n; i++) {
         cout<<i  <<"\n"
             <<"  x: " <<G(i)[0]-G0(i)[0] <<"\n"
             <<"  y: " <<G(i)[1]-G0(i)[1] <<"\n"
             <<"  z: " <<G(i)[2]-G0(i)[2] <<"\n"
             <<"  d2: " <<L(i)-L0(i) <<"\n"
             <<endl;
       }
       cout<<endl <<endl;
       APP_ABORT("end of test 1");
  */
  Walker_t::Buffer_t wbuffer;
  wbuffer.clear();
  log=registerData(P,wbuffer);
//     log = msd->evaluate(P,G,L);
  log0 = evaluate(P,G0,L0);
  PosType dr;
  dr[0] = 0.1;
  dr[1]=0.05;
  dr[2] = -0.01;
  PosType newpos(P.makeMove(iat,dr));
  app_log() <<"Testing ratio(P,dG,dL). \n";
  G=0;
  G0=0;
  L=0;
  L0=0;
//     log = msd->ratio(P,iat,G,L);
  log0 = ratio(P,iat,G0,L0);
  cout<<"Psi: " <<log <<"   " <<log0 <<"   " <<log/log0 <<endl;
  for(int i=0; i<n; i++)
  {
    cout<<i  <<"\n"
        <<"  x: " <<G(i)[0]-G0(i)[0] <<"  " <<G(i)[0]   <<"\n"
        <<"  y: " <<G(i)[1]-G0(i)[1] <<"  " <<G(i)[1] <<"\n"
        <<"  z: " <<G(i)[2]-G0(i)[2] <<"  " <<G(i)[2] <<"\n"
        <<"  d2: " <<L(i)-L0(i) <<"  " <<L(i) <<"\n"
        <<endl;
  }
  cout<<endl <<endl;
  APP_ABORT("After MultiSlaterDeterminantFast::testMSD()");
}

void MultiSlaterDeterminantFast::updateCoeff(RealType R){
  // update determinant coeffients with ion position

  // !!! hard-code coefficient interpolation
  RealType mych20[] = 
{
  -0.956538,-0.0997944,0.070268,-0.0648965,-0.0407198,-0.0243953,0.0243953,-0.0487907,0.00840395,0.00840395,0.0168079,-0.0234341,-0.0234341,-0.0213796,-0.0213796,-0.0202011,0.00830108,0.00830108,-0.00830108,-0.00830108,0.0166022,0.0166022,0.02033,0.02033,0.00469331,-0.00469331,0.00938663,-0.0183559,-0.0183559,-0.0133941,-0.0133941,-0.0139433,-0.0139433,0.0148545,0.0148545,0.0169626,0.0169626,-0.0143823,-0.0118947,-0.0118947,0.00340511,-0.00170256,0.00170256,-0.00170256,0.00170256,-0.00340511,0.00340511,-0.00340511,0.00681023,0.0105784,0.0105784,0.0211569,-0.0074983,0.0074983,0.0149966,-0.00680168,0.00680168,-0.0136034,0.0205972,0.0205972,0.00582225,-0.00291113,0.00291113,-0.00291113,0.00291113,-0.00582225,0.00582225,-0.00582225,0.0116445,-0.0083734,-0.0083734,-0.0169486,0.00809171,0.00809171,0.00809171,0.00809171,-0.0112299,-0.0112299,-0.00701585,-0.00701585,-0.0140317,-0.0119935,-0.010752,-0.010752,-0.00316951,0.00316951,0.00633903,-0.00686289,-0.00686289,0.010894,0.010894,-0.00697884,-0.00697884,-0.00534885,0.00534885,0.00844805,0.00844805,0.00545893,-0.00545893,-0.0109179,-0.0138126,
  0.00387314,0.00387314,-0.00271256,-0.00271256,0.00271256,0.00271256,-0.00542512,-0.00542512,-0.00972768,-0.0113759,-0.0104834,-0.00204338,0.00204338,0.00408676,0.0136913,0.000447654,-0.000447654,-0.000895309,0.00423399,-0.00423399,0.00846798,0.00416214,-0.00208107,0.00208107,0.00208107,-0.00208107,0.00416214,0.0152539,0.00385381,0.00385381,0.00770763,-0.00971754,0.00548338,0.00548338,0.00972216,-0.00245447,0.00245447,-0.00490894,0.00361912,-0.00361912,0.00723826,0.00657623,0.00657623,0.00397124,-0.00397124,0.0079425,0.00185768,-0.00092884,0.00092884,-0.00092884,0.00092884,-0.00185768,0.00185768,-0.00185768,0.00371536,-0.00615077,0.00822004,0.00822004,-0.00513809,0.00513809,0.0102762,0.00395013,0.00395013,0.00790027,-0.0015996,0.000799802,-0.000799802,0.000799802,-0.000799802,0.0015996,-0.0015996,0.0015996,-0.0031992,-0.00697282,-0.00548458,-0.00793559,0.00793559,0.00405394,0.00405394,0.00478325,-0.00239163,0.00239163,0.00239163,-0.00239163,0.00478325,-0.00654511,-0.00654511,7.65655e-05,-7.65655e-05,-0.000153131,0.0138027,0.000465793,0.000465793,-0.000465793,-0.000465793,0.000931586,0.000931586,0.00779886,0.00779886,0.00480364,
  0.00480364,-0.00163229,0.00163229,0.00326458,-0.00789706,-0.00461175,-0.00461175,-0.00922351,0.00400209,-0.00400209,-0.00800419,0.00301866,-0.00301866,0.00603733,-0.00101484,-0.00101484,0.00101484,0.00101484,-0.00202968,-0.00202968,-0.000616966,0.000616966,0.00123393,-0.00579446,-0.00223955,-0.00223955,-0.00223955,-0.00223955,-0.00197386,-0.00197386,0.000864451,0.000864451,-0.000864451,-0.000864451,0.0017289,0.0017289,0.00486028,0.00486028,0.00666613,-0.00666613,-0.0144238,-0.00321441,0.00321441,0.00642883,-0.00243289,0.00121644,-0.00121644,-0.00121644,0.00121644,-0.00243289,0.00377325,0.00377325,-0.00201313,0.00201313,0.00402626,-0.000254032,0.000254032,-0.00498991,-0.000580832,0.00222693,-0.00222693,0.00445387,-0.00393928,-0.0074823,-0.0074823,-0.00122566,-0.00122566,0.00122566,0.00122566,-0.00245133,-0.00245133,0.00279609,-0.00279609,-0.00559219,-0.00134955,-0.00134955,0.00134955,0.00134955,-0.0026991,-0.0026991,-0.00329603,-0.00329603,0.00147612,-0.00147612,0.00295225,-0.00285608,-0.00208145,-0.00208145,-0.003203,0.00123752,-0.00123752,0.00247504,-0.00115663,-0.00115663,0.00115663,0.00115663,-0.00231326,-0.00231326,-0.0048064,-0.00302771,
  -0.00302771,0.00261399,0.00261399,0.00522799,-0.00189266,0.00189266,-0.00378532,0.000816396,0.000816396,-0.000816396,-0.000816396,0.00163279,0.00163279,0.00185218,-0.000926088,0.000926088,0.000926088,-0.000926088,0.00185218,-0.00386162,-0.00386162,0.00889961,0.00297926,-0.00297926,-0.00595854,-0.00278364,0.00278364,-0.00556728,-0.000615873,-0.000307937,0.000307937,-0.000307937,0.000307937,0.000615873,0.000615873,-0.000615873,-0.00123175,0.00179692,-0.00179692,0.00359384,-0.00108652,-0.00108652,0.000698212,0.000698212,0.000698212,0.000698212,-0.00100748,-0.00100748,-0.00100748,-0.00100748,0.00246432,0.00246432,-0.00246432,-0.00246432,0.00492865,0.00492865,0.00361321,0.00361321,-0.00387647,0.000850908,-0.000850908,-0.00170182,-0.00245523,0.00274268,0.00274268,-0.00428656,-0.00428656,0.00172362,0.00172362,0.00344724,-0.00098747,0.00098747,-0.00197494,-0.0032867,-0.0018966,0.0018966,-0.0037932,-0.000951686,0.000951686,0.00190337,-0.00382581,0.00194324,0.00194324,0.00388648,0.00177688,-0.00177688,-0.00355376,-0.00338845,-0.00229446,0.00114723,-0.00114723,-0.00114723,0.00114723,-0.00229446,-0.00137346,0.00137346,-0.00274692,-0.0020755,-0.0020755,-0.00684759,
  0.00684759,0.0012913,0.0012913,0.0012913,0.0012913,-0.000769853,-0.000769853,-0.00237076,-0.00237076,-0.00144999,0.00144999,-0.00289998,-0.00185573,-0.00185573,-0.00185573,-0.00185573,0.00172366,-0.00172366,-0.00344732,0.00150922,0.000754614,-0.000754614,0.000754614,-0.000754614,-0.00150922,-0.00150922,0.00150922,0.00301845,0.00228741,0.00228741,-0.00381566,-0.00381566,-0.000784264,-0.000784264,0.000784264,0.000784264,-0.00156853,-0.00156853,0.00112193,-0.00112193,0.00224386,-0.00406245,0.00129594,0.00129594,0.00129594,0.00129594,-0.00105859,-0.00273237,0.00258676,0.00258676,0.00258676,0.00258676,-0.000429814,-0.00187089,-0.00187089,0.00196143,-0.00191085,-0.00191085,-0.00155416,-0.00155416,-0.00155416,-0.00155416,0.00798381,-0.00798381,0.00193561,0.00193561,-0.00181851,0.00181851,0.00363702,0.00126517,0.000632588,-0.000632588,0.000632588,-0.000632588,-0.00126517,-0.00126517,0.00126517,0.00253035,0.000366493,0.000366493,0.0019935,0.0019935,-0.00143629,-0.00491069,-0.00491069,0.0019239,0.0019239,-0.000824457,-0.000824457,-0.00212587,-0.00212587,-0.00425174,-0.00105896,-0.00105896,-0.00105896,-0.00105896,-0.00145224,-0.00145224,-0.00290448,0.00108442,
  0.00054221,-0.00054221,0.00054221,-0.00054221,-0.00108442,-0.00108442,0.00108442,0.00216884,-0.00395638,-0.00395638,0.00204159,0.00204159,-0.00169268,0.00169268,-0.00338537,-0.00116793,-0.00109054,0.00109054,-0.00218108,-0.00197676,-0.00197676,0.000713528,0.000713528,-0.000713528,-0.000713528,0.00142706,0.00142706,0.00108676,0.00108676,0.00217351,0.000382833,0.000382833,-0.000382833,-0.000382833,0.000765666,0.000765666,0.00195272,0.00195272,0.00390544,0.00466503,0.00466503,-0.000810307,-0.000405155,0.000405155,-0.000405155,0.000405155,0.000810307,0.000810307,-0.000810307,-0.00162062,-0.00157283,0.00157283,-0.00314566,0.00101451,-0.00101451,-0.00202903,-0.00119835,-0.00119835,-0.00355643,-0.00355643,-0.00265652,-0.00172965,-0.00172965,0.00185964,0.00185964,-0.00179456,-0.00331868,0.00331868,-0.000298403,-0.000298403,-0.00129323,0.00249165,0.00249165,0.00498331,-0.00234547,-0.00234547,-0.0014806,-0.0014806,-0.00101837,-9.2677e-05,-9.2677e-05,9.2677e-05,9.2677e-05,-0.000185354,-0.000185354,0.000721602,-0.000721602,-0.00144321,0.00108748,0.00108748,0.00275028,0.00275028,0.000323249,0.000323249,0.00257938,0.00257938,-0.00011565,-5.78253e-05,5.78253e-05,-5.78253e-05,
  5.78253e-05,0.00011565,0.00011565,-0.00011565,-0.000231301,-0.0016329,-0.0016329,-0.000426221,0.000426221,0.000852443,-0.00246471,-0.00246471,-0.00492943,0.000630747,0.000630747,-0.000630747,-0.000630747,0.00126149,0.00126149,0.000627024,0.000627024,0.00124132,0.00124132,-0.00129171,-0.00140443,0.000702216,-0.000702216,0.000702216,-0.000702216,0.00140443,-0.00140443,0.00140443,-0.00280886,-0.000240007,-0.000240007,0.000240007,0.000240007,-0.000480014,-0.000480014,-1.58112e-05,-1.58112e-05,1.58112e-05,-1.58112e-05,1.58112e-05,1.58112e-05,-1.58112e-05,1.58112e-05,1.58112e-05,4.74338e-05,-0.000993863,0.000993863,0.00124956,0.00124956,0.00152286,0.00152286,-0.000497582,0.000310118,-0.000155059,0.000155059,-0.000155059,0.000155059,-0.000310118,0.000310118,-0.000310118,0.000620237,-0.0012808,0.00333848,0.00333848,0.00237162,0.00237162,2.4549e-05,1.22746e-05,-1.22746e-05,1.22746e-05,-1.22746e-05,-2.4549e-05,-2.4549e-05,2.4549e-05,4.90981e-05,-0.000710472,0.000710472,-0.00142095,-0.000578698,-0.000578698,-0.000345737,-0.000172869,0.000172869,0.000172869,-0.000172869,-0.000345737,0.00101373,0.00101373,-0.00147656,-0.00147656,-0.00167819,-0.00167819,-0.00358004,0.00358004,0.00716008,0.00308971,
  0.00308971,0.00124624,0.00124624,0.0022323,0.0022323,-0.00356418,-4.45436e-05,-4.45436e-05,4.45436e-05,4.45436e-05,-8.90873e-05,-8.90873e-05,0.00169604,0.00169604,-0.000902524,-0.000902524,-0.000702533,-0.000702533,-0.000654222,-0.000654222,-4.17307e-05,-4.17307e-05
};
  
  RealType mych225[] = 
{
  -0.956538,-0.0982993,0.0932378,-0.0660692,-0.0469631,-0.0241544,0.0241544,-0.0483089,0.00944105,0.00944105,0.0188821,-0.0271198,-0.0271198,-0.0226402,-0.0226402,-0.0193681,0.00782195,0.00782195,-0.00782195,-0.00782195,0.0156439,0.0156439,0.0204492,0.0204492,0.00405354,-0.00405354,0.00810709,-0.019117,-0.019117,-0.0149198,-0.0149198,-0.0113135,-0.0113135,0.0143962,0.0143962,0.0178877,0.0178877,-0.0138716,-0.0116838,-0.0116838,0.00440998,-0.002205,0.002205,-0.002205,0.002205,-0.00440998,0.00440998,-0.00440998,0.00881997,0.00975916,0.00975916,0.0195184,-0.00698603,0.00698603,0.0139721,-0.00667938,0.00667938,-0.0133588,0.0181814,0.0181814,0.00525466,-0.00262734,0.00262734,-0.00262734,0.00262734,-0.00525466,0.00525466,-0.00525466,0.0105093,-0.00761994,-0.00761994,-0.0165496,0.00832048,0.00832048,0.00832048,0.00832048,-0.0117428,-0.0117428,-0.00697655,-0.00697655,-0.0139531,-0.0106699,-0.00753963,-0.00753963,-0.00585657,0.00585657,0.0117132,-0.00643085,-0.00643085,0.0149402,0.0149402,-0.00682099,-0.00682099,-0.0055485,0.0055485,0.00933482,0.00933482,0.00664357,-0.00664357,-0.0132872,-0.0120297,
  0.00204615,0.00204615,-0.00272349,-0.00272349,0.00272349,0.00272349,-0.00544698,-0.00544698,-0.00904162,-0.00906239,-0.00879279,-0.00108791,0.00108791,0.00217582,0.0143817,4.71543e-05,-4.71543e-05,-9.43088e-05,0.00345565,-0.00345565,0.00691132,0.00390266,-0.00195133,0.00195133,0.00195133,-0.00195133,0.00390266,0.014723,0.00395697,0.00395697,0.00791396,-0.0112765,0.00451305,0.00451305,0.00835104,-0.00204872,0.00204872,-0.00409744,0.00325652,-0.00325652,0.00651305,0.00609387,0.00609387,0.0040223,-0.0040223,0.00804461,0.00223775,-0.00111888,0.00111888,-0.00111888,0.00111888,-0.00223775,0.00223775,-0.00223775,0.00447551,-0.00515622,0.00796067,0.00796067,-0.00518275,0.00518275,0.0103655,0.00390874,0.00390874,0.0078175,-0.00155762,0.000778811,-0.000778811,0.000778811,-0.000778811,0.00155762,-0.00155762,0.00155762,-0.00311524,-0.00538144,-0.00573521,-0.00789928,0.00789928,0.00397355,0.00397355,0.0053678,-0.0026839,0.0026839,0.0026839,-0.0026839,0.0053678,-0.00645674,-0.00645674,-0.000863048,0.000863048,0.0017261,0.0132942,0.000998809,0.000998809,-0.000998809,-0.000998809,0.00199762,0.00199762,0.00883238,0.00883238,0.00440105,
  0.00440105,-0.00160499,0.00160499,0.00320999,-0.00791044,-0.00447324,-0.00447324,-0.00894648,0.00729231,-0.00729231,-0.0145846,0.00234982,-0.00234982,0.00469966,-0.000986604,-0.000986604,0.000986604,0.000986604,-0.00197321,-0.00197321,-0.00121438,0.00121438,0.00242877,-0.00703836,-0.00252283,-0.00252283,-0.00252283,-0.00252283,-0.00238538,-0.00238538,0.000727593,0.000727593,-0.000727593,-0.000727593,0.00145519,0.00145519,0.0043847,0.0043847,-0.0129596,0.0129596,-0.0111635,-0.00340405,0.00340405,0.00680811,-0.00275124,0.00137562,-0.00137562,-0.00137562,0.00137562,-0.00275124,0.0030333,0.0030333,-0.00220115,0.00220115,0.0044023,-0.00801632,0.00801632,-0.00450464,-7.42713e-06,0.00168337,-0.00168337,0.00336675,-0.00386412,-0.00409521,-0.00409521,-0.00153862,-0.00153862,0.00153862,0.00153862,-0.00307725,-0.00307725,0.00322213,-0.00322213,-0.00644428,-0.00131935,-0.00131935,0.00131935,0.00131935,-0.0026387,-0.0026387,-0.00168171,-0.00168171,0.00136491,-0.00136491,0.00272982,-0.00920113,-0.00281841,-0.00281841,-0.00139964,0.00126684,-0.00126684,0.00253368,-0.00103482,-0.00103482,0.00103482,0.00103482,-0.00206964,-0.00206964,-0.00555547,-0.00378182,
  -0.00378182,0.00229632,0.00229632,0.00459264,-0.00222752,0.00222752,-0.00445506,0.000679506,0.000679506,-0.000679506,-0.000679506,0.00135901,0.00135901,0.00236429,-0.00118215,0.00118215,0.00118215,-0.00118215,0.00236429,-0.00387231,-0.00387231,0.0108272,0.00178948,-0.00178948,-0.00357897,-0.00286181,0.00286181,-0.00572362,-0.000127332,-6.36663e-05,6.36663e-05,-6.36663e-05,6.36663e-05,0.000127332,0.000127332,-0.000127332,-0.000254665,0.00168771,-0.00168771,0.00337543,-0.00160774,-0.00160774,0.000814653,0.000814653,0.000814653,0.000814653,-0.00102598,-0.00102598,-0.00102598,-0.00102598,0.0025492,0.0025492,-0.0025492,-0.0025492,0.00509841,0.00509841,0.00335635,0.00335635,-0.00549849,0.0019903,-0.0019903,-0.0039806,-0.00299516,0.00216393,0.00216393,-0.00454603,-0.00454603,0.00172211,0.00172211,0.00344423,-0.00113393,0.00113393,-0.00226786,-0.00283492,-0.00180333,0.00180333,-0.00360667,-0.00109754,0.00109754,0.00219508,-0.00324517,0.00205393,0.00205393,0.00410786,0.00172038,-0.00172038,-0.00344077,-0.00237549,-0.00263438,0.00131719,-0.00131719,-0.00131719,0.00131719,-0.00263438,-0.00132607,0.00132607,-0.00265215,-0.00214207,-0.00214207,0.00698467,
  -0.00698467,0.00134273,0.00134273,0.00134273,0.00134273,-0.00148195,-0.00148195,-0.00208378,-0.00208378,-0.00114136,0.00114136,-0.00228273,-0.00129762,-0.00129762,-0.00129762,-0.00129762,0.00185857,-0.00185857,-0.00371714,0.00159345,0.000796726,-0.000796726,0.000796726,-0.000796726,-0.00159345,-0.00159345,0.00159345,0.0031869,0.00210325,0.00210325,-0.00353516,-0.00353516,-0.000681217,-0.000681217,0.000681217,0.000681217,-0.00136243,-0.00136243,0.00105626,-0.00105626,0.00211253,-0.00545501,0.000427365,0.000427365,0.000427365,0.000427365,-0.00135367,-0.0031044,0.00313756,0.00313756,0.00313756,0.00313756,-0.00172289,-0.00239183,-0.00239183,0.00157564,-0.00171957,-0.00171957,-0.00131717,-0.00131717,-0.00131717,-0.00131717,0.00417335,-0.00417335,0.00151895,0.00151895,-0.00123308,0.00123308,0.00246617,0.00121564,0.000607824,-0.000607824,0.000607824,-0.000607824,-0.00121564,-0.00121564,0.00121564,0.00243129,0.000885955,0.000885955,0.001455,0.001455,-0.00210544,-0.00399471,-0.00399471,0.00224013,0.00224013,-0.000439071,-0.000439071,-0.00359822,-0.00359822,-0.00719645,-0.000467738,-0.000467738,-0.000467738,-0.000467738,-0.00154437,-0.00154437,-0.00308873,0.000798269,
  0.000399136,-0.000399136,0.000399136,-0.000399136,-0.000798269,-0.000798269,0.000798269,0.00159654,-0.0032592,-0.0032592,0.00172116,0.00172116,-0.00180338,0.00180338,-0.00360676,-0.0021236,-0.000755811,0.000755811,-0.00151162,-0.00281705,-0.00281705,0.000738928,0.000738928,-0.000738928,-0.000738928,0.00147786,0.00147786,0.000648169,0.000648169,0.00129634,0.000193031,0.000193031,-0.000193031,-0.000193031,0.000386063,0.000386063,0.00212067,0.00212067,0.00424135,0.00470131,0.00470131,-0.000644844,-0.000322423,0.000322423,-0.000322423,0.000322423,0.000644844,0.000644844,-0.000644844,-0.00128969,-0.0012752,0.0012752,-0.00255041,0.000728642,-0.000728642,-0.00145729,-0.00171516,-0.00171516,-0.00336387,-0.00336387,-0.00234446,-0.00156007,-0.00156007,0.00194746,0.00194746,-0.00146673,-0.000539651,0.000539651,0.00053944,0.00053944,-0.00213051,0.00284864,0.00284864,0.00569728,-0.00111807,-0.00111807,-0.00157715,-0.00157715,-0.00110496,0.000527848,0.000527848,-0.000527848,-0.000527848,0.0010557,0.0010557,0.00112775,-0.00112775,-0.00225551,0.00118218,0.00118218,0.00202063,0.00202063,0.00143137,0.00143137,0.00296597,0.00296597,-0.000558604,-0.000279303,0.000279303,-0.000279303,
  0.000279303,0.000558604,0.000558604,-0.000558604,-0.00111721,-0.000955159,-0.000955159,-0.000574601,0.000574601,0.0011492,-0.000835028,-0.000835028,-0.00167006,0.000811382,0.000811382,-0.000811382,-0.000811382,0.00162276,0.00162276,0.000824027,0.000824027,0.000393994,0.000393994,-0.000803782,-0.00145966,0.000729831,-0.000729831,0.000729831,-0.000729831,0.00145966,-0.00145966,0.00145966,-0.00291932,-0.000486959,-0.000486959,0.000486959,0.000486959,-0.000973918,-0.000973918,0.000539715,0.000539715,-0.000539715,0.000539715,-0.000539715,-0.000539715,0.000539715,-0.000539715,-0.000539715,-0.00161915,-0.0027409,0.0027409,0.00178177,0.00178177,0.00177389,0.00177389,-0.00131566,-0.000114294,5.71473e-05,-5.71473e-05,5.71473e-05,-5.71473e-05,0.000114294,-0.000114294,0.000114294,-0.000228589,-0.00141017,0.00309887,0.00309887,0.00426398,0.00426398,-5.819e-05,-2.90951e-05,2.90951e-05,-2.90951e-05,2.90951e-05,5.819e-05,5.819e-05,-5.819e-05,-0.00011638,-0.000557547,0.000557547,-0.0011151,-0.000901317,-0.000901317,0.000213297,0.000106649,-0.000106649,-0.000106649,0.000106649,0.000213297,0.00111729,0.00111729,-0.000226797,-0.000226797,-0.00243008,-0.00243008,-0.00210061,0.00210061,0.00420123,0.00162139,
  0.00162139,0.00235796,0.00235796,0.00241742,0.00241742,-0.00261122,-0.000412234,-0.000412234,0.000412234,0.000412234,-0.000824469,-0.000824469,0.00117833,0.00117833,-0.00142276,-0.00142276,-0.000595436,-0.000595436,-0.000864669,-0.000864669,-0.000512024,-0.000512024
};

  vector<RealType> ch20(mych20,mych20+sizeof(mych20)/sizeof(RealType));
  vector<RealType> ch225(mych225,mych225+sizeof(mych225)/sizeof(RealType));

  RealType Ro=2.116493107;
  RealType RminusRo=R-Ro;

  // C contains current coefficients
  vector<RealType>::iterator it(C.begin()),last(C.end());
  vector<RealType>::iterator ch20_it(ch20.begin()),ch20_last(ch20.end());
  vector<RealType>::iterator ch225_it(ch225.begin()),ch225_last(ch225.end());
  while(it != last){ // linear interpolation
    // Ci(R)= (Ci(2.25)-Ci(2.0))/0.25 * (R-2.25) + Ci(2.25)
    (*it) = (*ch225_it-*ch20_it)/0.25 * (RminusRo+Ro-2.25)+*ch225_it;
    it++; ch20_it++; ch225_it++;
  }
}

OrbitalBase::ValueType MultiSlaterDeterminantFast::evaluate(ParticleSet& P
    , ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
{
  EvaluateTimer.start();
  //for(int i=0; i<Dets.size(); i++)
  //  Dets[i]->evaluateForWalkerMove(P);
  Dets[0]->evaluateForWalkerMove(P);
  Dets[1]->evaluateForWalkerMove(P);
  // can this change over time??? I don't know yet
  ValueVector_t& detValues_up = Dets[0]->detValues;
  ValueVector_t& detValues_dn = Dets[1]->detValues;
  GradMatrix_t& grads_up = Dets[0]->grads;
  GradMatrix_t& grads_dn = Dets[1]->grads;
  ValueMatrix_t& lapls_up = Dets[0]->lapls;
  ValueMatrix_t& lapls_dn = Dets[1]->lapls;
  int N1 = Dets[0]->FirstIndex;
  int N2 = Dets[1]->FirstIndex;
  int NP1 = Dets[0]->NumPtcls;
  int NP2 = Dets[1]->NumPtcls;
  psiCurrent=0.0;
  myG=0.0;
  myL=0.0;
  vector<int>::iterator upC(C2node_up.begin()),dnC(C2node_dn.begin());
  vector<RealType>::iterator it(C.begin()),last(C.end());
  //int idx=0; // !!!!!!!! debuging coefficient interpolation
  while(it != last)
  {
    /*if (idx<3){
      cout << idx << ": " << (*it) << endl;
    }
    idx++;*/
    
    psiCurrent += (*it)*detValues_up[*upC]*detValues_dn[*dnC];
    for(int k=0,n=N1; k<NP1; k++,n++)
    {
      myG(n) += (*it)*grads_up(*upC,k)*detValues_dn[*dnC];
      myL(n) += (*it)*lapls_up(*upC,k)*detValues_dn[*dnC];
    }
    for(int k=0,n=N2; k<NP2; k++,n++)
    {
      myG(n) += (*it)*grads_dn(*dnC,k)*detValues_up[*upC];
      myL(n) += (*it)*lapls_dn(*dnC,k)*detValues_up[*upC];
    }
    it++;
    upC++;
    dnC++;
  }
  ValueType psiinv = 1.0/psiCurrent;
  myG *= psiinv;
  myL *= psiinv;
  G += myG;
  for(int i=0; i<L.size(); i++)
    L(i) += myL[i] - dot(myG[i],myG[i]);
  EvaluateTimer.stop();
  return psiCurrent;
}

OrbitalBase::RealType MultiSlaterDeterminantFast::evaluateLog(ParticleSet& P
    , ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
{
  ValueType psi = evaluate(P,G,L);
  return LogValue = evaluateLogAndPhase(psi,PhaseValue);
}

OrbitalBase::RealType MultiSlaterDeterminantFast::evaluateLog(ParticleSet& P,
    ParticleSet::ParticleGradient_t& G,
    ParticleSet::ParticleLaplacian_t& L,
    PooledData<RealType>& buf,
    bool fillBuffer )
{
  if(fillBuffer)
  {
    Dets[0]->evaluateForWalkerMove(P);
    Dets[0]->copyToDerivativeBuffer(P,buf);
    Dets[1]->evaluateForWalkerMove(P);
    Dets[1]->copyToDerivativeBuffer(P,buf);
  }
  else
  {
    Dets[0]->copyFromDerivativeBuffer(P,buf);
    Dets[1]->copyFromDerivativeBuffer(P,buf);
  }
  // can this change over time??? I don't know yet
  ValueVector_t& detValues_up = Dets[0]->detValues;
  ValueVector_t& detValues_dn = Dets[1]->detValues;
  GradMatrix_t& grads_up = Dets[0]->grads;
  GradMatrix_t& grads_dn = Dets[1]->grads;
  ValueMatrix_t& lapls_up = Dets[0]->lapls;
  ValueMatrix_t& lapls_dn = Dets[1]->lapls;
  int N1 = Dets[0]->FirstIndex;
  int N2 = Dets[1]->FirstIndex;
  int NP1 = Dets[0]->NumPtcls;
  int NP2 = Dets[1]->NumPtcls;
  psiCurrent=0.0;
  myG=0.0;
  myL=0.0;
  vector<int>::iterator upC(C2node_up.begin()),dnC(C2node_dn.begin());
  vector<RealType>::iterator it(C.begin()),last(C.end());
  while(it != last)
  {
    psiCurrent += (*it)*detValues_up[*upC]*detValues_dn[*dnC];
    for(int k=0,n=N1; k<NP1; k++,n++)
    {
      myG(n) += (*it)*grads_up(*upC,k)*detValues_dn[*dnC];
      myL(n) += (*it)*lapls_up(*upC,k)*detValues_dn[*dnC];
    }
    for(int k=0,n=N2; k<NP2; k++,n++)
    {
      myG(n) += (*it)*grads_dn(*dnC,k)*detValues_up[*upC];
      myL(n) += (*it)*lapls_dn(*dnC,k)*detValues_up[*upC];
    }
    it++;
    upC++;
    dnC++;
  }
  ValueType psiinv = 1.0/psiCurrent;
  myG *= psiinv;
  myL *= psiinv;
  G += myG;
  for(int i=0; i<L.size(); i++)
    L(i) += myL[i] - dot(myG[i],myG[i]);
  return evaluateLogAndPhase(psiCurrent,PhaseValue);
}


OrbitalBase::GradType MultiSlaterDeterminantFast::evalGrad(ParticleSet& P, int iat)
{
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: evalGrad not implemented. \n");
  }
  GradType grad_iat;
  if(DetID[iat] == 0)
  {
    Dets[0]->evaluateGrads(P,iat);
    ValueVector_t& detValues_up = Dets[0]->detValues;
    ValueVector_t& detValues_dn = Dets[1]->detValues;
    GradMatrix_t& grads_up = Dets[0]->grads;
    int N1 = Dets[0]->FirstIndex;
    ValueType psi=0.0;
    vector<int>::iterator upC(C2node_up.begin()),dnC(C2node_dn.begin());
    vector<RealType>::iterator it(C.begin()),last(C.end());
    while(it != last)
    {
      psi += (*it)*detValues_up[*upC]*detValues_dn[*dnC];
      grad_iat += (*it)*grads_up(*upC,iat-N1)*detValues_dn[*dnC];
      it++;
      upC++;
      dnC++;
    }
    grad_iat *= 1.0/psi;
    return grad_iat;
  }
  else
  {
    Dets[1]->evaluateGrads(P,iat);
    ValueType psi=0.0;
    ValueVector_t& detValues_up = Dets[0]->detValues;
    ValueVector_t& detValues_dn = Dets[1]->detValues;
    GradMatrix_t& grads_dn = Dets[1]->grads;
    int N2 = Dets[1]->FirstIndex;
    vector<int>::iterator upC(C2node_up.begin()),dnC(C2node_dn.begin());
    vector<RealType>::iterator it(C.begin()),last(C.end());
    while(it != last)
    {
      psi += (*it)*detValues_up[*upC]*detValues_dn[*dnC];
      grad_iat += (*it)*grads_dn(*dnC,iat-N2)*detValues_up[*upC];
      it++;
      upC++;
      dnC++;
    }
    grad_iat *= 1.0/psi;
    return grad_iat;
  }
}

OrbitalBase::ValueType MultiSlaterDeterminantFast::ratioGrad(ParticleSet& P
    , int iat, GradType& grad_iat)
{
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: ratioGrad not implemented. \n");
  }
  UpdateMode=ORB_PBYP_PARTIAL;
  if(DetID[iat] == 0)
  {
    RatioGradTimer.start();
    Ratio1GradTimer.start();
    Dets[0]->evaluateDetsAndGradsForPtclMove(P,iat);
    Ratio1GradTimer.stop();
    ValueVector_t& detValues_up = Dets[0]->new_detValues;
    ValueVector_t& detValues_dn = Dets[1]->detValues;
    GradMatrix_t& grads_up = Dets[0]->new_grads;
    int N1 = Dets[0]->FirstIndex;
    vector<int>::iterator upC(C2node_up.begin()),dnC(C2node_dn.begin());
    vector<RealType>::iterator it(C.begin()),last(C.end());
    ValueType psiNew=0.0;
    GradType dummy;
    it=C.begin();
    last=C.end();
    while(it != last)
    {
      psiNew += (*it)*detValues_up[*upC]*detValues_dn[*dnC];
      dummy += (*it)*grads_up(*upC,iat-N1)*detValues_dn[*dnC];
      it++;
      upC++;
      dnC++;
    }
    grad_iat+=dummy/psiNew;
    curRatio = psiNew/psiCurrent;
    RatioGradTimer.stop();
    return curRatio;
  }
  else
  {
    RatioGradTimer.start();
    Ratio1GradTimer.start();
    Dets[1]->evaluateDetsAndGradsForPtclMove(P,iat);
    Ratio1GradTimer.stop();
    ValueVector_t& detValues_up = Dets[0]->detValues;
    ValueVector_t& detValues_dn = Dets[1]->new_detValues;
    GradMatrix_t& grads_dn = Dets[1]->new_grads;
    int N2 = Dets[1]->FirstIndex;
    vector<int>::iterator upC(C2node_up.begin()),dnC(C2node_dn.begin());
    vector<RealType>::iterator it(C.begin()),last(C.end());
    ValueType psiNew=0.0;
    GradType dummy;
    while(it != last)
    {
      psiNew += (*it)*detValues_up[*upC]*detValues_dn[*dnC];
      dummy += (*it)*grads_dn(*dnC,iat-N2)*detValues_up[*upC];
      it++;
      upC++;
      dnC++;
    }
    grad_iat+=dummy/psiNew;
    curRatio = psiNew/psiCurrent;
    RatioGradTimer.stop();
    return curRatio;
  }
}


// This routine need work, sloppy for now
OrbitalBase::ValueType  MultiSlaterDeterminantFast::ratio(ParticleSet& P, int iat
    , ParticleSet::ParticleGradient_t& dG,ParticleSet::ParticleLaplacian_t& dL)
{
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: ratio(P,dG,dL) not implemented. \n");
  }
  UpdateMode=ORB_PBYP_ALL;
  if(DetID[iat] == 0)
  {
    RatioAllTimer.start();
    /*
          P.acceptMove(iat);
          Dets[0]->evaluateForWalkerMove(P);
          ValueVector_t& detValues_up = Dets[0]->detValues;
          ValueVector_t& detValues_dn = Dets[1]->detValues;
          GradMatrix_t& grads_up = Dets[0]->grads;
          GradMatrix_t& grads_dn = Dets[1]->grads;
          ValueMatrix_t& lapls_up = Dets[0]->lapls;
          ValueMatrix_t& lapls_dn = Dets[1]->lapls;
    */
//*
    Ratio1AllTimer.start();
    Dets[0]->evaluateAllForPtclMove(P,iat);
    Ratio1AllTimer.stop();
    ValueVector_t& detValues_up = Dets[0]->new_detValues;
    ValueVector_t& detValues_dn = Dets[1]->detValues;
    GradMatrix_t& grads_up = Dets[0]->new_grads;
    GradMatrix_t& grads_dn = Dets[1]->grads;
    ValueMatrix_t& lapls_up = Dets[0]->new_lapls;
    ValueMatrix_t& lapls_dn = Dets[1]->lapls;
//*/
    int N1 = Dets[0]->FirstIndex;
    int N2 = Dets[1]->FirstIndex;
    int NP1 = Dets[0]->NumPtcls;
    int NP2 = Dets[1]->NumPtcls;
    ValueType psiNew=0.0;
    // myG,myL should contain current grad and lapl
    vector<int>::iterator upC(C2node_up.begin()),dnC(C2node_dn.begin());
    vector<RealType>::iterator it(C.begin()),last(C.end());
    myG_temp=0.0;
    myL_temp=0.0;
    while(it != last)
    {
      psiNew += (*it)*detValues_up[*upC]*detValues_dn[*dnC];
      for(int k=0,n=N1; k<NP1; k++,n++)
      {
        myG_temp(n) += (*it)*grads_up(*upC,k)*detValues_dn[*dnC];
        myL_temp(n) += (*it)*lapls_up(*upC,k)*detValues_dn[*dnC];
      }
      for(int k=0,n=N2; k<NP2; k++,n++)
      {
        myG_temp(n) += (*it)*grads_dn(*dnC,k)*detValues_up[*upC];
        myL_temp(n) += (*it)*lapls_dn(*dnC,k)*detValues_up[*upC];
      }
      it++;
      upC++;
      dnC++;
    }
    ValueType psiNinv=1.0/psiNew;
    myG_temp *= psiNinv;
    myL_temp *= psiNinv;
    dG += myG_temp-myG;
    for(int i=0; i<dL.size(); i++)
      dL(i) += myL_temp[i] - myL[i] - dot(myG_temp[i],myG_temp[i]) + dot(myG[i],myG[i]);
    curRatio = psiNew/psiCurrent;
    RatioAllTimer.stop();
    return curRatio;
  }
  else
  {
    RatioAllTimer.start();
    Ratio1AllTimer.start();
    Dets[1]->evaluateAllForPtclMove(P,iat);
    Ratio1AllTimer.stop();
    ValueVector_t& detValues_up = Dets[0]->detValues;
    ValueVector_t& detValues_dn = Dets[1]->new_detValues;
    GradMatrix_t& grads_up = Dets[0]->grads;
    GradMatrix_t& grads_dn = Dets[1]->new_grads;
    ValueMatrix_t& lapls_up = Dets[0]->lapls;
    ValueMatrix_t& lapls_dn = Dets[1]->new_lapls;
    int N1 = Dets[0]->FirstIndex;
    int N2 = Dets[1]->FirstIndex;
    int NP1 = Dets[0]->NumPtcls;
    int NP2 = Dets[1]->NumPtcls;
    ValueType psiNew=0.0;
    // myG,myL should contain current grad and lapl
    vector<int>::iterator upC(C2node_up.begin()),dnC(C2node_dn.begin());
    vector<RealType>::iterator it(C.begin()),last(C.end());
    myG_temp=0.0;
    myL_temp=0.0;
    while(it != last)
    {
      psiNew += (*it)*detValues_up[*upC]*detValues_dn[*dnC];
      for(int k=0,n=N1; k<NP1; k++,n++)
      {
        myG_temp(n) += (*it)*grads_up(*upC,k)*detValues_dn[*dnC];
        myL_temp(n) += (*it)*lapls_up(*upC,k)*detValues_dn[*dnC];
      }
      for(int k=0,n=N2; k<NP2; k++,n++)
      {
        myG_temp(n) += (*it)*grads_dn(*dnC,k)*detValues_up[*upC];
        myL_temp(n) += (*it)*lapls_dn(*dnC,k)*detValues_up[*upC];
      }
      it++;
      upC++;
      dnC++;
    }
    ValueType psiNinv=1.0/psiNew;
    myG_temp *= psiNinv;
    myL_temp *= psiNinv;
    dG += myG_temp-myG;
    for(int i=0; i<dL.size(); i++)
      dL(i) += myL_temp[i] - myL[i] - dot(myG_temp[i],myG_temp[i]) + dot(myG[i],myG[i]);
    curRatio = psiNew/psiCurrent;
    RatioAllTimer.stop();
    return curRatio;
  }
}

// use ci_node for this routine only
OrbitalBase::ValueType MultiSlaterDeterminantFast::ratio(ParticleSet& P, int iat)
{
// debug
//    testMSD(P,iat);
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: ratio not implemented. \n");
  }
  UpdateMode=ORB_PBYP_RATIO;
  if(DetID[iat] == 0)
  {
    RatioTimer.start();
    Ratio1Timer.start();
    Dets[0]->evaluateDetsForPtclMove(P,iat);
    Ratio1Timer.stop();
    ValueVector_t& detValues_up = Dets[0]->new_detValues;
    ValueVector_t& detValues_dn = Dets[1]->detValues;
    ValueType psiNew=0.0;
    vector<int>::iterator upC(C2node_up.begin()),dnC(C2node_dn.begin());
    vector<RealType>::iterator it(C.begin()),last(C.end());
    while(it != last)
      psiNew += (*(it++))*detValues_up[*(upC++)]*detValues_dn[*(dnC++)];
    curRatio = psiNew/psiCurrent;
    RatioTimer.stop();
    return curRatio;
  }
  else
  {
    RatioTimer.start();
    Ratio1Timer.start();
    Dets[1]->evaluateDetsForPtclMove(P,iat);
    Ratio1Timer.stop();
    ValueVector_t& detValues_up = Dets[0]->detValues;
    ValueVector_t& detValues_dn = Dets[1]->new_detValues;
    ValueType psiNew=0.0;
    vector<int>::iterator upC(C2node_up.begin()),dnC(C2node_dn.begin());
    vector<RealType>::iterator it(C.begin()),last(C.end());
    while(it != last)
      psiNew += (*(it++))*detValues_up[*(upC++)]*detValues_dn[*(dnC++)];
    curRatio = psiNew/psiCurrent;
    RatioTimer.stop();
    return curRatio;
  }
}

void MultiSlaterDeterminantFast::acceptMove(ParticleSet& P, int iat)
{
// this should depend on the type of update, ratio / ratioGrad
// for now is incorrect fot ratio(P,iat,dG,dL) updates
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: acceptMove not implemented. \n");
  }
// update psiCurrent,myG_temp,myL_temp
  AccRejTimer.start();
  psiCurrent *= curRatio;
  curRatio=1.0;
  Dets[DetID[iat]]->acceptMove(P,iat);
  switch(UpdateMode)
  {
  case ORB_PBYP_ALL:
    // ratio(P,iat,dG,dL)
    myG = myG_temp;
    myL = myL_temp;
    break;
  default:
    break;
  }
  AccRejTimer.stop();
//    Dets[0]->evaluateForWalkerMove(P);
//    Dets[1]->evaluateForWalkerMove(P);
  // can this change over time??? I don't know yet
  /*
      ValueVector_t& detValues_up = Dets[0]->detValues;
      ValueVector_t& detValues_dn = Dets[1]->detValues;
      GradMatrix_t& grads_up = Dets[0]->grads;
      GradMatrix_t& grads_dn = Dets[1]->grads;
      ValueMatrix_t& lapls_up = Dets[0]->lapls;
      ValueMatrix_t& lapls_dn = Dets[1]->lapls;
      int N1 = Dets[0]->FirstIndex;
      int N2 = Dets[1]->FirstIndex;
      int NP1 = Dets[0]->NumPtcls;
      int NP2 = Dets[1]->NumPtcls;

      ValueType psi=0.0;
      myG_temp=0.0;
      myL_temp=0.0;
      vector<int>::iterator upC(C2node_up.begin()),dnC(C2node_dn.begin());
      vector<RealType>::iterator it(C.begin()),last(C.end());
      while(it != last) {
        psi += (*it)*detValues_up[*upC]*detValues_dn[*dnC];
        for(int k=0,n=N1; k<NP1; k++,n++) {
          myG_temp(n) += (*it)*grads_up(*upC,k)*detValues_dn[*dnC];
          myL_temp(n) += (*it)*lapls_up(*upC,k)*detValues_dn[*dnC];
        }
        for(int k=0,n=N2; k<NP2; k++,n++) {
          myG_temp(n) += (*it)*grads_dn(*dnC,k)*detValues_up[*upC];
          myL_temp(n) += (*it)*lapls_dn(*dnC,k)*detValues_up[*upC];
        }
        it++;upC++;dnC++;
      }
      ValueType psiinv = 1.0/psi;
      myG_temp *= psiinv;
      myL_temp *= psiinv;
  */
}

void MultiSlaterDeterminantFast::restore(int iat)
{
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: restore not implemented. \n");
  }
  AccRejTimer.start();
  Dets[DetID[iat]]->restore(iat);
  curRatio=1.0;
  AccRejTimer.stop();
}

void MultiSlaterDeterminantFast::update(ParticleSet& P
                                        , ParticleSet::ParticleGradient_t& dG, ParticleSet::ParticleLaplacian_t& dL
                                        , int iat)
{
  APP_ABORT("IMPLEMENT MultiSlaterDeterminantFast::update");
}

OrbitalBase::RealType MultiSlaterDeterminantFast::evaluateLog(ParticleSet& P,BufferType& buf)
{
  Dets[0]->evaluateLog(P,buf);
  Dets[1]->evaluateLog(P,buf);
  buf.put(psiCurrent);
  buf.put(myL.first_address(), myL.last_address());
  buf.put(FirstAddressOfG,LastAddressOfG);
  return LogValue = evaluateLogAndPhase(psiCurrent,PhaseValue);
}

OrbitalBase::RealType MultiSlaterDeterminantFast::registerData(ParticleSet& P, BufferType& buf)
{
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: restore not implemented. \n");
  }
  Dets[0]->registerData(P,buf);
  Dets[1]->registerData(P,buf);
  LogValue = evaluateLog(P,P.G,P.L);
  FirstAddressOfG = &myG[0][0];
  LastAddressOfG = FirstAddressOfG + P.getTotalNum()*DIM;
  buf.add(psiCurrent);
  buf.add(myL.first_address(), myL.last_address());
  buf.add(FirstAddressOfG,LastAddressOfG);
// debug, erase
//    msd->registerData(P,buf);
  return LogValue;
}

// this routine does not initialize the data, just reserves the space
void MultiSlaterDeterminantFast::registerDataForDerivatives(ParticleSet& P, BufferType& buf, int storageType)
{
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: registerDataForDerivatives not implemented. \n");
  }
  Dets[0]->registerDataForDerivatives(P,buf,storageType);
  Dets[1]->registerDataForDerivatives(P,buf,storageType);
}

// FIX FIX FIX
OrbitalBase::RealType MultiSlaterDeterminantFast::updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch)
{
  UpdateTimer.start();
  Dets[0]->updateBuffer(P,buf,fromscratch);
  Dets[1]->updateBuffer(P,buf,fromscratch);
  //Dets[0]->updateBuffer(P,buf,true);
  //Dets[1]->updateBuffer(P,buf,true);
  // can this change over time??? I don't know yet
  ValueVector_t& detValues_up = Dets[0]->detValues;
  ValueVector_t& detValues_dn = Dets[1]->detValues;
  GradMatrix_t& grads_up = Dets[0]->grads;
  GradMatrix_t& grads_dn = Dets[1]->grads;
  ValueMatrix_t& lapls_up = Dets[0]->lapls;
  ValueMatrix_t& lapls_dn = Dets[1]->lapls;
  int N1 = Dets[0]->FirstIndex;
  int N2 = Dets[1]->FirstIndex;
  int NP1 = Dets[0]->NumPtcls;
  int NP2 = Dets[1]->NumPtcls;
  psiCurrent=0.0;
  myG=0.0;
  myL=0.0;
  vector<int>::iterator upC(C2node_up.begin()),dnC(C2node_dn.begin());
  vector<RealType>::iterator it(C.begin()),last(C.end());
  while(it != last)
  {
    psiCurrent += (*it)*detValues_up[*upC]*detValues_dn[*dnC];
    for(int k=0,n=N1; k<NP1; k++,n++)
    {
      myG(n) += (*it)*grads_up(*upC,k)*detValues_dn[*dnC];
      myL(n) += (*it)*lapls_up(*upC,k)*detValues_dn[*dnC];
    }
    for(int k=0,n=N2; k<NP2; k++,n++)
    {
      myG(n) += (*it)*grads_dn(*dnC,k)*detValues_up[*upC];
      myL(n) += (*it)*lapls_dn(*dnC,k)*detValues_up[*upC];
    }
    it++;
    upC++;
    dnC++;
  }
  ValueType psiinv = 1.0/psiCurrent;
  myG *= psiinv;
  myL *= psiinv;
  P.G += myG;
  for(int i=0; i<P.L.size(); i++)
    P.L(i) += myL[i] - dot(myG[i],myG[i]);
  buf.put(psiCurrent);
  buf.put(myL.first_address(), myL.last_address());
  buf.put(FirstAddressOfG,LastAddressOfG);
  UpdateTimer.stop();
  return LogValue = evaluateLogAndPhase(psiCurrent,PhaseValue);
}

void MultiSlaterDeterminantFast::copyFromBuffer(ParticleSet& P, BufferType& buf)
{
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: copyFromBuffer not implemented. \n");
  }
  Dets[0]->copyFromBuffer(P,buf);
  Dets[1]->copyFromBuffer(P,buf);
  buf.get(psiCurrent);
  buf.get(myL.first_address(), myL.last_address());
  buf.get(FirstAddressOfG,LastAddressOfG);
}


void MultiSlaterDeterminantFast::checkInVariables(opt_variables_type& active)
{
  if(Optimizable)
  {
    if(myVars.size())
      active.insertFrom(myVars);
    else
      Optimizable=false;
  }
}

void MultiSlaterDeterminantFast::checkOutVariables(const opt_variables_type& active)
{
  if(Optimizable)
    myVars.getIndex(active);
}

/** resetParameters with optVariables
 *
 * USE_resetParameters
 */
void MultiSlaterDeterminantFast::resetParameters(const opt_variables_type& active)
{
  if(Optimizable)
  {
    if(usingCSF)
    {
      for(int i=0; i<CSFcoeff.size()-1; i++)
      {
        int loc=myVars.where(i);
        if(loc>=0)
          CSFcoeff[i+1]=myVars[i]=active[loc];
      }
      int cnt=0;
      for(int i=0; i<DetsPerCSF.size(); i++)
      {
        for(int k=0; k<DetsPerCSF[i]; k++)
        {
          C[cnt] = CSFcoeff[i]*CSFexpansion[cnt];
          cnt++;
        }
      }
      //for(int i=0; i<Dets.size(); i++) Dets[i]->resetParameters(active);
    }
    else
    {
      for(int i=0; i<C.size()-1; i++)
      {
        int loc=myVars.where(i);
        if(loc>=0)
          C[i+1]=myVars[i]=active[loc];
      }
      //for(int i=0; i<Dets.size(); i++) Dets[i]->resetParameters(active);
    }
  }
}
void MultiSlaterDeterminantFast::reportStatus(ostream& os)
{
}

//   OrbitalBasePtr MultiSlaterDeterminantFast::makeClone(ParticleSet& tqp) const
//   {
//      APP_ABORT("IMPLEMENT OrbitalBase::makeClone");
//      return 0;
//   }

void MultiSlaterDeterminantFast::evaluateDerivatives(ParticleSet& P,
    const opt_variables_type& optvars,
    vector<RealType>& dlogpsi,
    vector<RealType>& dhpsioverpsi)
{
  bool recalculate(false);
  for (int k=0; k<myVars.size(); ++k)
  {
    int kk=myVars.where(k);
    if (kk<0)
      continue;
    if (optvars.recompute(kk))
      recalculate=true;
  }
// need to modify for CSF later on, right now assume Slater Det basis
  if (recalculate)
  {
    if(usingCSF)
    {
      if(laplSum_up.size() == 0)
        laplSum_up.resize(Dets[0]->detValues.size());
      if(laplSum_dn.size() == 0)
        laplSum_dn.resize(Dets[1]->detValues.size());
      // assume that evaluateLog has been called in opt routine before
      //   Dets[0]->evaluateForWalkerMove(P);
      //   Dets[1]->evaluateForWalkerMove(P);
      ValueVector_t& detValues_up = Dets[0]->detValues;
      ValueVector_t& detValues_dn = Dets[1]->detValues;
      GradMatrix_t& grads_up = Dets[0]->grads;
      GradMatrix_t& grads_dn = Dets[1]->grads;
      ValueMatrix_t& lapls_up = Dets[0]->lapls;
      ValueMatrix_t& lapls_dn = Dets[1]->lapls;
      int N1 = Dets[0]->FirstIndex;
      int N2 = Dets[1]->FirstIndex;
      int NP1 = Dets[0]->NumPtcls;
      int NP2 = Dets[1]->NumPtcls;
// myG,myL should already be calculated
      int n = P.getTotalNum();
      ValueType psiinv = 1.0/psiCurrent;
      ValueType lapl_sum=0.0;
      ValueType gg=0.0, ggP=0.0;
      myG_temp=0.0;
      int num=laplSum_up.size();
      ValueVector_t::iterator it(laplSum_up.begin());
      ValueVector_t::iterator last(laplSum_up.end());
      ValueType* ptr0 = lapls_up[0];
      while(it != last)
      {
        (*it)=0.0;
        for(int k=0; k<nels_up; k++,ptr0++)
          (*it) += *ptr0;
        it++;
      }
      it=laplSum_dn.begin();
      last=laplSum_dn.end();
      ptr0 = lapls_dn[0];
      while(it != last)
      {
        (*it)=0.0;
        for(int k=0; k<nels_dn; k++,ptr0++)
          (*it) += *ptr0;
        it++;
      }
      for(int i=0; i<C.size(); i++)
      {
        int upC = C2node_up[i];
        int dnC = C2node_dn[i];
        ValueType tmp1 = C[i]*detValues_dn[dnC]*psiinv;
        ValueType tmp2 = C[i]*detValues_up[upC]*psiinv;
        lapl_sum += tmp1*laplSum_up[upC]+tmp2*laplSum_dn[dnC];
        for(int k=0,j=N1; k<NP1; k++,j++)
          myG_temp[j] += tmp1*grads_up(upC,k);
        for(int k=0,j=N2; k<NP2; k++,j++)
          myG_temp[j] += tmp2*grads_dn(dnC,k);
      }
      gg=ggP=0.0;
      for(int i=0; i<n; i++)
      {
        gg += dot(myG_temp[i],myG_temp[i])-dot(P.G[i],myG_temp[i]);
      }
//       for(int i=0; i<C.size(); i++){
      num=CSFcoeff.size()-1;
      int cnt=0;
//        this one is not optable
      cnt+=DetsPerCSF[0];
      int ip(1);
      for(int i=0; i<num; i++,ip++)
      {
        int kk=myVars.where(i);
        if (kk<0)
        {
          cnt+=DetsPerCSF[ip];
          continue;
        }
        ValueType cdet=0.0,q0=0.0,v1=0.0,v2=0.0;
        for(int k=0; k<DetsPerCSF[ip]; k++)
        {
          int upC = C2node_up[cnt];
          int dnC = C2node_dn[cnt];
          ValueType tmp1=CSFexpansion[cnt]*detValues_dn[dnC]*psiinv;
          ValueType tmp2=CSFexpansion[cnt]*detValues_up[upC]*psiinv;
          cdet+=CSFexpansion[cnt]*detValues_up[upC]*detValues_dn[dnC]*psiinv;
          q0 += (tmp1*laplSum_up[upC] + tmp2*laplSum_dn[dnC]);
          for(int l=0,j=N1; l<NP1; l++,j++)
            v1 += tmp1*(dot(P.G[j],grads_up(upC,l))-dot(myG_temp[j],grads_up(upC,l)) );
          for(int l=0,j=N2; l<NP2; l++,j++)
            v2 += tmp2*(dot(P.G[j],grads_dn(dnC,l))-dot(myG_temp[j],grads_dn(dnC,l)));
          cnt++;
        }
        convert(cdet,dlogpsi[kk]);
        ValueType dhpsi =  -0.5*(q0-cdet*lapl_sum)
                           -cdet*gg-v1-v2;
        //ValueType dhpsi =  -0.5*(tmp1*laplSum_up[upC]+tmp2*laplSum_dn[dnC]
        //                         -cdet*lapl_sum)
        //                   -cdet*gg-(tmp1*v1+tmp2*v2);
        convert(dhpsi,dhpsioverpsi[kk]);
      }
    }
    else
      //usingCSF
    {
      if(laplSum_up.size() == 0)
        laplSum_up.resize(Dets[0]->detValues.size());
      if(laplSum_dn.size() == 0)
        laplSum_dn.resize(Dets[1]->detValues.size());
      // assume that evaluateLog has been called in opt routine before
      //   Dets[0]->evaluateForWalkerMove(P);
      //   Dets[1]->evaluateForWalkerMove(P);
      ValueVector_t& detValues_up = Dets[0]->detValues;
      ValueVector_t& detValues_dn = Dets[1]->detValues;
      GradMatrix_t& grads_up = Dets[0]->grads;
      GradMatrix_t& grads_dn = Dets[1]->grads;
      ValueMatrix_t& lapls_up = Dets[0]->lapls;
      ValueMatrix_t& lapls_dn = Dets[1]->lapls;
      int N1 = Dets[0]->FirstIndex;
      int N2 = Dets[1]->FirstIndex;
      int NP1 = Dets[0]->NumPtcls;
      int NP2 = Dets[1]->NumPtcls;
      int n = P.getTotalNum();
      ValueType psiinv = 1.0/psiCurrent;
      ValueType lapl_sum=0.0;
      ValueType gg=0.0, ggP=0.0;
      myG_temp=0.0;
      int num=laplSum_up.size();
      ValueVector_t::iterator it(laplSum_up.begin());
      ValueVector_t::iterator last(laplSum_up.end());
      ValueType* ptr0 = lapls_up[0];
      while(it != last)
      {
        (*it)=0.0;
        for(int k=0; k<nels_up; k++,ptr0++)
          (*it) += *ptr0;
        it++;
      }
      it=laplSum_dn.begin();
      last=laplSum_dn.end();
      ptr0 = lapls_dn[0];
      while(it != last)
      {
        (*it)=0.0;
        for(int k=0; k<nels_dn; k++,ptr0++)
          (*it) += *ptr0;
        it++;
      }
      for(int i=0; i<C.size(); i++)
      {
        int upC = C2node_up[i];
        int dnC = C2node_dn[i];
        ValueType tmp1 = C[i]*detValues_dn[dnC]*psiinv;
        ValueType tmp2 = C[i]*detValues_up[upC]*psiinv;
        lapl_sum += tmp1*laplSum_up[upC]+tmp2*laplSum_dn[dnC];
        for(int k=0,j=N1; k<NP1; k++,j++)
          myG_temp[j] += tmp1*grads_up(upC,k);
        for(int k=0,j=N2; k<NP2; k++,j++)
          myG_temp[j] += tmp2*grads_dn(dnC,k);
      }
      gg=ggP=0.0;
      for(int i=0; i<n; i++)
      {
        gg += dot(myG_temp[i],myG_temp[i])-dot(P.G[i],myG_temp[i]);
      }
      for(int i=1; i<C.size(); i++)
      {
        int kk=myVars.where(i-1);
        if (kk<0)
          continue;
        int upC = C2node_up[i];
        int dnC = C2node_dn[i];
        ValueType cdet=detValues_up[upC]*detValues_dn[dnC]*psiinv;
        ValueType tmp1=detValues_dn[dnC]*psiinv;
        ValueType tmp2=detValues_up[upC]*psiinv;
        convert(cdet,dlogpsi[kk]);
        ValueType v1=0.0,v2=0.0;
        for(int k=0,j=N1; k<NP1; k++,j++)
          v1 += (dot(P.G[j],grads_up(upC,k))-dot(myG_temp[j],grads_up(upC,k)) );
        for(int k=0,j=N2; k<NP2; k++,j++)
          v2 += (dot(P.G[j],grads_dn(dnC,k))-dot(myG_temp[j],grads_dn(dnC,k)));
        ValueType dhpsi =  -0.5*(tmp1*laplSum_up[upC]+tmp2*laplSum_dn[dnC]
                                 -cdet*lapl_sum)
                           -cdet*gg-(tmp1*v1+tmp2*v2);
        convert(dhpsi,dhpsioverpsi[kk]);
      }
    } // usingCSF
  }
}

void MultiSlaterDeterminantFast::registerTimers()
{
  RatioTimer.reset();
  RatioGradTimer.reset();
  RatioAllTimer.reset();
  Ratio1Timer.reset();
  Ratio1GradTimer.reset();
  Ratio1AllTimer.reset();
  UpdateTimer.reset();
  EvaluateTimer.reset();
  AccRejTimer.reset();
  TimerManager.addTimer (&RatioTimer);
  TimerManager.addTimer (&RatioGradTimer);
  TimerManager.addTimer (&RatioAllTimer);
  TimerManager.addTimer (&Ratio1Timer);
  TimerManager.addTimer (&Ratio1GradTimer);
  TimerManager.addTimer (&Ratio1AllTimer);
  TimerManager.addTimer (&UpdateTimer);
  TimerManager.addTimer (&EvaluateTimer);
  TimerManager.addTimer (&AccRejTimer);
}


}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 3416 $   $Date: 2008-12-07 11:34:49 -0600 (Sun, 07 Dec 2008) $
 * $Id: MultiSlaterDeterminantFast.cpp 3416 2008-12-07 17:34:49Z jnkim $
 ***************************************************************************/
