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
/** @file QMCUpdateBase
 * @brief Declare QMCUpdateBase class
 */
#ifndef QMCPLUSPLUS_QMCUPDATE_BASE_H
#define QMCPLUSPLUS_QMCUPDATE_BASE_H
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/ParticleSet.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "QMCHamiltonians/NonLocalTOperator.h"
#include "QMCDrivers/SimpleFixedNodeBranch.h"
//#define ENABLE_COMPOSITE_ESTIMATOR
//#include "Estimators/CompositeEstimators.h"
#include "Estimators/EstimatorManager.h"

namespace qmcplusplus
{


/** @ingroup QMC
 * @brief Base class for update methods for each step
 *
 * QMCUpdateBase provides the common functions to update all the walkers for each time step.
 * Derived classes should implement advanceWalkers to complete a step.
 */
class QMCUpdateBase: public QMCTraits
{

public:

  typedef MCWalkerConfiguration::Walker_t Walker_t;
  typedef MCWalkerConfiguration::iterator WalkerIter_t;
  typedef SimpleFixedNodeBranch           BranchEngineType;

  ///If true, terminate the simulation
  bool BadState;
  ///number of steps per measurement
  int nSubSteps;
  ///MaxAge>0 indicates branch is done
  IndexType MaxAge;
  ///counter for number of moves accepted
  IndexType nAccept;
  ///counter for number of moves rejected
  IndexType nReject;
  ///Total number of the steps when all the particle moves are rejected.
  IndexType nAllRejected;
  ///Total number of node crossings per block
  IndexType nNodeCrossing;
  ///Total numer of non-local moves accepted
  IndexType NonLocalMoveAccepted;
  ///timestep
  RealType Tau;

  int NumIons;

  ///Amount of steps
  int tauCount;
  ///How frequently the ions are moved ( min(tauCountFreq,nStepsPerBlock) )       
  vector<int> tauCountFreq;
  int nStepsPerBlock;
  // Frequency of all particles is equal to one
  //bool FrequencyEqOne;
  // a list of coefficients for each ion to customize the ion wave function
  vector< vector<RealType> > ionCoeff;

  vector<RealType> nuclei_dist_step;
  vector<RealType> ion_rc;
  ParticleSet::ParticlePos_t ionRo;

  RealType UniformGrid_h;
  RealType nucleiCoeff;

  /// Constructor.
  //QMCUpdateBase(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
  //              RandomGenerator_t& rg);

  QMCUpdateBase(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
                RandomGenerator_t& rg);


  ///Alt Constructor.
  
  QMCUpdateBase(MCWalkerConfiguration& w, TrialWaveFunction& psi, TrialWaveFunction& guide, QMCHamiltonian& h,
                RandomGenerator_t& rg);
  
  ///destructor
  virtual ~QMCUpdateBase();

  inline RealType acceptRatio() const
  {
    return static_cast<RealType>(nAccept)/static_cast<RealType>(nAccept+nReject);
  }

  /** reset the QMCUpdateBase parameters
   * @param brancher engine which handles branching
   *
   * Update time-step variables to move walkers
   */
  void resetRun(BranchEngineType* brancher, EstimatorManager* est);

  inline RealType getTau()
  {
    //SpeciesSet tspecies(W.getSpeciesSet());
    //int massind=tspecies.addAttribute("mass");
    //RealType mass = tspecies(massind,0);
    //return m_tauovermass*mass;
    return Tau;
  }

  inline void setTau(RealType t)
  {
    //SpeciesSet tspecies(W.getSpeciesSet());
    //int massind=tspecies.addAttribute("mass");
    //RealType mass = tspecies(massind,0);
    //RealType oneovermass = 1.0/mass;
    //RealType oneoversqrtmass = std::sqrt(oneovermass);
//  //     Tau=brancher->getTau();
//  //     assert (Tau==i);
    //m_tauovermass = i/mass;
    Tau=t;
    m_tauovermass = t*MassInvS[0];
    m_oneover2tau = 0.5/(m_tauovermass);
    m_sqrttau = std::sqrt(m_tauovermass);
  }

  inline void getLogs(std::vector<RealType>& logs)
  {
    Psi.getLogs(logs);
  }

  //inline RealType IonKinetic(std::vector<RealType> wfs, RealType h, int edge_x,int edge_y,int edge_z)

  inline RealType IonKinetic(vector<RealType> wfs, RealType h)
  {
    //bool NoneOnEdge;
    RealType m_proton = 1836.15;
    
    return -GridLaplacian(wfs,h)/2.0/m_proton;
    
    
  }

  inline RealType GridLaplacian(vector<RealType> wfs, RealType h)
  {

    RealType x1=wfs[1];
    RealType x3=wfs[2];
    RealType y1=wfs[3];
    RealType y3=wfs[4];
    RealType z1=wfs[5];
    RealType z3=wfs[6];
    RealType f=wfs[7];
    RealType h2=h;//GridDisplacement
    h2 = h2*h2;//GridDisplacement^2    
    return -(x1+x3+y1+y3+z1+z3-6.0*f)/h2;

  }

  inline PosType GridGradient(vector<RealType> wfs, RealType h)
  {

    RealType x1=wfs[1];
    RealType x3=wfs[2];
    RealType y1=wfs[3];
    RealType y3=wfs[4];
    RealType z1=wfs[5];
    RealType z3=wfs[6];

    PosType dr;

    dr(0) = (x3-x1)/h/2.0;
    dr(1) = (y3-y1)/h/2.0;
    dr(2) = (z3-z1)/h/2.0;

    return dr;
  }


  inline RealType IonKineticEnergy(RealType wfs[][2], RealType f,RealType h, RealType ion_mass)
  {
    //bool NoneOnEdge;
    //RealType m_proton = 1836.15;
    
    //return std::max(-GridLaplacianArray(wfs,f,h)/2.0/ion_mass/f,0.0);
    return -GridLaplacianArray(wfs,f,h)/2.0/ion_mass/f;
    
    
  }

  inline RealType IonKineticEnergy3(RealType wfs1[][2], RealType wfs2[][2], RealType f1, RealType f2,RealType h, RealType ion_mass)
  {    
    PosType d1=GridGradientArray(wfs1,h);
    PosType d2=GridGradientArray(wfs2,h);
    RealType f1div = std::max(f1,1.0e-100);
    RealType f2div = std::max(f2,1.0e-100);
    RealType mixTerm=0.0;
    for (int i=0;i<3;++i)
      mixTerm += d1[i]*d2[i];
    mixTerm = mixTerm*2.0/f1div/f2div;

    return -(GridLaplacianArray(wfs1,f1,h)/f1div+GridLaplacianArray(wfs2,f2,h)/f2div+mixTerm)/2.0/ion_mass;    
  }


  inline RealType IonKineticEnergy3_2(RealType wfs1[][2], RealType f1, RealType h, RealType ion_mass, int particle, ParticleSet::ParticlePos_t ionR)
  {    
    PosType d1=GridGradientArray(wfs1,h);
    PosType d2=nuclei_wfs_gradient(ionR,particle);
    RealType f1div = std::max(f1,1.0e-100);
    RealType mixTerm=0.0;
    for (int i=0;i<3;++i)
      mixTerm += d1[i]*d2[i];
    mixTerm = mixTerm*2.0/f1div;

    return -(GridLaplacianArray(wfs1,f1,h)/f1div+nuclei_wfs_laplacian(ionR,particle)+mixTerm)/2.0/ion_mass;    
  }

  /*inline RealType pairGaussian(PosType ion1, PosType ion2, RealType coeff, RealType rc)
  {
    RealType dr = 0.0;
    for (int i=0;i<3;++i)
      dr += (ion1[i]-ion2[i])*(ion1[i]-ion2[i]);

    dr = std::sqrt(dr);

    return std::exp(-(dr-rc)*(dr-rc)*coeff);
  } */

  inline RealType directionalGaussian(PosType ion1, PosType iono, int particle)
  { // directional Gaussian around initial ion position
    //  ion1 and iono better be the current and orginal positions of particle
    RealType dr2 = 0.0;
    for (int coord=0;coord<3;coord++){
      dr2 += ionCoeff[particle][coord] * (ion1[coord]-iono[coord])*(ion1[coord]-iono[coord]);
    }
    return std::exp(-dr2);
  }

  inline RealType nuclei_wfs(ParticleSet::ParticlePos_t ionR)
  {
    RealType n_wfs = 1.0;
    for (int i=0;i<ionR.size();++i){ if (i==ION0) {
	    n_wfs *= directionalGaussian(ionR[i],ionRo[i],i);
    } }

    return n_wfs;
  }
  
  inline PosType nuclei_wfs_gradient(ParticleSet::ParticlePos_t ionR, int particle)
  {
    PosType dr(0.0);
    RealType dist;

    for (int i=0;i<3;i++){
      dr[i] = -2*ionCoeff[particle][i]*(ionR[particle][i]-ionRo[particle][i]);
    }

    return dr;
  }

  inline RealType nuclei_wfs_laplacian(ParticleSet::ParticlePos_t ionR, int particle)
  {
    PosType dr(0.0);
    RealType lap=0.0;
    for (int i=0;i<3;i++){
        RealType coeff = ionCoeff[particle][i];
        lap += -2*coeff+4*coeff*coeff*(ionR[particle][i]-ionRo[particle][i])*(ionR[particle][i]-ionRo[particle][i]);
    }
    return lap;
  }
  
  inline RealType GridLaplacianArray(RealType wfs[][2], RealType f,RealType h, int dim=3)
  {
    
    RealType dr2;
    RealType h2 = h*h;
    
    dr2 = 0.0;
    for (int i=0;i<dim;++i)
      dr2 += (wfs[i][1]+wfs[i][0]-2.0*f)/h2;
    
    return dr2;

  }

  inline RealType IonKineticEnergy2(RealType wfs[][3], RealType f,RealType h, int edge[],RealType ion_mass)
  {
    //bool NoneOnEdge;
    //RealType m_proton = 1836.15;
    
    return -GridLaplacianArray2(wfs,f,h,edge)/2.0/ion_mass/f;
    
    
  }
  

  inline RealType GridLaplacianArray2(RealType wfs[][3], RealType f,RealType h, int edge[],int dim=3)
  {

    RealType dr2;
    RealType h2 = h*h;
    
    dr2 = 0.0;
    for (int i=0;i<dim;++i)
      {
	if (edge[i]==0)
	  dr2 += (wfs[i][1]+wfs[i][0]-2.0*f)/h2;
	else // 2.0*f-5.0*wfs[i][x +/- h]+4.0*wfs[i][x +/- 2h]-wfs[i][x +/- 3h]
	  dr2 += (2.0*f-5.0*wfs[i][0]+4.0*wfs[i][1]-wfs[i][2])/h2;
      }
    
    return dr2;
    
  }


  inline PosType ionDrift(RealType wfs[][2], RealType h, RealType f, RealType Dtau)
  {
    RealType fdiv = std::max(f,1.0e-100);
    return 2.0*GridGradientArray(wfs,h)/fdiv*Dtau;
  }


  inline PosType GridGradientArray(RealType wfs[][2], RealType h, int dim=3)
  {
    PosType dr;

    for (int i=0;i<dim;++i)
      dr(i) = (wfs[i][1]-wfs[i][0])/h/2.0; //wfs[i][x+h]-wfs[i][x-h]

    return dr;
  }

  inline PosType GridGradientArray2(RealType wfs[][3], int edge[], RealType h)
  {

    PosType dr;

    for (int i=0;i<3;++i)
      {	
	if (edge[i]==0) //wfs[i][x+h]-wfs[i][x-h]
	  dr(i) = (wfs[i][1]-wfs[i][0])/h/2.0;
	else if (edge[i]==-1) //-wfs[i][x]+4.0*wfs[i][x+h]-wfs[i][x+2h]
	  dr(i) = (-wfs[i][0]+4.0*wfs[i][1]-wfs[i][2])/2.0/h;
	else // (edge==1) // wfs[i][x]-4.0*wfs[i][x-h]+wfs[i][x-2h]
	  dr(i) = (wfs[i][0]-4.0*wfs[i][1]+wfs[i][2])/2.0/h;
      }

    return dr;
  }

  
  




  //JNKIM: move implementation to QMCUpdateBase
  void tauFrequecies(int steps_in, int ip);

  ///** start a run */
  void startRun(int blocks, bool record);
  /** stop a run */
  void stopRun();
  /** reset the trial energy */
  void resetEtrial(RealType et);

  /** prepare to start a block
   * @param steps number of steps within the block
   */
  void startBlock(int steps);

  /** stop a block
   */
  void stopBlock(bool collectall=true);

  /** set the multiplicity of the walkers to branch */
  void setMultiplicity(WalkerIter_t it, WalkerIter_t it_end);

  /** set the multiplicity of the walkers to branch */
  void setReleasedNodeMultiplicity(WalkerIter_t it, WalkerIter_t it_end);

  /** initialize Walker buffers for PbyP update
   */
  virtual void initWalkersForPbyP(WalkerIter_t it, WalkerIter_t it_end);

  /** initalize Walker for walker update
   */
  virtual void initWalkers(WalkerIter_t it, WalkerIter_t it_end);

  /** update Walker buffers for PbyP update
   */
  void updateWalkers(WalkerIter_t it, WalkerIter_t it_end);

  /** simple routine to test the performance
   */
  void benchMark(WalkerIter_t it, WalkerIter_t it_end, int ip);

  /**  process options
   */
  bool put(xmlNodePtr cur);

  inline void accumulate(WalkerIter_t it, WalkerIter_t it_end)
  {
    Estimators->accumulate(W,it,it_end);
  }

  ///move a walker, all-particle (waler) move, using drift
  void advanceWalker(Walker_t& thisWalker);
  ///move a walker, by particle-by-particle move using fast drift
  void advancePbyP(Walker_t& thisWalker);

  /** advance walkers executed at each step
   *
   * Derived classes implement how to move walkers and accept/reject
   * moves.
   */
  virtual void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)=0;
  virtual RealType advanceWalkerForEE(Walker_t& w1, vector<PosType>& dR, vector<int>& iats, vector<int>& rs, vector<RealType>& ratios)
  {
    return 0.0;
  };
//       virtual RealType advanceWalkerForCSEE(Walker_t& w1, vector<PosType>& dR, vector<int>& iats, vector<int>& rs, vector<RealType>& ratios, vector<RealType>& weights, vector<RealType>& logs ) {return 0.0;};
  virtual void setLogEpsilon(RealType eps) {};
//       virtual void advanceCSWalkers(vector<TrialWaveFunction*>& pclone, vector<MCWalkerConfiguration*>& wclone, vector<QMCHamiltonian*>& hclone, vector<RandomGenerator_t*>& rng, vector<RealType>& c_i){};

  ///normalization offset for cs type runs.
  RealType csoffset;

//       virtual void estimateNormWalkers(vector<TrialWaveFunction*>& pclone
//     , vector<MCWalkerConfiguration*>& wclone
//     , vector<QMCHamiltonian*>& hclone
//     , vector<RandomGenerator_t*>& rng
//     , vector<RealType>& ratio_i_0){};
  int RMC_checkIndex(int N, int NMax)
  {
    if(N<0)
      return N+NMax;
    else
      if (N>=NMax)
        return N-NMax;
      else
        return N;
  }

  void RMC_checkWalkerBounds(WalkerIter_t& it, WalkerIter_t first, WalkerIter_t last)
  {
    if (it>=last)
      it-=(last-first);
    else
      if (it<first)
        it+=(last-first);
  }

  inline RealType logBackwardGF(const ParticleSet::ParticlePos_t& displ)
  {
    RealType t=0.5/Tau;
    RealType logGb=0.0;
    for(int iat=0; iat<W.getTotalNum();++iat)
      logGb += t*MassInvP[iat]*dot(displ[iat],displ[iat]);
    return -logGb;
  }

protected:
  ///update particle-by-particle
  bool UpdatePbyP;
  ///use T-moves
  bool UseTMove;
  ///number of particles
  IndexType NumPtcl;

  ///Time-step factor \f$ 1/(2\Tau)\f$
  RealType m_oneover2tau;
  ///Time-step factor \f$ \sqrt{\Tau}\f$
  RealType m_sqrttau;
  ///tau/mass
  RealType m_tauovermass;
  ///maximum displacement^2
  RealType m_r2max;
  ///walker ensemble
  MCWalkerConfiguration& W;
  ///trial function
  TrialWaveFunction& Psi;
  ///guide function
  TrialWaveFunction& Guide;
  ///Hamiltonian
  QMCHamiltonian& H;
  ///random number generator
  RandomGenerator_t& RandomGen;

  ///branch engine
  BranchEngineType* branchEngine;
  ///estimator
  EstimatorManager* Estimators;
  ///parameters
  ParameterSet myParams;
  ///1/Mass per species
  vector<RealType> MassInvS;
  ///1/Mass per particle
  vector<RealType> MassInvP;
  ///sqrt(tau/Mass) per particle
  vector<RealType> SqrtTauOverMass;
  ///non local operator
  NonLocalTOperator nonLocalOps;
  ///temporary storage for drift
  ParticleSet::ParticlePos_t drift;
  ///temporary storage for random displacement
  ParticleSet::ParticlePos_t deltaR;

  ParticleSet::ParticlePos_t ionPos;
  ///storage for differential gradients for PbyP update
  ParticleSet::ParticleGradient_t G, dG;
  ///storage for differential laplacians for PbyP update
  ParticleSet::ParticleLaplacian_t L, dL;

  /** evaluate the ratio of scaled velocity and velocity
   * @param g gradient
   * @param gscaled scaled gradient
   * @return the ratio
   */
  RealType getNodeCorrection(const ParticleSet::ParticleGradient_t& g, ParticleSet::ParticlePos_t& gscaled);

  ///copy constructor
  QMCUpdateBase(const QMCUpdateBase& a);

  /** a VMC step to randomize awalker
   */
  void randomize(Walker_t& awalker);

private:

  ///set default parameters
  void setDefaults();
  /// Copy operator (disabled).
  QMCUpdateBase& operator=(const QMCUpdateBase&)
  {
    return *this;
  }
};
}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1592 $   $Date: 2007-01-04 16:48:00 -0600 (Thu, 04 Jan 2007) $
 * $Id: QMCUpdateBase.h 1592 2007-01-04 22:48:00Z jnkim $
 ***************************************************************************/
