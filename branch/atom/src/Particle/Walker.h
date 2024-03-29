//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_WALKER_H
#define QMCPLUSPLUS_WALKER_H

#include "OhmmsPETE/OhmmsMatrix.h"
#include "Utilities/PooledData.h"
#ifdef QMC_CUDA
#include "Utilities/PointerPool.h"
#include "CUDA/gpu_vector.h"
#endif
#include <assert.h>
#include <deque>
namespace qmcplusplus
{

/** an enum denoting index of physical properties
 *
 * LOCALPOTENTIAL should be always the last enumeation
 * When a new enum is needed, modify ParticleSet::initPropertyList to match the list
 */
enum {LOGPSI=0,       /*!< log(fabs(psi)) instead of square of the many-body wavefunction \f$|\Psi|^2\f$ */
      SIGN,           /*!< value of the many-body wavefunction \f$\Psi(\{R\})\f$ */
      UMBRELLAWEIGHT, /*!< sum of wavefunction ratios for multiple H and Psi */
      R2ACCEPTED,     /*!< r^2 for accepted moves */
      R2PROPOSED,     /*!< r^2 for proposed moves */
      DRIFTSCALE,     /*!< scaling value for the drift */
      ALTERNATEENERGY,  /*!< alternatelocal energy, the sum of all the components */
      LOCALENERGY,    /*!< local energy, the sum of all the components */
      LOCALPOTENTIAL, /*!< local potential energy = local energy - kinetic energy */
      NUMPROPERTIES   /*!< the number of properties */
     };

/** A container class to represent a walker.
 *
 * A walker stores the particle configurations {R}  and a property container.
 * RealTypehe template (P)articleSet(A)ttribute is a generic container  of position types.
 * RealTypehe template (G)radient(A)ttribute is a generic container of gradients types.
 * Data members for each walker
 * - ID : identity for a walker. default is 0.
 * - Age : generation after a move is accepted.
 * - Weight : weight to take the ensemble averages
 * - Multiplicity : multiplicity for branching. Probably can be removed.
 * - Properties  : 2D container. RealTypehe first index corresponds to the H/Psi index and second index >=NUMPROPERTIES.
 * - DataSet : anonymous container.
 */
template<typename t_traits, typename p_traits>
struct Walker
{
  enum {DIM=t_traits::DIM};
  /** typedef for real data type */
  typedef typename t_traits::RealType RealType;
  /** typedef for value data type. */
  typedef typename t_traits::ValueType ValueType;
  /** array of particles */
  typedef typename p_traits::ParticlePos_t ParticlePos_t;
  /** array of gradients */
  typedef typename p_traits::ParticleGradient_t ParticleGradient_t;
  /** array of laplacians */
  typedef typename p_traits::ParticleLaplacian_t ParticleLaplacian_t;

  ///typedef for the property container, fixed size
  typedef Matrix<RealType>      PropertyContainer_t;
  typedef PooledData<RealType>  Buffer_t;

  ///id reserved for forward walking
  long ID;
  ///id reserved for forward walking
  long ParentID;
  ///DMCgeneration
  int Generation;
  ///Age of this walker age is incremented when a walker is not moved after a sweep
  int Age;
  ///Age of this walker age is incremented when a walker is not moved after a sweep
  int ReleasedNodeAge;
  ///Weight of the walker
  RealType Weight;
  ///Weight of the walker
  RealType ReleasedNodeWeight;
  /** Number of copies for branching
   *
   * When Multiplicity = 0, this walker will be destroyed.
   */
  RealType Multiplicity;

  /**the configuration vector (3N-dimensional vector to store
     the positions of all the particles for a single walker)*/
  ParticlePos_t R;

  ParticlePos_t ionPos;

  /** \f$ \nabla_i d\log \Psi for the i-th particle */
  ParticleGradient_t G;
  /** \f$ \nabla^2_i d\log \Psi for the i-th particle */
  ParticleLaplacian_t L;
  /////drift of the walker \f$ Drift({\bf R}) = \tau v_{drift}({\bf R}) \f$
  //ParticlePos_t Drift;

  ///scalar properties of a walker
  PropertyContainer_t  Properties;

  ///Property history vector
  vector<vector<RealType> >  PropertyHistory;
  vector<int> PHindex;

  ///buffer for the data for particle-by-particle update
  Buffer_t DataSet;

  ///buffer for the constant data in the evaluation of
  //analytical derivatives during linear optimization, e.g. MultiDeterminants
  Buffer_t DataSetForDerivatives;

  /// Data for GPU-vectorized versions
#ifdef QMC_CUDA
  static int cuda_DataSize;
  typedef gpu::device_vector<CUDA_PRECISION> cuda_Buffer_t;
  cuda_Buffer_t cuda_DataSet;
  // Note that R_GPU has size N+1.  The last element contains the
  // proposed position for single-particle moves.
  gpu::device_vector<TinyVector<CUDA_PRECISION,OHMMS_DIM> > R_GPU, Grad_GPU;
  gpu::device_vector<CUDA_PRECISION> Lap_GPU, Rhok_GPU;
  int k_species_stride;
  inline void resizeCuda(int size, int num_species, int num_k)
  {
    cuda_DataSize = size;
    cuda_DataSet.resize(size);
    int N = R.size();
    R_GPU.resize(N);
    Grad_GPU.resize(N);
    Lap_GPU.resize(N);
    // For GPU coallescing
    k_species_stride = ((2*num_k + 15)/16) * 16;
    if (num_k)
      Rhok_GPU.resize (num_species * k_species_stride);
  }
  inline CUDA_PRECISION* get_rhok_ptr ()
  {
    return Rhok_GPU.data();
  }
  inline CUDA_PRECISION* get_rhok_ptr (int isp)
  {
    return Rhok_GPU.data() + k_species_stride * isp;
  }

#endif

  ///create a walker for n-particles
  inline explicit Walker(int nptcl=0, int nions=0)
#ifdef QMC_CUDA
    :cuda_DataSet("Walker::walker_buffer"), R_GPU("Walker::R_GPU"),
     Grad_GPU("Walker::Grad_GPU"), Lap_GPU("Walker::Lap_GPU"),
     Rhok_GPU("Walker::Rhok_GPU")
#endif
  {
    ID=0;
    ParentID=0;
    Generation=0;
    Age=0;
    Weight=1.0;
    Multiplicity=1.0;
    ReleasedNodeWeight=1.0;
    ReleasedNodeAge=0;
    Properties.resize(1,NUMPROPERTIES);
    if(nptcl>0)
      resize(nptcl);
    if(nions>0)
      ions_resize(nions);
    Properties=0.0;
  }


  inline int addPropertyHistory(int leng)
  {
    int newL = PropertyHistory.size();
    vector<RealType> newVecHistory=vector<RealType>(leng,0.0);
    PropertyHistory.push_back(newVecHistory);
    PHindex.push_back(0);
    return newL;
  }

  inline void deletePropertyHistory()
  {
    PropertyHistory.erase(PropertyHistory.begin(), PropertyHistory.end());
  }

  inline void resetPropertyHistory()
  {
    for (int i=0; i<PropertyHistory.size(); i++)
    {
      PHindex[i]=0;
      for (int k=0; k<PropertyHistory[i].size(); k++)
      {
        PropertyHistory[i][k]=0.0;
      }
    }
  }

  inline void addPropertyHistoryPoint(int index, RealType data)
  {
    PropertyHistory[index][PHindex[index]]=(data);
    PHindex[index]++;
    if (PHindex[index]==PropertyHistory[index].size())
      PHindex[index]=0;
//       PropertyHistory[index].pop_back();
  }

  inline RealType getPropertyHistorySum(int index, int endN)
  {
    RealType mean=0.0;
    typename vector<RealType>::const_iterator phStart;
    phStart=PropertyHistory[index].begin()+PHindex[index];
    for (int i=0; i<endN; phStart++,i++)
    {
      if (phStart>=PropertyHistory[index].end())
        phStart -= PropertyHistory[index].size();
      mean+= (*phStart);
    }
    return mean ;
  }

  inline ~Walker() { }

  ///assignment operator
  inline Walker& operator=(const Walker& a)
  {
    if (this != &a)
      makeCopy(a);
    return *this;
  }

  ///return the number of particles per walker
  inline int size() const
  {
    return R.size();
  }

  ///resize for n particles
  inline void resize(int nptcl)
  {
    R.resize(nptcl);
    G.resize(nptcl);
    L.resize(nptcl);
#ifdef QMC_CUDA
    R_GPU.resize(nptcl);
    Grad_GPU.resize(nptcl);
    Lap_GPU.resize(nptcl);
#endif
    //Drift.resize(nptcl);
  }
  
  inline void ions_resize(int nions)
  {
    ionPos.resize(nions);
  }

  ///copy the content of a walker
  inline void makeCopy(const Walker& a)
  {
    ID=a.ID;
    ParentID=a.ParentID;
    Generation=a.Generation;
    Age=a.Age;
    Weight=a.Weight;
    Multiplicity=a.Multiplicity;
    ReleasedNodeWeight=a.ReleasedNodeWeight;
    ReleasedNodeAge=a.ReleasedNodeAge;
    if (R.size()!=a.R.size())
      resize(a.R.size());
    if (ionPos.size()!=a.ionPos.size())
      ions_resize(a.ionPos.size());
    ionPos = a.ionPos;
    R = a.R;
    G = a.G;
    L = a.L;
    //Drift = a.Drift;
    Properties.copy(a.Properties);
    DataSet=a.DataSet;
    if (PropertyHistory.size()!=a.PropertyHistory.size())
      PropertyHistory.resize(a.PropertyHistory.size());
    for (int i=0; i<PropertyHistory.size(); i++)
      PropertyHistory[i]=a.PropertyHistory[i];
    PHindex=a.PHindex;
#ifdef QMC_CUDA
    cuda_DataSet = a.cuda_DataSet;
    R_GPU = a.R_GPU;
    Grad_GPU = a.Grad_GPU;
    Lap_GPU = a.Lap_GPU;
#endif
  }

  //return the address of the values of Hamiltonian terms
  inline RealType* restrict getPropertyBase()
  {
    return Properties.data();
  }

  //return the address of the values of Hamiltonian terms
  inline const RealType* restrict getPropertyBase() const
  {
    return Properties.data();
  }

  ///return the address of the i-th properties
  inline RealType* restrict getPropertyBase(int i)
  {
    return Properties[i];
  }

  ///return the address of the i-th properties
  inline const RealType* restrict getPropertyBase(int i) const
  {
    return Properties[i];
  }


  /** reset the property of a walker
   *@param logpsi \f$\log |\Psi|\f$
   *@param sigN  sign of the trial wavefunction
   *@param ene the local energy
   *
   *Assign the values and reset the age
   * but leave the weight and multiplicity
   */
  inline void resetProperty(RealType logpsi, RealType sigN, RealType ene)
  {
    Age=0;
    //Weight=1.0;
    Properties(LOGPSI)=logpsi;
    Properties(SIGN)=sigN;
    Properties(LOCALENERGY) = ene;
  }

  inline void resetReleasedNodeProperty(RealType localenergy, RealType alternateEnergy, RealType altR)
  {
    Properties(ALTERNATEENERGY)=alternateEnergy;
    Properties(LOCALENERGY) = localenergy;
    Properties(SIGN) = altR;
  }
  inline void resetReleasedNodeProperty(RealType localenergy, RealType alternateEnergy)
  {
    Properties(ALTERNATEENERGY)=alternateEnergy;
    Properties(LOCALENERGY) = localenergy;
  }
  /** reset the property of a walker
   * @param logpsi \f$\log |\Psi|\f$
   * @param sigN  sign of the trial wavefunction
   * @param ene the local energy
   * @param r2a \f$r^2\f$ for the accepted moves
   * @param r2p \f$r^2\f$ for the proposed moves
   * @param vq \f$\bar{V}/V\f$ scaling to control node divergency in JCP 93
   *
   *Assign the values and reset the age
   * but leave the weight and multiplicity
   */
  inline void resetProperty(RealType logpsi, RealType sigN, RealType ene, RealType r2a, RealType r2p, RealType vq)
  {
    Age=0;
    Properties(LOGPSI)=logpsi;
    Properties(SIGN)=sigN;
    Properties(LOCALENERGY) = ene;
    Properties(R2ACCEPTED) = r2a;
    Properties(R2PROPOSED) = r2p;
    Properties(DRIFTSCALE) = vq;
  }

  /** marked to die
       *
       * Multiplicity and weight are set to zero.
       */
  inline void willDie()
  {
    Multiplicity=0;
    Weight=0.0;
  }

  /** reset the walker weight, multiplicity and age */
  inline void reset()
  {
    Age=0;
    Multiplicity=1.0e0;
    Weight=1.0e0;
  }

  inline void resizeProperty(int n, int m)
  {
    Properties.resize(n,m);
  }


  /** byte size for a packed message
   *
   * ID, Age, Properties, R, Drift, DataSet is packed
   */
  inline int byteSize()
  {
    int numPH(0);
    for (int iat=0; iat<PropertyHistory.size(); iat++)
      numPH += PropertyHistory[iat].size();
    int bsize =
      2*sizeof(long)+3*sizeof(int)+ PHindex.size()*sizeof(int)
      +(Properties.size()+DataSet.size()+ numPH + 1)*sizeof(RealType)
      +R.size()*(DIM*sizeof(RealType)+(DIM+1)*sizeof(ValueType))//;//R+G+L
      +ionPos.size()*(DIM*sizeof(RealType));
    //+R.size()*(DIM*2*sizeof(RealType)+(DIM+1)*sizeof(ValueType));//R+Drift+G+L
#ifdef QMC_CUDA
    bsize += 3 *sizeof (int); // size and N and M
    bsize += cuda_DataSize               * sizeof(CUDA_PRECISION); // cuda_DataSet
    bsize += R.size()        * OHMMS_DIM * sizeof(CUDA_PRECISION); // R_GPU
    bsize += R.size()        * OHMMS_DIM * sizeof(CUDA_PRECISION); // Grad_GPU
    bsize += R.size()        * 1         * sizeof(CUDA_PRECISION); // Lap_GPU
    bsize += Rhok_GPU.size()             * sizeof(CUDA_PRECISION); // Lap_GPU
#endif
    return bsize;
  }

  template<class Msg>
  inline Msg& putMessage(Msg& m)
  {
    const int nat=R.size();
    const int nions=ionPos.size();
    m << ID << ParentID << Generation << Age << ReleasedNodeAge << ReleasedNodeWeight;
    m.Pack(&(R[0][0]),nat*OHMMS_DIM);
    m.Pack(&(ionPos[0][0]),nions*OHMMS_DIM);
#if defined(QMC_COMPLEX)
    m.Pack(reinterpret_cast<RealType*>(&(G[0][0])),nat*OHMMS_DIM*2);
    m.Pack(reinterpret_cast<RealType*>(L.first_address()),nat*2);
#else
    m.Pack(&(G[0][0]),nat*OHMMS_DIM);
    m.Pack(L.first_address(),nat);
#endif
    m.Pack(Properties.data(),Properties.size());
    m.Pack(DataSet.data(),DataSet.size());
    //Properties.putMessage(m);
    //DataSet.putMessage(m);
    for (int iat=0; iat<PropertyHistory.size(); iat++)
      m.Pack(&(PropertyHistory[iat][0]),PropertyHistory[iat].size());
    m.Pack(&(PHindex[0]),PHindex.size());
#ifdef QMC_CUDA
    // Pack GPU data
    std::vector<CUDA_PRECISION> host_data, host_rhok;
    std::vector<TinyVector<CUDA_PRECISION,OHMMS_DIM> > R_host;
    std::vector<CUDA_PRECISION> host_lapl;
    cuda_DataSet.copyFromGPU(host_data);
    R_GPU.copyFromGPU(R_host);
    int size = host_data.size();
    int N = R_host.size();
    m.Pack(size);
    m.Pack(N);
    m.Pack(&(host_data[0]), host_data.size());
    m.Pack(&(R_host[0][0]), OHMMS_DIM*R_host.size());
    Grad_GPU.copyFromGPU(R_host);
    m.Pack(&(R_host[0][0]), OHMMS_DIM*R_host.size());
    Lap_GPU.copyFromGPU(host_lapl);
    m.Pack(&(host_lapl[0]), host_lapl.size());
    Rhok_GPU.copyFromGPU(host_rhok);
    int M = host_rhok.size();
    m.Pack(M);
    m.Pack(&(host_rhok[0]), host_rhok.size());
#endif
    return m;
  }

  template<class Msg>
  inline Msg& getMessage(Msg& m)
  {
    const int nat=R.size();
    const int nions=ionPos.size();
    m>>ID >> ParentID >> Generation >> Age >> ReleasedNodeAge >> ReleasedNodeWeight;
    m.Unpack(&(R[0][0]),nat*OHMMS_DIM);
    m.Unpack(&(ionPos[0][0]),nions*OHMMS_DIM);
#if defined(QMC_COMPLEX)
    m.Unpack(reinterpret_cast<RealType*>(&(G[0][0])),nat*OHMMS_DIM*2);
    m.Unpack(reinterpret_cast<RealType*>(L.first_address()),nat*2);
#else
    m.Unpack(&(G[0][0]),nat*OHMMS_DIM);
    m.Unpack(L.first_address(),nat);
#endif
    m.Unpack(Properties.data(),Properties.size());
    m.Unpack(DataSet.data(),DataSet.size());
    //Properties.getMessage(m);
    //DataSet.getMessage(m);
    for (int iat=0; iat<PropertyHistory.size(); iat++)
      m.Unpack(&(PropertyHistory[iat][0]),PropertyHistory[iat].size());
    m.Unpack(&(PHindex[0]),PHindex.size());
#ifdef QMC_CUDA
    // Pack GPU data
    std::vector<CUDA_PRECISION> host_data, host_rhok;
    std::vector<TinyVector<CUDA_PRECISION,OHMMS_DIM> > R_host;
    std::vector<CUDA_PRECISION> host_lapl;
    int size, N;
    m.Unpack(size);
    m.Unpack(N);
    host_data.resize(size);
    R_host.resize(N);
    host_lapl.resize(N);
    m.Unpack(&(host_data[0]), size);
    cuda_DataSet = host_data;
    m.Unpack(&(R_host[0][0]), OHMMS_DIM*N);
    R_GPU = R_host;
    m.Unpack(&(R_host[0][0]), OHMMS_DIM*N);
    Grad_GPU = R_host;
    m.Unpack(&(host_lapl[0]), N);
    Lap_GPU = host_lapl;
    int M;
    m.Unpack(M);
    host_rhok.resize(M);
    m.Unpack(&(host_rhok[0]), M);
    Rhok_GPU = host_rhok;
#endif
    return m;
  }

};

template<class RealType, class PA>
ostream& operator<<(ostream& out, const Walker<RealType,PA>& rhs)
{
  copy(rhs.Properties.begin(), rhs.Properties.end(),
       ostream_iterator<double>(out," "));
  out << endl;
  out << rhs.R;
  return out;
}
}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 5820 $   $Date: 2013-05-03 14:58:44 -0500 (Fri, 03 May 2013) $
 * $Id: Walker.h 5820 2013-05-03 19:58:44Z jnkim $
 ***************************************************************************/
