//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
/** @file BasisSetBase.h
 * @brief Declaration of a base class of BasisSet
 */
#ifndef QMCPLUSPLUS_ORBITALSETTRAITS_H
#define QMCPLUSPLUS_ORBITALSETTRAITS_H

#include "Configuration.h"
#include "type_traits/scalar_traits.h"
#include "Optimize/VariableSet.h"

namespace qmcplusplus
{
/** dummy class for templated classes
 */
struct DummyGrid
{
  inline void locate(double r) {}
  DummyGrid* makeClone() const
  {
    return new DummyGrid;
  }
};

typedef TinyVector<int,4> QuantumNumberType;

enum {q_n=0,q_l,q_m, q_s};

/** trait class to handel a set of Orbitals
 */
template<typename T>
struct OrbitalSetTraits//: public OrbitalTraits<T>
{
  enum {DIM=OHMMS_DIM};
  typedef typename scalar_traits <T>::real_type RealType;
  typedef typename scalar_traits <T>::value_type ValueType;
  typedef int                            IndexType;
  typedef TinyVector<RealType,DIM>       PosType;
  typedef TinyVector<ValueType,DIM>      GradType;
  typedef Tensor<ValueType,DIM>          HessType;
  typedef Tensor<ValueType,DIM>          TensorType;
  typedef TinyVector<Tensor<ValueType,DIM>,DIM> GradHessType;
  typedef Vector<IndexType>     IndexVector_t;
  typedef Vector<ValueType>     ValueVector_t;
  typedef Matrix<ValueType>     ValueMatrix_t;
  typedef Vector<GradType>      GradVector_t;
  typedef Matrix<GradType>      GradMatrix_t;
  typedef Vector<HessType>      HessVector_t;
  typedef Matrix<HessType>      HessMatrix_t;
  typedef Vector<GradHessType>  GradHessVector_t;
  typedef Matrix<GradHessType>  GradHessMatrix_t;

};

///typedef for a set of variables that are varied during an optimization
typedef optimize::VariableSet  opt_variables_type;
///typedef for a set of variables that can be varied
typedef optimize::VariableSet::variable_map_type variable_map_type;


template<typename T> inline double evaluatePhase(T sign_v)
{
  return (sign_v>0)?0.0:M_PI;
}

template<>
inline double evaluatePhase(const std::complex<double>& psi)
{
  return std::arg(psi);
}

/** evaluate the log(|psi|) and phase
 * @param psi real/complex value
 * @param phase phase of psi
 * @return log(|psi|)
 */
template<class T>
inline T evaluateLogAndPhase(const T psi, T& phase)
{
  if(psi<0.0)
  {
    phase= M_PI;
    return std::log(-psi);
  }
  else
  {
    phase = 0.0;
    return std::log(psi);
  }
}

template<class T>
inline T
evaluateLogAndPhase(const std::complex<T>& psi, T& phase)
{
  phase = std::arg(psi);
  if(phase<0.0)
    phase += 2.0*M_PI;
  return std::log( std::abs(psi) );
//      return 0.5*std::log(psi.real()*psi.real()+psi.imag()*psi.imag());
  //return std::log(psi);
}

inline double evaluatePhase(const double psi)
{
  return (psi<numeric_limits<double>::epsilon())?M_PI:0.0;
}

inline double evaluatePhase(const std::complex<double>& psi)
{
  return std::arg(psi);
}

}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 19:14:53 -0500 (Thu, 25 Apr 2013) $
 * $Id: OrbitalSetTraits.h 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/
