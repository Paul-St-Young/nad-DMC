//////////////////////////////////////////////////////////////////
// (c) Copyright 2007- by Jeongnim Kim
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
/** @file QMCLinearOptimize.h
 * @brief Definition of QMCDriver which performs VMC and optimization.
 */
#ifndef QMCPLUSPLUS_QMCLINEAROPTIMIZATION_VMCSINGLE_H
#define QMCPLUSPLUS_QMCLINEAROPTIMIZATION_VMCSINGLE_H

#include "QMCDrivers/QMCDriver.h"
#include "Optimize/OptimizeBase.h"
#include "QMCApp/WaveFunctionPool.h"
#include "Numerics/LinearFit.h"


namespace qmcplusplus
{

///forward declaration of a cost function
class QMCCostFunctionBase;
class HamiltonianPool;

/** @ingroup QMCDrivers
 * @brief Implements wave-function optimization
 *
 * Optimization by correlated sampling method with configurations
 * generated from VMC.
 */

class QMCLinearOptimize: public QMCDriver
{
public:

  ///Constructor.
  QMCLinearOptimize(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                    QMCHamiltonian& h, HamiltonianPool& hpool, WaveFunctionPool& ppool);

  ///Destructor
  virtual ~QMCLinearOptimize();

  ///Run the Optimization algorithm.
  virtual bool run()=0;
  ///process xml node
  virtual bool put(xmlNodePtr cur);
  void resetComponents(xmlNodePtr cur);
  ///add a configuration file to the list of files
  void addConfiguration(const string& a);
  void setWaveFunctionNode(xmlNodePtr cur)
  {
    wfNode=cur;
  }

  vector<RealType> optdir, optparm;
  ///index to denote the partition id
  int PartID;
  ///total number of partitions that will share a set of configuratons
  int NumParts;
  ///total number of Warmup Blocks
  int WarmupBlocks;
  ///total number of Warmup Blocks
  int NumOfVMCWalkers;
  ///Number of iterations maximum before generating new configurations.
  int Max_iterations;
  ///need to know HamiltonianPool to use OMP
  HamiltonianPool& hamPool;
  ///target cost function to optimize
  QMCCostFunctionBase* optTarget;
  ///Dimension of matrix and number of parameters
  int N,numParams;
  ///vmc engine
  QMCDriver* vmcEngine;
  ///xml node to be dumped
  xmlNodePtr wfNode;
  ///xml node for optimizer
  xmlNodePtr optNode;
  ///list of files storing configurations
  vector<string> ConfigFile;

  RealType param_tol;

  inline bool tooLow(RealType safeValue, RealType CurrentValue)
  {
    RealType lowestCostAllowed=std::min(safeValue-0.03*std::abs(safeValue),safeValue-3.0);
    if (CurrentValue<lowestCostAllowed)
      return true;
    else
      return false;
  }

  bool fitMappedStabilizers(vector<std::pair<RealType,RealType> >& mappedStabilizers, RealType& XS, RealType& val, RealType tooBig );

  inline int CubicFormula (double a, double b, double c, double d,
                           double &x1, double &x2, double &x3)
  {
    double A = b/a;
    double B = c/a;
    double C = d/a;
    double Q = (A*A - 3.0*B)/9.0;
    double R = (2.0*A*A*A - 9.0*A*B + 27.0*C)/54.0;
    //cerr << "Q = " << Q << " R = " << R << "\n";
    if ((R*R) < (Q*Q*Q))
    {
      double theta = std::acos(R/std::sqrt(Q*Q*Q));
      double twosqrtQ = 2.0*std::sqrt(Q);
      double third = 1.0/3.0;
      double thirdA = third * A;
      x1 = -twosqrtQ*std::cos(third*theta) - thirdA;
      x2 = -twosqrtQ*std::cos(third*(theta + 2.0*M_PI)) - thirdA;
      x3 = -twosqrtQ*std::cos(third*(theta - 2.0*M_PI)) - thirdA;
      return 3;
    }
    else
    {
      double D = -Q*Q*Q + R*R;
      double u = cbrt(-R + std::sqrt(D));
      double v = cbrt(-R - std::sqrt(D));
      double y1 = u+v;
      x1 = y1 - A/3.0;
      return 1;
    }
  }

  inline RealType QuarticMinimum (vector<RealType> &coefs)
  {
    double a, b, c, d;
    a = 4.0*coefs[4];
    b = 3.0*coefs[3];
    c = 2.0*coefs[2];
    d = coefs[1];
    double x1, x2, x3;
    int numroots = CubicFormula (a, b, c, d, x1, x2, x3);
    if (numroots == 1)
      return x1;
    else
    {
      double v1 = coefs[0] + coefs[1]*x1 + coefs[2]*x1*x1 + coefs[3]*x1*x1*x1
                  + coefs[4]*x1*x1*x1*x1;
      double v2 = coefs[0] + coefs[1]*x2 + coefs[2]*x2*x2 + coefs[3]*x2*x2*x2
                  + coefs[4]*x2*x2*x2*x2;
      double v3 = coefs[0] + coefs[1]*x3 + coefs[2]*x3*x3 + coefs[3]*x3*x3*x3
                  + coefs[4]*x3*x3*x3*x3;
      if (v1 < v2 && v1 < v3)
        return x1;
      if (v2 < v1 && v2 < v3)
        return x2;
      if (v3 < v1 && v3 < v2)
        return x3;
      return x1;
    }
  }

  ///common operation to start optimization, used by the derived classes
  void start();
  ///common operation to finish optimization, used by the derived classes
  void finish();
  //asymmetric generalized EV
  RealType getLowestEigenvector(Matrix<RealType>& A, Matrix<RealType>& B, vector<RealType>& ev);
  //asymmetric EV
  RealType getLowestEigenvector(Matrix<RealType>& A, vector<RealType>& ev);
  RealType getSplitEigenvectors(int first, int last, Matrix<RealType>& FullLeft, Matrix<RealType>& FullRight, vector<RealType>& FullEV, vector<RealType>& LocalEV, string CSF_Option, bool& CSF_scaled);
  void getNonLinearRange(int& first, int& last);
  void orthoScale(std::vector<RealType>& dP, Matrix<RealType>& S);
  bool nonLinearRescale( vector<RealType>& dP, Matrix<RealType>& S);
  RealType getNonLinearRescale( vector<RealType>& dP, Matrix<RealType>& S);
  void generateSamples();
  void add_timers(vector<NewTimer*>& timers);
  vector<NewTimer*> myTimers;
  Timer t1;
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 5900 $   $Date: 2013-07-10 13:51:52 -0500 (Wed, 10 Jul 2013) $
 * $Id: QMCLinearOptimize.h 5900 2013-07-10 18:51:52Z jnkim $
 ***************************************************************************/
