//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim and Jordan Vincent
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
#include "SQD/HFConfiguration.h"
#include "SQD/SphericalPotential/ZOverRPotential.h"
#include "SQD/SphericalPotential/HarmonicPotential.h"
#include "SQD/SphericalPotential/StepPotential.h"
#include "SQD/SphericalPotential/RegularLinearTransform.h"
#include "SQD/SphericalPotential/RegularLogTransform.h"
#include "SQD/SphericalPotential/NuclearLinearTransform.h"
#include "SQD/SphericalPotential/NuclearLogTransform.h"
#include "SQD/SphericalPotential/NuclearRelLogTransform.h"
#include "Numerics/Numerov.h"
#include "Numerics/RadialFunctorUtility.h"
#include "SQD/HartreeFock.h"

namespace ohmmshf
{

/** Solve a HF eigen problem for a spherically-symmetric external potential.
 * @param norb the number of eigen vectors to be obtained
 */
template<typename Transform_t>
inline void
HartreeFock::run(int norb)
{
  typedef Numerov<Transform_t, RadialOrbital_t> Numerov_t;
  value_type Vtotal,KEnew, KEold,E;
  value_type lowerbound, upperbound;
  vector<value_type> energy(Pot.size());
  eigVal.resize(norb);
  int iter = 0;
  Vtotal = Pot.evaluate(Psi,energy,norb);
  Pot.mix(0.0);
  KEnew = Pot.calcKE(Psi,0,norb);
  string label("spdf");
  std::ofstream log_stream(LogFileName.c_str());
  log_stream.precision(8);
  do
  {
    KEold = KEnew;
    value_type eigsum = 0.0;
    //loop over the orbitals
    for(int ob=0; ob < norb; ob++)
    {
      //set the number of nodes of the eigen vector
      Pot.V[ob].setNumOfNodes(Pot.getNumOfNodes(Psi.N[ob],Psi.L[ob]));
      //calculate the lower and upper bounds for the eigenvalue
      Pot.EnergyBound(lowerbound,upperbound);
      //set up the transformer
      Transform_t es(Pot.V[ob], Psi.N[ob], Psi.L[ob],Psi.CuspParam,
                     Pot.getMass());
      //initialize the numerov solver
      Numerov_t numerov(es,Psi(ob));
      //calculate the eigenvalue and the corresponding orbital
      eigsum += (eigVal[ob] =
                   numerov.solve(lowerbound, upperbound, eig_tol));
      log_stream << Psi.N[ob]<< label[Psi.L[ob]] << '\t'
                 << eigVal[ob] << endl;
    }
    log_stream << endl;
    //normalize the orbitals
    Psi.normalize(norb);
    //restrict the orbitals
    Psi.applyRestriction(norb);
    //calculate the new kinetic energy
    KEnew = Pot.calcKE(Psi,eigsum,norb);
    //the total energy
    //  E = KEnew + Vtotal;
    //for the new orbitals Psi, calculate the new SCF potentials
    Vtotal = Pot.evaluate(Psi,energy,norb);
    //calculate the total energy
    E = KEnew + Vtotal;
    //restrict the potential
    Pot.applyRestriction(Psi);
    //mix the new SCF potential with the old
    Pot.mix(ratio);
    log_stream.precision(10);
    log_stream << "Iteration #" << iter+1 << endl;
    log_stream << "KE    = " << setw(15) << KEnew
               << "  PE     = " << setw(15) << Vtotal << endl;
    log_stream << "PE/KE = " << setw(15) << Vtotal/KEnew
               << "  Energy = " << setw(15) << E << endl;
    log_stream << endl;
    iter++;
    //continue the loop until the kinetic energy converges
  }
  while(fabs(KEnew-KEold)>scf_tol && iter<maxiter);
  log_stream << "V_External = " << energy[0] << endl;
  log_stream << "V_Hartree = "  << energy[1] << endl;
  log_stream << "V_Exchange = " << energy[2] << endl;
  log_stream << "E_tot = " << E << endl;
}

/** Instantiate a Transformation function based on the potential and grid type and call run.
 */
bool HartreeFock::solve()
{
  int norb = Psi.size();
  if(PotType == "nuclear")
  {
    if(GridType == "linear")
    {
      run<NuclearLinearTransform<RadialOrbital_t> >(norb);
    }
    else if(GridType == "log")
    {
      run<NuclearLogTransform<RadialOrbital_t> >(norb);
    }
  }
  else if(PotType == "nuclear_scalar_rel")
  {
    if(GridType == "log")
    {
      run<NuclearRelLogTransform<RadialOrbital_t> >(norb);
    }
  }
  else
  {
    if(GridType == "linear")
    {
      run<RegularLinearTransform<RadialOrbital_t> >(norb);
    }
    else if(GridType == "log")
    {
      run<RegularLogTransform<RadialOrbital_t> >(norb);
    }
  }
  return true;
}


}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 5892 $   $Date: 2013-07-01 14:09:54 -0500 (Mon, 01 Jul 2013) $
 * $Id: HartreeFock.cpp 5892 2013-07-01 19:09:54Z jnkim $
 ***************************************************************************/
