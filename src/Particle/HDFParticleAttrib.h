//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
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
#ifndef OHMMS_HDF_PARTICLEATTRIBIO_H
#define OHMMS_HDF_PARTICLEATTRIBIO_H
#include "OhmmsData/HDFAttribIO.h"

namespace qmcplusplus
{
/** Specialization for ParticleAttrib<int> */
template<>
struct HDFAttribIO<ParticleAttrib<int> >: public HDFAttribIOBase
{

  typedef ParticleAttrib<int> ArrayType_t;
  ArrayType_t&  ref;

  HDFAttribIO<ArrayType_t>(ArrayType_t& a):ref(a) { }

  void write(hid_t grp, const char* name)
  {
    hsize_t dim = ref.size();
    hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
    hid_t dataset =
      H5Dcreate(grp, name, H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
    hid_t ret =
      H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&ref[0]);
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }

  void read(hid_t grp, const char* name)
  {
    hid_t h1 = H5Dopen(grp, name);
    hid_t ret = H5Dread(h1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ref[0]);
    H5Dclose(h1);
  }
};

/** Specialization for ParticleAttrib<double> */
template<>
struct HDFAttribIO<ParticleAttrib<double> >: public HDFAttribIOBase
{

  typedef ParticleAttrib<double> ArrayType_t;
  ArrayType_t&  ref;

  HDFAttribIO<ArrayType_t>(ArrayType_t& a):ref(a) { }

  void write(hid_t grp, const char* name)
  {
    hsize_t dim = ref.size();
    hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
    hid_t dataset =
      H5Dcreate(grp, name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
    hid_t ret =
      H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&ref[0]);
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }

  void read(hid_t grp, const char* name)
  {
    hid_t h1 = H5Dopen(grp, name);
    hid_t ret = H5Dread(h1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ref[0]);
    H5Dclose(h1);
  }

};

/** Specialization for ParticleAttrib<TinyVector<double,3> > > */
template<>
struct HDFAttribIO<ParticleAttrib<TinyVector<double,3> > >: public HDFAttribIOBase
{

  typedef TinyVector<double,3> SingleParticlePos_t;
  typedef ParticleAttrib<SingleParticlePos_t> ArrayType_t;

  ArrayType_t&  ref;

  HDFAttribIO<ArrayType_t>(ArrayType_t& a):ref(a) { }

  ~HDFAttribIO<ParticleAttrib<TinyVector<double,3> > > () { }

  void write(hid_t grp, const char* name)
  {
    hsize_t dims[2];
    dims[0] = ref.size();
    dims[1] = 3;
    hid_t dataspace  = H5Screate_simple(2, dims, NULL);
    hid_t dataset =
      H5Dcreate(grp, name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
    hid_t ret =
      H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&(ref[0][0]));
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }

  void read(hid_t grp, const char* name)
  {
    hid_t h1 = H5Dopen(grp, name);
    hid_t ret = H5Dread(h1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                        &(ref[0][0]));
    H5Dclose(h1);
  }
};

}
// // specialization for 2D vector type, just for fun
// template<>
// struct HDFParticleAttrib<TinyVector<double,2> > : public HDFParticleAttribBase {

//   typedef TinyVector<double,2> SingleParticlePos_t;
//   typedef ParticleAttrib<SingleParticlePos_t> ArrayType_t;

//   ArrayType_t&  ref;

//   // possibly using trait class make sense but my guess is that it won't affect the memory
//   // nor performance
//   hid_t PosID;

//   HDFParticleAttrib<TinyVector<double,2> > (ArrayType_t& a):ref(a) {
//     PosID = H5Tcreate(H5T_COMPOUND, 2*sizeof(double));
//     H5Tinsert(PosID,"x",offsetof(SingleParticlePos_t, X[0]), H5T_NATIVE_DOUBLE);
//     H5Tinsert(PosID,"y",offsetof(SingleParticlePos_t, X[1]), H5T_NATIVE_DOUBLE);
//   }

//   void write(hid_t  grp, const char* name) {
//     hsize_t dim = ref.size();
//     hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
//     hid_t dataset =
//       H5Dcreate(grp, name, PosID, dataspace, H5P_DEFAULT);
//     hid_t ret =
//       H5Dwrite(dataset, PosID, H5S_ALL, H5S_ALL, H5P_DEFAULT, ref.begin());
//     H5Sclose(dataspace);
//     H5Dclose(dataset);
//   }

//   void read(hid_t  grp, const char* name) {
//     hid_t h1 = H5Dopen(grp, name);
//     hid_t ret = H5Dread(h1, PosID, H5S_ALL, H5S_ALL, H5P_DEFAULT, ref.begin());
//     H5Dclose(h1);
//   }
// };
#endif
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 19:14:53 -0500 (Thu, 25 Apr 2013) $
 * $Id: HDFParticleAttrib.h 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/
