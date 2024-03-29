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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-

/**@file ParticleTags.cpp
 *@brief initialization of the static data members of ParticleTags
 */
#include "ParticleTags.h"

//Carry over from old implementation
// #ifdef OHMMS_USE_OLD_PARTICLETYPE
// #define POSSTR       "PositionType"
// #define INDEXSTR     "IndexType"
// #define SCALARSTR    "ScalarType"
// #define TENSORSTR    "TensorType"
// #define POSITION_TAG "positions"
// #define ID_TAG       "id"
// #define IONID_TAG    "ionid"
// #else
// #define POSSTR       "posArray"
// #define INDEXSTR     "indexArray"
// #define SCALARSTR    "scalarArray"
// #define TENSORSTR    "tensorArray"
// #define STRINGSTR    "stringArray"
// #define POSITION_TAG "position"
// #define ID_TAG       "id"
// #define IONID_TAG    "ionid"
// #endif

std::string ParticleTags::null_tag="null";

std::string ParticleTags::indextype_tag="indexArray";
std::string ParticleTags::scalartype_tag="scalarArray";
std::string ParticleTags::stringtype_tag="stringArray";
std::string ParticleTags::postype_tag="posArray";
std::string ParticleTags::gradtype_tag="gradArray";
std::string ParticleTags::laptype_tag="lapArray";
std::string ParticleTags::tensortype_tag="tensorArray";
std::string ParticleTags::xmoltype_tag="xmolArray";

std::string ParticleTags::position_tag="position";
std::string ParticleTags::id_tag="id";
std::string ParticleTags::ionid_tag="ionid";
std::string ParticleTags::trajectory_tag="trajectory";
std::string ParticleTags::force_tag="f";
std::string ParticleTags::velocity_tag="v";
std::string ParticleTags::energy_tag="e";
std::string ParticleTags::sumbc_tag="sumbc";

std::string ParticleTags::root_tag="particleset";
std::string ParticleTags::attrib_tag="attrib";
std::string ParticleTags::name_tag="name";
std::string ParticleTags::datatype_tag="datatype";
std::string ParticleTags::condition_tag="condition";
std::string ParticleTags::size_tag="size";
std::string ParticleTags::format_tag="format";
std::string ParticleTags::role_tag="role";

/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 19:14:53 -0500 (Thu, 25 Apr 2013) $
 * $Id: ParticleTags.cpp 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/
