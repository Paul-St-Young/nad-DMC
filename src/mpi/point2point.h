#ifndef QMCPLUSPLUS_BOOSTADAPTOR_POINT_OPERATIONS_H
#define QMCPLUSPLUS_BOOSTADAPTOR_POINT_OPERATIONS_H

#include <type_traits/container_proxy.h>
#include <mpi/mpi_datatype.h>
//#if defined(HAVE_MPI)
//#include <boost/mpi/operations.hpp>
//#endif

namespace qmcplusplus
{

namespace mpi
{

#if defined(HAVE_MPI)

template<typename T>
inline void send(const communicator& comm, T& in,  int dest, int tag)
{
  container_proxy<T> tmp(in);
  MPI_Datatype type_id=get_mpi_datatype(*tmp.data());
  MPI_Send(tmp.data(),tmp.size(),type_id,dest,tag,comm);
}

template<typename T>
inline void recv(const communicator& comm, T& out,  int source, int tag)
{
  status s;
  container_proxy<T> tmp(out);
  MPI_Datatype type_id=get_mpi_datatype(*tmp.data());
  MPI_Recv(tmp.data(),tmp.size(),type_id,source,tag,comm,&s);
}

template<typename T>
inline request isend(const communicator& comm, T& in,  int dest, int tag)
{
  request r;
  container_proxy<T> tmp(in);
  MPI_Datatype type_id=get_mpi_datatype(*tmp.data());
  int ierr=MPI_Isend(tmp.data(),tmp.size(),type_id,dest,tag,comm,&r);
  return r;
}

template<typename T>
inline request irecv(const communicator& comm, T& out,  int source, int tag)
{
  request r;
  container_proxy<T> tmp(out);
  MPI_Datatype type_id=get_mpi_datatype(*tmp.data());
  int ierr=MPI_Irecv(tmp.data(),tmp.size(),type_id,source,tag,comm,&r);
  return r;
}

template<typename T>
inline void send(const communicator& comm, T* x,  int n,  int dest, int tag)
{
  MPI_Send(x,n,get_mpi_datatype(*x),dest,tag,comm);
}

template<typename T>
inline void recv(const communicator& comm, T* x,  int n,  int source, int tag)
{
  status s;
  int ierr=MPI_Recv(x,n,get_mpi_datatype(*x),source,tag,comm,&s);
}

template<typename T>
inline request isend(const communicator& comm, T* x,  int n,  int dest, int tag)
{
  request r;
  MPI_Isend(x,n,get_mpi_datatype(*x),dest,tag,comm,&r);
  return r;
}

template<typename T>
inline request irecv(const communicator& comm, T* x,  int n,  int source, int tag)
{
  request r;
  MPI_Irecv(x,n,get_mpi_datatype(*x),source,tag,comm,&r);
  return r;
}

#else
template<typename T>
inline void send(const communicator& comm, T& in,  int dest, int tag)
{
}

template<typename T>
inline void recv(const communicator& comm, T& out,  int source, int tag)
{
}

template<typename T>
inline request isend(const communicator& comm, T& in,  int dest, int tag)
{
  return request;
}
template<typename T>
inline request irecv(const communicator& comm, T& in,  int source, int tag)
{
  return request;
}
template<typename T>
inline void send(const communicator& comm, T* x,  int n,  int dest, int tag)
{
}

template<typename T>
inline void recv(const communicator& comm, T* x,  int n,  int source, int tag)
{
}
template<typename T>
inline request isend(const communicator& comm, T* x,  int n,  int dest, int tag)
{
  return request;
}

template<typename T>
inline request irecv(const communicator& comm, T* x,  int n,  int source, int tag)
{
  return request;
}
#endif
}
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2458 $   $Date: 2008-02-20 10:45:51 -0500 (Wed, 20 Feb 2008) $
 * $Id: CommOperators.h 2458 2008-02-20 15:45:51Z jnkim $
 ***************************************************************************/