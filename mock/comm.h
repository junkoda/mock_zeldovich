#ifndef COMM_H
#define COMM_H 1

#include <iostream>
#include <typeinfo>
#include <mpi.h>
#include "msg.h"
#include "error.h"

enum CommStatus {comm_uninitialised=0, comm_parallel=1, comm_finalised=2};

void comm_mpi_init(int* pargc, char*** pargv);
void comm_mpi_finalise();
void comm_mpi_msg();

int  comm_this_node();
int  comm_n_nodes();

void comm_abort();
void comm_barrier();

void comm_bcast_int(int* p_int, int count);
void comm_bcast_double(double* p_double, int count);

CommStatus comm_status();


static inline MPI_Datatype mpi_datatype(const std::type_info& type_id)
{
  if(type_id == typeid(int))
    return MPI_INT;
  else if(type_id == typeid(long))
    return MPI_LONG;
  else if(type_id == typeid(long long))
    return MPI_LONG_LONG;
  else if(type_id == typeid(float))
    return MPI_FLOAT;
  else if(type_id == typeid(double))
    return MPI_DOUBLE;

  msg_printf(msg_fatal, "Error: unknown data type\n");
  throw RuntimeError();

  return 0;
}


template<class T> void comm_bcast(T x)
{
  MPI_Bcast(&x, 1, mpi_datatype(typeid(T)), 0, MPI_COMM_WORLD);
}

template<class T> T comm_sum(T x)
{
  T x_reduced;
  MPI_Allreduce(&x, &x_reduced, 1, mpi_datatype(typeid(T)),
		MPI_SUM, MPI_COMM_WORLD);
  //MPI_Reduce(&x, &x_reduced,1,  mpi_datatype(typeid(T)), MPI_SUM, 0, MPI_COMM_WORLD);
  //MPI_Bcast(&x_reduced, 1, mpi_datatype(typeid(T)), 0, MPI_COMM_WORLD);

  return x_reduced;
}

template<class T> T comm_max(T x)
{
  T x_reduced;
  MPI_Allreduce(&x, &x_reduced, 1, mpi_datatype(typeid(T)),
		MPI_MAX, MPI_COMM_WORLD);

  return x_reduced;
}


template<class T> T comm_partial_sum(T x)
{
  T x_reduced;
  MPI_Scan(&x, &x_reduced, 1, mpi_datatype(typeid(T)),
	   MPI_SUM, MPI_COMM_WORLD);

  return x_reduced;
}


#endif
