//
// Functions for MPI communications
//

#include <cassert>
#include "msg.h"
#include "comm.h"

namespace {
  int this_node= -1;
  int n_nodes= 0;
  int parallel_level= 0;
    // 0: no MPI
    // 1: MPI_THREAD_SINGLE, only one thread per MPI node
    // 2: MPI_THREAD_FUNNELED, only the thread that called MPI_Init_thread will
    //    make MPI calls.

  CommStatus mpi_status= comm_uninitialised;
}

//
// Initialisation
//

void comm_mpi_init(int* p_argc, char*** p_argv)
{
#ifdef _OPENMP
  int thread_level;
  MPI_Init_thread(p_argc, p_argv, MPI_THREAD_FUNNELED, &thread_level);
  if(thread_level >= MPI_THREAD_FUNNELED)
    parallel_level= 2;
  else
    parallel_level= 1;
#else
  int ret= MPI_Init(p_argc, p_argv); assert(ret == MPI_SUCCESS);
  parallel_level= 1;
#endif

  mpi_status= comm_parallel;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &this_node);
  MPI_Comm_size(MPI_COMM_WORLD, &n_nodes);
}

void comm_mpi_msg()
{
  msg_printf(msg_verbose,
	     "MPI initialised with n_nodes= %d, parallel_level= %d.\n",
	     n_nodes, parallel_level);

}

void comm_mpi_finalise()
{
  msg_printf(msg_verbose, "Finishing program with MPI_Finalize()\n");
  mpi_status= comm_finalised;

  MPI_Finalize();
}

int comm_this_node(void)
{
  return this_node;
}

int comm_n_nodes(void)
{
  return n_nodes;
}

void comm_abort(void)
{
  MPI_Abort(MPI_COMM_WORLD, 1);
}

void comm_barrier()
{
  MPI_Barrier(MPI_COMM_WORLD);
}


void comm_bcast_int(int* p_int, int count)
{
  MPI_Bcast(p_int, count, MPI_INT, 0, MPI_COMM_WORLD);
}

void comm_bcast_double(double* p_double, int count)
{
  MPI_Bcast(p_double, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

CommStatus comm_status()
{
  return mpi_status;
}
