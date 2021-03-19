/*!
 * \file   include/MFEMMGIS/Manager.hxx

 * \brief
 * Manage start and death of a given simulation.
 *
 * \author Guillaume Latu
 * \date   19/03/2020
 */

#ifndef LIB_MFEM_MGIS_MANAGER_HXX
#define LIB_MFEM_MGIS_MANAGER_HXX

#ifdef DO_USE_MPI
#include <mpi.h>
#warning "MPI activated"
#else
#warning "NO MPI activated"
#define MPI_COMM_WORLD 0
#define MPI_SUCCESS 0
#define MPI_Finalize() 
#define MPI_Comm_rank(args...) (0)
#define MPI_Comm_size(args...) (1)
#define MPI_Init(args...) (MPI_SUCCESS)
#define MPI_Finalized(flag) (*(flag)=0, MPI_SUCCESS)
#define MPI_Initialized(flag) (*(flag)=1, MPI_SUCCESS)
#endif /* DO_USE_MPI */

namespace mfem_mgis {
  namespace Manager {
    void Init();
    void Init(int *argc, char ***argv);
    void Broadcast(int *argc, char ***argv);
    bool Initialized();
    bool Finalized();
    void Finalize();
  }
} // namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_MANAGER_HXX */
