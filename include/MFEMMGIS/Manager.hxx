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
//#warning "NO MPI activated"
//#define MPI_COMM_WORLD 0
//#define MPI_SUCCESS 0
//#define MPI_Comm_rank(args...) (0)
//#define MPI_Comm_size(args...) (1)
#endif /* DO_USE_MPI */

namespace mfem_mgis {
  namespace Manager {

    template <bool parallel>
    void Init();

    template <bool parallel>
    void Init(int *argc, char ***argv);

    template <bool parallel>
    void Broadcast(int *argc, char ***argv);

    template <bool parallel>
    bool Initialized();

    template <bool parallel>
    bool Finalized();

    template <bool parallel>
    void Finalize();
  }
} // namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_MANAGER_HXX */
