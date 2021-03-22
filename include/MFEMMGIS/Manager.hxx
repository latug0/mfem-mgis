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

#ifdef MFEM_USE_MPI
#include <mpi.h>
#else
#define MPI_COMM_WORLD 0
//#define MPI_SUCCESS 0
//#define MPI_Comm_rank(args...) (0)
//#define MPI_Comm_size(args...) (1)
#endif /* MFEM_USE_MPI */
#include <iostream>

namespace mfem_mgis {

  class ManagerBase {
  public :
    virtual void Init() const = 0;
    virtual void Init(int *argc, char ***argv) const = 0;
    virtual void Finalize() const = 0;
    virtual void FinalizeExit (const std::string &message={}, const int &errcode=EXIT_FAILURE) const = 0;
    virtual void Abort (const std::string &message={}, const int &errcode=EXIT_FAILURE) const = 0;
    virtual ~ManagerBase() = default;
    static bool isparallel;
  };
  
  template <bool parallel>
  class Manager : public ManagerBase {
  public:
    Manager();
    Manager(int *argc, char ***argv);
    void Init() const override;
    void Init(int *argc, char ***argv) const override;
    void Finalize() const override;
    void FinalizeExit (const std::string &message={}, const int &errcode=EXIT_FAILURE) const override;
    void Abort (const std::string &message={}, const int &errcode=EXIT_FAILURE) const override;
    //TODO: void Broadcast(int *argc, char ***argv);
    ~Manager() = default;

    static constexpr bool isparallel = parallel;

  protected :
    bool Initialized() const;
    bool Finalized() const;
  };
} // namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_MANAGER_HXX */
