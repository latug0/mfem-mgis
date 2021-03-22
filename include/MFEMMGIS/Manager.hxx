/*!
 * \file   include/MFEMMGIS/Manager.hxx

 * \brief class designed to handling start and death of a given simulation, taking
 *        in charge the MPI parallelism.
 *
 * \author Guillaume Latu
 * \date   19/03/2020
 */

#ifndef LIB_MFEM_MGIS_MANAGER_HXX
#define LIB_MFEM_MGIS_MANAGER_HXX

#ifdef MFEM_USE_MPI
#include <mpi.h>
#else
#define MPI_COMM_WORLD 
//#define MPI_SUCCESS 0
//#define MPI_Comm_rank(args...) (0)
//#define MPI_Comm_size(args...) (1)
#endif /* MFEM_USE_MPI */
#include <iostream>

namespace mfem_mgis {

  class ManagerBase {
  public :

    /*!
     * \brief Initialize Manager and distributed parallelism (if needed).
     */
    virtual void Init() const = 0;

    /*!
     * \brief Initialize Manager and distributed parallelism (if needed).
     *
     * Compared the `Init` function without parameters, this function
     * handles and broadcast parameters `argc` and `argv` to all MPI
     * processes (if the setting is not serial).
     *
     * \param[in] argc: number of command line parameters
     * \param[in] argv: array of strings that stores input parameters
     */
    virtual void Init(int *argc, char ***argv) const = 0;

    /*!
     * \brief Terminate the Manager and close MPI properly.
     */
    virtual void Finalize() const = 0;

    /*!
     * \brief Terminate the Manager and close MPI properly with a final call to exit.
     *
     * Optional arguments permit to print a `message` on console, and to
     * exit the program with an integer return code `errcode`; 
     *
     * \param[in] message: message to be printed at screen
     * \param[in] errcode: return code
     */
    virtual void FinalizeExit (const std::string &message={}, const int &errcode=EXIT_FAILURE) const = 0;

    /*!
     * \brief Terminate the Manager and close MPI properly with a final call to exit.
     *
     * Optional arguments permit to print a `message` on console, and to
     * exit the program with an integer return code `errcode`.
     *
     * \param[in] message: message to be printed at screen
     * \param[in] errcode: return code
     */
    virtual void Abort (const std::string &message={}, const int &errcode=EXIT_FAILURE) const = 0;
    
    virtual ~ManagerBase() = default;
    /* !
       \brief member variable that stores if Manager handles MPI parallelism or not.
    */
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
