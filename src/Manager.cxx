/*!
 * \file   src/Manager.cxx

 * \brief
 * Manage start and death of a given simulation.
 *
 * \author Guillaume Latu
 * \date   19/03/2020
 */
#include <exception>
#include <utility>
#include <iostream>
#include "MFEMMGIS/Manager.hxx"
using namespace std;

namespace mfem_mgis {

  namespace Manager {

    void Init() {
      Init(NULL,NULL);
    }
    
    void Init(int *argc, char ***argv) {
      if (Finalized())
	throw std::runtime_error("MPI session has already been finalized");
      if (not Initialized()) {
	int retcode = (MPI_Init(argc, const_cast<char ***>(argv)));
	if (retcode != MPI_SUCCESS)
	  throw std::runtime_error("MPI session failed to initialize");
      }
    }

//    void Broadcast(int *argc, const char ***argv) {
//      /* TODO */
//    }

    bool Initialized() {
      int flag;
      auto const error = MPI_Initialized(&flag);
      if (error != MPI_SUCCESS) {
	throw std::runtime_error("Error while calling MPI_Initialized");
      }
      return static_cast<bool>(flag);
    }
    
    bool Finalized() {
      int finalized;
      MPI_Finalized(&finalized);
      return finalized;
    }
    
    void Finalize() {
      if (Finalized() or not Initialized()) return;
      MPI_Finalize();
    }
  }
} // namespace mfem_mgis
