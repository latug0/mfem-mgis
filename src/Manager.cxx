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

namespace mfem_mgis {

  namespace Manager {

    template<>
    bool Initialized<true>() {
      int flag;
      auto const error = MPI_Initialized(&flag);
      if (error != MPI_SUCCESS) {
	throw std::runtime_error("Error while calling MPI_Initialized");
      }
      return static_cast<bool>(flag);
    }

    template<>
    bool Initialized<false>() {
      return true;
    }

    template<>
    bool Finalized<true>() {
      int finalized;
      MPI_Finalized(&finalized);
      return finalized;
    }
    
    template<>
    bool Finalized<false>() {
      return false;
    }

    template<>
    void Init<true>(int *argc, char ***argv) {
      if (Finalized<true>())
	throw std::runtime_error("MPI session has already been finalized");
      if (not Initialized<true>()) {
	int retcode = (MPI_Init(argc, const_cast<char ***>(argv)));
	if (retcode != MPI_SUCCESS)
	  throw std::runtime_error("MPI session failed to initialize");
      }
    }

    template<>
    void Init<false>(int *argc, char ***argv) {}

    template<>
    void Init<true>() {
      Init<true>(NULL,NULL);
    }
    template<>
    void Init<false>() {}
    
    
    //    void Broadcast(int *argc, const char ***argv) {
    //      /* TODO */
    //    }
    
    
    template<>
    void Finalize<true>() {
      if (Finalized<true>() or not Initialized<true>())
	return;
      MPI_Finalize();
    }

    template<>
    void Finalize<false>() { }
  }
} // namespace mfem_mgis
