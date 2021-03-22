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
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Manager.hxx"

namespace mfem_mgis {

#ifdef MFEM_USE_MPI
  template<>
  Manager<true>::Manager() {};

  template<>
  bool Manager<true>::Initialized() const {
    int flag;
    auto const error = MPI_Initialized(&flag);
    if (error != MPI_SUCCESS) {
      throw std::runtime_error("Error while calling MPI_Initialized");
    }
    return static_cast<bool>(flag);
  }
  
  template<>
  bool Manager<true>::Finalized() const {
    int finalized;
    MPI_Finalized(&finalized);
    return finalized;
  }

  template<>
  void Manager<true>::Init(int *argc, char ***argv)  const {
    if (Finalized())
      throw std::runtime_error("MPI session has already been finalized");
    if (not Initialized()) {
      int retcode = (MPI_Init(argc, const_cast<char ***>(argv)));
      if (retcode != MPI_SUCCESS)
	throw std::runtime_error("MPI session failed to initialize");
    }
  }

  template<>
  Manager<true>::Manager(int *argc, char ***argv) {
    Init(argc, argv);
  }
  
  template<>
  void Manager<true>::Init()  const{
    Init(NULL,NULL);
  }
  
  template<>
  void Manager<true>::Finalize() const {
      if (Finalized() or not Initialized())
	return;
      MPI_Finalize();
  }

  template<>
  Manager<true>::~Manager() = default;
  

#endif /* MFEM_USE_MPI */
    
  template<>
  Manager<false>::Manager() {};

  template<>
  bool Manager<false>::Initialized() const {
    return true;
  }
  
  template<>
  bool Manager<false>::Finalized() const {
    return false;
  }
  
  template<>
  void Manager<false>::Init(int *argc, char ***argv) const {
    (void)argc;
    (void)argv;
  }

  template<>
  Manager<false>::Manager(int *argc, char ***argv) {
    Init(argc, argv);
  }
  
  template<>
  void Manager<false>::Init() const {}
    
  template<>
  void Manager<false>::Finalize() const { }

  template<>
  Manager<false>::~Manager() = default;
 
} // namespace mfem_mgis
