/*!
 * \file   src/FiniteElementDiscretization.cxx
 * \brief
 * \author Thomas Helfer
 * \date 16/12/2020
 */

#include <utility>
#include <mfem/fem/fespace.hpp>
#include <mfem/fem/pfespace.hpp>
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"

namespace mfem_mgis {

  FiniteElementDiscretization::FiniteElementDiscretization(
      std::shared_ptr<mMesh> m,
      std::shared_ptr<const mfem::FiniteElementCollection> c,
      const size_type d)
      : mesh(std::move(m)), fec(std::move(c)) {
    this->fe_space = std::make_unique<mFiniteElementSpace>(
        this->mesh.get(), this->fec.get(), d);
  }  // end of FiniteElementDiscretization

  FiniteElementDiscretization::FiniteElementDiscretization(
      std::shared_ptr<mMesh> m,
      std::shared_ptr<const mfem::FiniteElementCollection> c,
      std::unique_ptr<mFiniteElementSpace> s)
      : mesh(std::move(m)), fec(std::move(c)), fe_space(std::move(s)) {
    if (this->mesh.get() != this->fe_space->GetMesh()) {
      mgis::raise(
          "FiniteElementDiscretization::FiniteElementDiscretization: "
          "mesh pointer don't match the mesh on which the finite element space "
          "is built");
    }
  }  // end of FiniteElementDiscretization

  mMesh& FiniteElementDiscretization::getMesh() {
    return *(this->mesh);
  }  // end of getMesh

  const mMesh& FiniteElementDiscretization::getMesh() const {
    return *(this->mesh);
  }  // end of getMesh

  mFiniteElementSpace&
  FiniteElementDiscretization::getFiniteElementSpace() {
    return *(this->fe_space);
  }  // end of getFiniteElementSpace

  const mFiniteElementSpace&
  FiniteElementDiscretization::getFiniteElementSpace() const {
    return *(this->fe_space);
  }  // end of getFiniteElementSpace

  const mfem::FiniteElementCollection&
  FiniteElementDiscretization::getFiniteElementCollection() const {
    return *(this->fec);
  }  // end of getFiniteElementCollection

  FiniteElementDiscretization::~FiniteElementDiscretization() = default;

} // end of namespace mfem_mgis
