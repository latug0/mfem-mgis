/*!
 * \file   src/PartialQuadratureSpace.cxx
 * \brief
 * \author Thomas Helfer
 * \date   8/06/2020
 */

#include <cmath>
#include <iterator>
#include <algorithm>
#include <mfem/fem/fespace.hpp>
#include <mfem/fem/pfespace.hpp>
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"

namespace mfem_mgis {

  PartialQuadratureSpace::PartialQuadratureSpace(
      const FiniteElementDiscretization& fed,
      const size_type m,
      const std::function<const mfem::IntegrationRule&(
          const mfem::FiniteElement&, const mfem::ElementTransformation&)>&
          integration_rule_selector)
      : fe_discretization(fed), id(m) {
    const auto& fespace = this->fe_discretization.getFiniteElementSpace();
    this->ng = size_type{};
    for (size_type i = 0; i != fespace.GetNE(); ++i) {
      if (fespace.GetAttribute(i) != m) {
        continue;
      }
      const auto& fe = *(fespace.GetFE(i));
      const auto& tr = *(fespace.GetElementTransformation(i));
      this->offsets[i] = this->ng;
      const auto& ir = integration_rule_selector(fe, tr);
      this->ng += ir.GetNPoints();
    }
  }  // end of PartialQuadratureSpace::PartialQuadratureSpace

  size_type PartialQuadratureSpace::getNumberOfElements() const {
    return this->offsets.size();
  }  // end of PartialQuadratureSpace::getNumberOfElements

  size_type PartialQuadratureSpace::getNumberOfIntegrationPoints() const {
    return this->ng;
  }  // end of PartialQuadratureSpace::getNumberOfIntegrationPoints

  size_type PartialQuadratureSpace::getOffset(const size_type i) const {
    const auto p = this->offsets.find(i);
    if (p == this->offsets.end()) {
      mgis::raise(
          "PartialQuadratureSpace::getOffset: "
          "invalid element number '" +
          std::to_string(i) + "' for material '" + std::to_string(id) + "'");
    }
    return p->second;
  }  // end of PartialQuadratureSpace::getOffset

  PartialQuadratureSpace::~PartialQuadratureSpace() = default;

}  // end of namespace mfem_mgis
