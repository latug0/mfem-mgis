#include <algorithm>
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MFEMMGIS/IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator.hxx"

namespace mfem_mgis {

  inline void
  IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator::
      updateGradients(mgis::span<real> &g,
                      const mfem::Vector &u,
                      const mfem::DenseMatrix &dN,
                      const size_type ni) noexcept {
    const auto nnodes = dN.NumRows();
    const auto dNi_0 = dN(ni, 0);
    const auto dNi_1 = dN(ni, 1);
    const auto u_0 = u[ni];
    const auto u_1 = u[ni + nnodes];
    g[0] += dNi_0 * u_0;
    g[1] += u_1 * dNi_1;
    g[2] += 0;
    g[3] += dNi_1 * u_0;
    g[4] += u_1 * dNi_0;
  }  // end of updateGradients

  inline void
  IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator::
      updateInnerForces(mfem::Vector &Fe,
                        const mgis::span<const real> &s,
                        const mfem::DenseMatrix &dN,
                        const real w,
                        const size_type ni) const noexcept {
    const auto nnodes = dN.NumRows();
    const auto dNi_0 = dN(ni, 0);
    const auto dNi_1 = dN(ni, 1);
    const auto ni_0 = ni;
    const auto ni_1 = ni + nnodes;
    Fe[ni_0] += w * (s[0] * dNi_0 + dNi_1 * s[3]);
    Fe[ni_1] += w * (s[4] * dNi_0 + s[1] * dNi_1);
  }  // end of updateInnerForces

  inline void
  IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator::
      updateStiffnessMatrix(mfem::DenseMatrix &Ke,
                            const mgis::span<const real> &Kip,
                            const mfem::DenseMatrix &dN,
                            const real w,
                            const size_type ni) const noexcept {
    const auto nnodes = dN.NumRows();
    const auto dNi_0 = dN(ni, 0);
    const auto dNi_1 = dN(ni, 1);
    const auto ni_0 = ni;
    const auto ni_1 = ni + nnodes;
    for (size_type nj = 0; nj != nnodes; ++nj) {
      const auto dNj_0 = dN(nj, 0);
      const auto dNj_1 = dN(nj, 1);
      const auto nj_0 = nj;
      const auto nj_1 = nj + nnodes;
      Ke(ni_0, nj_0) += w * (dNj_0 * Kip[0] * dNi_0 + dNj_1 * dNi_1 * Kip[18] +
                             Kip[15] * dNj_0 * dNi_1 + Kip[3] * dNj_1 * dNi_0);
      Ke(ni_0, nj_1) += w * (dNj_0 * Kip[19] * dNi_1 + Kip[16] * dNj_1 * dNi_1 +
                             dNj_0 * Kip[4] * dNi_0 + dNj_1 * Kip[1] * dNi_0);
      Ke(ni_1, nj_0) += w * (dNj_1 * Kip[8] * dNi_1 + dNj_1 * Kip[23] * dNi_0 +
                             dNj_0 * dNi_0 * Kip[20] + dNj_0 * Kip[5] * dNi_1);
      Ke(ni_1, nj_1) += w * (dNj_0 * dNi_0 * Kip[24] + dNj_0 * dNi_1 * Kip[9] +
                             Kip[6] * dNj_1 * dNi_1 + dNj_1 * dNi_0 * Kip[21]);
    }  // end of for (size_type nj = 0; nj != nnodes; ++nj)
  }    // end of updateStiffnessMatrix

  IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator::
      IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator(
          const FiniteElementDiscretization &fed,
          const size_type m,
          std::unique_ptr<const Behaviour> b_ptr)
      : StandardBehaviourIntegratorCRTPBase<
            IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator>(
            buildQuadratureSpace(fed, m), std::move(b_ptr)) {
    if (this->b.symmetry != Behaviour::ISOTROPIC) {
      mgis::raise("invalid behaviour symmetry");
    }
  }  // end of
     // IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator

  void IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator::
      setRotationMatrix(const RotationMatrix2D &) {
    mgis::raise(
        "IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator::"
        "setRotationMatrix: invalid call");
  }

  void IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator::
      setRotationMatrix(const RotationMatrix3D &) {
    mgis::raise(
        "IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator::"
        "setRotationMatrix: invalid call");
  }

  IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator::
      RotationMatrix
      IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator::
          getRotationMatrix() const {
    return RotationMatrix{};
  }  // end of getRotationMatrix

  void IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator::
      rotateGradients(mgis::span<real>, const RotationMatrix &) {
  }  // end of rotateGradients

  mgis::span<const real>
  IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator::
      rotateThermodynamicForces(mgis::span<const real> s,
                                const RotationMatrix &) {
    return s;
  }

  void IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator::
      rotateTangentOperatorBlocks(mgis::span<real>, const RotationMatrix &) {}

  void IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator::
      computeInnerForces(mfem::Vector &Fe,
                         const mfem::FiniteElement &e,
                         mfem::ElementTransformation &tr,
                         const mfem::Vector &u) {
    this->implementComputeInnerForces(Fe, e, tr, u);
  }  // end of computeInnerForces

  void IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator::
      computeStiffnessMatrix(mfem::DenseMatrix &Ke,
                             const mfem::FiniteElement &e,
                             mfem::ElementTransformation &tr,
                             const mfem::Vector &) {
    this->implementComputeStiffnessMatrix(Ke, e, tr);
  }  // end of computeStiffnessMatrix

  const mfem::IntegrationRule &
  IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator::
      getIntegrationRule(const mfem::FiniteElement &el,
                         const mfem::ElementTransformation &Trans) {
    const auto order = 2 * Trans.OrderGrad(&el);
    return mfem::IntRules.Get(el.GetGeomType(), order);
  }

  std::shared_ptr<const PartialQuadratureSpace>
  IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator::
      buildQuadratureSpace(const FiniteElementDiscretization &fed,
                           const size_type m) {
    auto selector = [](const mfem::FiniteElement &e,
                       const mfem::ElementTransformation &tr)
        -> const mfem::IntegrationRule & {
      return getIntegrationRule(e, tr);
    };  // end of selector
    return std::make_shared<PartialQuadratureSpace>(fed, m, selector);
  }  // end of buildQuadratureSpace

  IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator::
      ~IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator() =
          default;

}  // end of namespace mfem_mgis
