#include <algorithm>
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MFEMMGIS/IsotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator.hxx"

namespace mfem_mgis {

  inline void
  IsotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
      updateGradients(mgis::span<real> &g,
                      const mfem::Vector &u,
                      const mfem::DenseMatrix &dN,
                      const size_type ni) noexcept {
    const auto nnodes = dN.NumRows();
    const auto dNi_0 = dN(ni, 0);
    const auto dNi_1 = dN(ni, 1);
    const auto dNi_2 = dN(ni, 2);
    const auto u_0 = u[ni];
    const auto u_1 = u[ni + nnodes];
    const auto u_2 = u[ni + 2 * nnodes];
    g[0] += u_0 * dNi_0;
    g[1] += u_1 * dNi_1;
    g[2] += u_2 * dNi_2;
    g[3] += dNi_1 * u_0;
    g[4] += u_1 * dNi_0;
    g[5] += dNi_2 * u_0;
    g[6] += u_2 * dNi_0;
    g[7] += u_1 * dNi_2;
    g[8] += dNi_1 * u_2;
  }  // end of updateGradients

  inline void
  IsotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
      updateInnerForces(mfem::Vector &Fe,
                        const mgis::span<const real> &s,
                        const mfem::DenseMatrix &dN,
                        const real w,
                        const size_type ni) const noexcept {
    const auto nnodes = dN.NumRows();
    const auto dNi_0 = dN(ni, 0);
    const auto dNi_1 = dN(ni, 1);
    const auto dNi_2 = dN(ni, 2);
    const auto ni_0 = ni;
    const auto ni_1 = ni + nnodes;
    const auto ni_2 = ni + 2 * nnodes;
    Fe[ni_0] += w * (s[0] * dNi_0 + dNi_1 * s[3] + s[5] * dNi_2);
    Fe[ni_1] += w * (dNi_2 * s[7] + s[1] * dNi_1 + s[4] * dNi_0);
    Fe[ni_2] += w * (s[2] * dNi_2 + dNi_1 * s[8] + s[6] * dNi_0);
  }  // end of updateInnerForces

  inline void
  IsotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
      updateStiffnessMatrix(mfem::DenseMatrix &Ke,
                            const mgis::span<const real> &Kip,
                            const mfem::DenseMatrix &dN,
                            const real w,
                            const size_type ni) const noexcept {
    const auto nnodes = dN.NumRows();
    const auto dNi_0 = dN(ni, 0);
    const auto dNi_1 = dN(ni, 1);
    const auto dNi_2 = dN(ni, 2);
    const auto ni_0 = ni;
    const auto ni_1 = ni + nnodes;
    const auto ni_2 = ni + 2 * nnodes;
    for (size_type nj = 0; nj != nnodes; ++nj) {
      const auto dNj_0 = dN(nj, 0);
      const auto dNj_1 = dN(nj, 1);
      const auto dNj_2 = dN(nj, 2);
      const auto nj_0 = nj;
      const auto nj_1 = nj + nnodes;
      const auto nj_2 = nj + 2 * nnodes;
      Ke(ni_0, nj_0) += w * (Kip[5] * dNi_0 * dNj_2 + dNj_1 * dNi_1 * Kip[30] +
                             dNi_2 * dNj_0 * Kip[45] + Kip[27] * dNj_0 * dNi_1 +
                             Kip[32] * dNj_2 * dNi_1 + dNi_2 * dNj_1 * Kip[48] +
                             dNj_1 * dNi_0 * Kip[3] + dNi_2 * Kip[50] * dNj_2 +
                             dNj_0 * Kip[0] * dNi_0);
      Ke(ni_0, nj_1) += w * (Kip[46] * dNi_2 * dNj_1 + dNj_1 * dNi_0 * Kip[1] +
                             Kip[34] * dNj_2 * dNi_1 + dNi_0 * dNj_2 * Kip[7] +
                             dNj_0 * Kip[31] * dNi_1 + Kip[28] * dNj_1 * dNi_1 +
                             dNi_2 * dNj_2 * Kip[52] + Kip[49] * dNi_2 * dNj_0 +
                             dNj_0 * dNi_0 * Kip[4]);
      Ke(ni_0, nj_2) += w * (dNj_0 * dNi_0 * Kip[6] + dNi_2 * Kip[47] * dNj_2 +
                             Kip[2] * dNi_0 * dNj_2 + Kip[29] * dNj_2 * dNi_1 +
                             dNj_0 * dNi_1 * Kip[33] + dNi_2 * Kip[53] * dNj_1 +
                             dNj_1 * Kip[35] * dNi_1 + dNi_2 * dNj_0 * Kip[51] +
                             Kip[8] * dNj_1 * dNi_0);
      Ke(ni_1, nj_0) += w * (Kip[12] * dNj_1 * dNi_1 + dNj_0 * Kip[9] * dNi_1 +
                             dNj_1 * dNi_0 * Kip[39] + Kip[68] * dNi_2 * dNj_2 +
                             dNi_2 * dNj_0 * Kip[63] + dNi_2 * Kip[66] * dNj_1 +
                             Kip[36] * dNj_0 * dNi_0 + dNi_0 * Kip[41] * dNj_2 +
                             Kip[14] * dNj_2 * dNi_1);
      Ke(ni_1, nj_1) += w * (dNi_2 * dNj_0 * Kip[67] + dNj_0 * Kip[40] * dNi_0 +
                             dNj_1 * Kip[10] * dNi_1 + Kip[43] * dNi_0 * dNj_2 +
                             dNj_0 * Kip[13] * dNi_1 + Kip[37] * dNj_1 * dNi_0 +
                             dNi_2 * dNj_1 * Kip[64] + dNi_2 * Kip[70] * dNj_2 +
                             Kip[16] * dNj_2 * dNi_1);
      Ke(ni_1, nj_2) += w * (dNj_1 * dNi_1 * Kip[17] + Kip[11] * dNj_2 * dNi_1 +
                             dNi_0 * Kip[38] * dNj_2 + dNi_2 * dNj_0 * Kip[69] +
                             Kip[65] * dNi_2 * dNj_2 + dNj_1 * Kip[44] * dNi_0 +
                             Kip[71] * dNi_2 * dNj_1 + dNj_0 * dNi_0 * Kip[42] +
                             dNj_0 * Kip[15] * dNi_1);
      Ke(ni_2, nj_0) += w * (dNj_0 * dNi_0 * Kip[54] + dNj_0 * Kip[72] * dNi_1 +
                             dNi_2 * dNj_0 * Kip[18] + dNi_2 * Kip[21] * dNj_1 +
                             Kip[75] * dNj_1 * dNi_1 + dNi_2 * dNj_2 * Kip[23] +
                             dNj_2 * dNi_1 * Kip[77] + dNj_1 * Kip[57] * dNi_0 +
                             Kip[59] * dNi_0 * dNj_2);
      Ke(ni_2, nj_1) += w * (dNi_2 * Kip[25] * dNj_2 + dNj_0 * Kip[76] * dNi_1 +
                             dNi_2 * dNj_0 * Kip[22] + dNj_1 * dNi_0 * Kip[55] +
                             dNi_0 * Kip[61] * dNj_2 + dNj_0 * dNi_0 * Kip[58] +
                             Kip[79] * dNj_2 * dNi_1 + dNi_2 * dNj_1 * Kip[19] +
                             dNj_1 * Kip[73] * dNi_1);
      Ke(ni_2, nj_2) += w * (Kip[56] * dNi_0 * dNj_2 + Kip[62] * dNj_1 * dNi_0 +
                             dNj_0 * Kip[78] * dNi_1 + Kip[24] * dNi_2 * dNj_0 +
                             dNj_1 * dNi_1 * Kip[80] + dNi_2 * dNj_1 * Kip[26] +
                             dNj_0 * Kip[60] * dNi_0 + Kip[74] * dNj_2 * dNi_1 +
                             dNi_2 * dNj_2 * Kip[20]);
    }  // end of for (size_type nj = 0; nj != nnodes; ++nj)
  }    // end of updateStiffnessMatrix

  IsotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
      IsotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator(
          const FiniteElementDiscretization &fed,
          const size_type m,
          std::unique_ptr<const Behaviour> b_ptr)
      : StandardBehaviourIntegratorCRTPBase<
            IsotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator>(
            buildQuadratureSpace(fed, m), std::move(b_ptr)) {
    if (this->b.symmetry != Behaviour::ISOTROPIC) {
      mgis::raise("invalid behaviour symmetry");
    }
  }  // end of
     // IsotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator

  IsotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
      RotationMatrix
      IsotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
          getRotationMatrix(const size_type) const {
    return RotationMatrix{};
  }  // end of getRotationMatrix

  void IsotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
      rotateGradients(mgis::span<real>, const RotationMatrix &) {
  }  // end of rotateGradients

  mgis::span<const real>
  IsotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
      rotateThermodynamicForces(mgis::span<const real> s,
                                const RotationMatrix &) {
    return s;
  }

  void IsotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
      rotateTangentOperatorBlocks(mgis::span<real>, const RotationMatrix &) {}

  void IsotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
      computeInnerForces(mfem::Vector &Fe,
                         const mfem::FiniteElement &e,
                         mfem::ElementTransformation &tr,
                         const mfem::Vector &u) {
    this->implementComputeInnerForces(Fe, e, tr, u);
  }  // end of computeInnerForces

  void IsotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
      computeStiffnessMatrix(mfem::DenseMatrix &Ke,
                             const mfem::FiniteElement &e,
                             mfem::ElementTransformation &tr,
                             const mfem::Vector &) {
    this->implementComputeStiffnessMatrix(Ke, e, tr);
  }  // end of computeStiffnessMatrix

  const mfem::IntegrationRule &
  IsotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
      getIntegrationRule(const mfem::FiniteElement &el,
                         const mfem::ElementTransformation &Trans) {
    const auto order = 2 * Trans.OrderGrad(&el);
    return mfem::IntRules.Get(el.GetGeomType(), order);
  }

  std::shared_ptr<const PartialQuadratureSpace>
  IsotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
      buildQuadratureSpace(const FiniteElementDiscretization &fed,
                           const size_type m) {
    auto selector = [](const mfem::FiniteElement &e,
                       const mfem::ElementTransformation &tr)
        -> const mfem::IntegrationRule & {
      return getIntegrationRule(e, tr);
    };  // end of selector
    return std::make_shared<PartialQuadratureSpace>(fed, m, selector);
  }  // end of buildQuadratureSpace

  IsotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
      ~IsotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator() =
          default;

}  // end of namespace mfem_mgis
