/*!
 * \file   src/NonLinearEvolutionProblemImplementationBase.cxx
 * \brief
 * \author Thomas Helfer
 * \date   11/12/2020
 */

#include <iostream>
#include <utility>
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/NewtonSolver.hxx"
#include "MFEMMGIS/IntegrationType.hxx"
#include "MFEMMGIS/SolverUtilities.hxx"
#include "MFEMMGIS/DirichletBoundaryCondition.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/MultiMaterialNonLinearIntegrator.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementationBase.hxx"

namespace mfem_mgis {

  const char* const NonLinearEvolutionProblemImplementationBase::
      UseMultiMaterialNonLinearIntegrator =
          "UseMultiMaterialNonLinearIntegrator";

  std::vector<std::string>
  NonLinearEvolutionProblemImplementationBase::getParametersList() {
    return {NonLinearEvolutionProblemImplementationBase::
                UseMultiMaterialNonLinearIntegrator};
  }  // end of getParametersList

  MultiMaterialNonLinearIntegrator* buildMultiMaterialNonLinearIntegrator(
      std::shared_ptr<FiniteElementDiscretization> fed,
      const Hypothesis h,
      const Parameters& p) {
    const auto* const n = NonLinearEvolutionProblemImplementationBase::
        UseMultiMaterialNonLinearIntegrator;
    if (contains(p, n)) {
      if (!get<bool>(p, n)) {
        return nullptr;
      }
    }
    return new MultiMaterialNonLinearIntegrator(fed, h);
  }  // end of buildMultiMaterialNonLinearIntegrator

  NonLinearEvolutionProblemImplementationBase::
      NonLinearEvolutionProblemImplementationBase(
          std::shared_ptr<FiniteElementDiscretization> fed,
          const Hypothesis h,
          const Parameters& p)
      : fe_discretization(fed),
        u0(getTrueVSize(*fed)),
        u1(getTrueVSize(*fed)),
        mgis_integrator(buildMultiMaterialNonLinearIntegrator(fed, h, p)),
        hypothesis(h) {
    this->u0 = real{0};
    this->u1 = real{0};
  }  // end of NonLinearEvolutionProblemImplementationBase

  FiniteElementDiscretization& NonLinearEvolutionProblemImplementationBase::
      getFiniteElementDiscretization() {
    return *(this->fe_discretization);
  }  // end of getFiniteElementDiscretization

  std::shared_ptr<FiniteElementDiscretization>
  NonLinearEvolutionProblemImplementationBase::
      getFiniteElementDiscretizationPointer() {
    return this->fe_discretization;
  }  // end of getFiniteElementDiscretization

  mfem::Vector& NonLinearEvolutionProblemImplementationBase::
      getUnknownsAtBeginningOfTheTimeStep() {
    return this->u0;
  }  // end of getUnknownsAtBeginningOfTheTimeStep

  const mfem::Vector& NonLinearEvolutionProblemImplementationBase::
      getUnknownsAtBeginningOfTheTimeStep() const {
    return this->u0;
  }  // end of getUnknownsAtBeginningOfTheTimeStep

  mfem::Vector&
  NonLinearEvolutionProblemImplementationBase::getUnknownsAtEndOfTheTimeStep() {
    return this->u1;
  }  // end of getUnknownsAtEndOfTheTimeStep

  const mfem::Vector&
  NonLinearEvolutionProblemImplementationBase::getUnknownsAtEndOfTheTimeStep()
      const {
    return this->u1;
  }  // end of getUnknownsAtEndOfTheTimeStep

  void NonLinearEvolutionProblemImplementationBase::revert() {
    this->u1 = this->u0;
    if (this->mgis_integrator != nullptr) {
      this->mgis_integrator->revert();
    }
  }  // end of revert

  void NonLinearEvolutionProblemImplementationBase::update() {
    this->u0 = this->u1;
    if (this->mgis_integrator != nullptr) {
      this->mgis_integrator->update();
    }
  }  // end of update

  static void checkMultiMaterialSupportEnabled(
      const char* const n, const MultiMaterialNonLinearIntegrator* const p) {
    if (p == nullptr) {
      std::string msg("NonLinearEvolutionProblemImplementationBase::");
      msg += n;
      msg += ": multi material support has been disabled";
      raise(msg);
    }
  }  // end of checkMultiMaterialSupportEnabled

  std::vector<size_type>
  NonLinearEvolutionProblemImplementationBase::getMaterialIdentifiers() const {
    checkMultiMaterialSupportEnabled("getMaterialIdentifiers",
                                     this->mgis_integrator);
    return this->mgis_integrator->getMaterialIdentifiers();
  }  // end of getMaterialIdentifiers

  void NonLinearEvolutionProblemImplementationBase::addBehaviourIntegrator(
      const std::string& n,
      const size_type m,
      const std::string& l,
      const std::string& b) {
    checkMultiMaterialSupportEnabled("addBehaviourIntegrator",
                                     this->mgis_integrator);
    this->mgis_integrator->addBehaviourIntegrator(n, m, l, b);
  }  // end of addBehaviourIntegrator

  const Material& NonLinearEvolutionProblemImplementationBase::getMaterial(
      const size_type m) const {
    checkMultiMaterialSupportEnabled("getMaterial", this->mgis_integrator);
    return this->mgis_integrator->getMaterial(m);
  }  // end of getMaterial

  Material& NonLinearEvolutionProblemImplementationBase::getMaterial(
      const size_type m) {
    checkMultiMaterialSupportEnabled("getMaterial", this->mgis_integrator);
    return this->mgis_integrator->getMaterial(m);
  }  // end of getMaterial

  const BehaviourIntegrator&
  NonLinearEvolutionProblemImplementationBase::getBehaviourIntegrator(
      const size_type m) const {
    checkMultiMaterialSupportEnabled("getBehaviourIntegrator",
                                     this->mgis_integrator);
    return this->mgis_integrator->getBehaviourIntegrator(m);
  }  // end of getBehaviourIntegrator

  BehaviourIntegrator&
  NonLinearEvolutionProblemImplementationBase::getBehaviourIntegrator(
      const size_type m) {
    checkMultiMaterialSupportEnabled("getBehaviourIntegrator",
                                     this->mgis_integrator);
    return this->mgis_integrator->getBehaviourIntegrator(m);
  }  // end of getBehaviourIntegrator

  void NonLinearEvolutionProblemImplementationBase::setTimeIncrement(
      const real dt) {
    if (this->mgis_integrator != nullptr) {
      this->mgis_integrator->setTimeIncrement(dt);
    }
  }  // end of setTimeIncrement

  void NonLinearEvolutionProblemImplementationBase::setMacroscopicGradients(
      const std::vector<real>& g) {
    checkMultiMaterialSupportEnabled("setMacroscopicGradients",
                                     this->mgis_integrator);
    this->mgis_integrator->setMacroscopicGradients(g);
  }  // end of setMacroscopicGradients

  void NonLinearEvolutionProblemImplementationBase::setSolverParameters(
      const Parameters& params) {
    mfem_mgis::setSolverParameters(*(this->solver), params);
  }  // end of setSolverParameters

  void NonLinearEvolutionProblemImplementationBase::setup(const real t,
                                                          const real dt) {
    if (this->initialization_phase) {
      if (!this->dirichlet_boundary_conditions.empty()) {
        auto fixed_dirichlet_dofs = std::vector<mfem_mgis::size_type>{};
        for (const auto& bc : this->dirichlet_boundary_conditions) {
          auto dofs = bc->getHandledDegreesOfFreedom();
          fixed_dirichlet_dofs.insert(fixed_dirichlet_dofs.end(), dofs.begin(),
                                      dofs.end());
        }
        this->markDegreesOfFreedomHandledByDirichletBoundaryConditions(
            fixed_dirichlet_dofs);
      }
    }
    this->initialization_phase = false;
    for (const auto& bc : this->dirichlet_boundary_conditions) {
      bc->updateImposedValues(this->u1, t + dt);
    }
    if (this->mgis_integrator != nullptr) {
      this->mgis_integrator->setup(t, dt);
    }
  }  // end of setup

  void NonLinearEvolutionProblemImplementationBase::updateLinearSolver(
      std::unique_ptr<LinearSolver> s) {
    this->linear_solver_preconditioner.reset();
    this->linear_solver = std::move(s);
    this->solver->setLinearSolver(*(this->linear_solver));
  }  // end of updateLinearSolver

  void NonLinearEvolutionProblemImplementationBase::updateLinearSolver(
      std::unique_ptr<LinearSolver> s,
      std::unique_ptr<LinearSolverPreconditioner> p) {
    if (p != nullptr) {
      auto* const isolver = dynamic_cast<IterativeSolver*>(s.get());
      if (isolver == nullptr) {
        mgis::raise(
            "NonLinearEvolutionProblemImplementationBase::updateLinearSolver: "
            "can't associate a preconditioner to a non iterative solver");
      }
      isolver->SetPreconditioner(*p);
      this->updateLinearSolver(std::move(s));
      this->linear_solver_preconditioner = std::move(p);
    } else {
      this->updateLinearSolver(std::move(s));
    }
  }  // end of updateLinearSolver

  void NonLinearEvolutionProblemImplementationBase::updateLinearSolver(
      LinearSolverHandler s) {
    this->updateLinearSolver(std::move(s.linear_solver),
                             std::move(s.preconditioner));
  }  // end of updateLinearSolver

  void NonLinearEvolutionProblemImplementationBase::addBoundaryCondition(
      std::unique_ptr<DirichletBoundaryCondition> bc) {
    this->dirichlet_boundary_conditions.push_back(std::move(bc));
  }  // end of addBoundaryCondition

  bool NonLinearEvolutionProblemImplementationBase::solve(const real t,
                                                          const real dt) {
    mfem::Vector zero;
    this->setTimeIncrement(dt);
    this->setup(t, dt);
    //    this->computePrediction(t, dt);
    this->solver->Mult(zero, this->u1);
    return this->solver->GetConverged();
  }  // end of solve

  void NonLinearEvolutionProblemImplementationBase::computePrediction(
      const real, const real) {
    //    constexpr auto it = IntegrationType::PREDICTION_ELASTIC_OPERATOR;
    //     std::cout << "ComputePrediction\n";
    //     mfem::Vector c;
    //     mfem::Vector r;
    //     r.SetSize(this->u1.Size());
    //     c.SetSize(this->u1.Size());
    //     if (!this->integrate(this->u1, it)) {
    //       return;
    //     }
    //     std::cout << "this->dirichlet_boundary_conditions: "
    //               << this->dirichlet_boundary_conditions.size() << '\n';
    //     for (const auto& bc : this->dirichlet_boundary_conditions) {
    //       std::cerr << "calling bc: " << bc.get() << '\n';
    //       bc->setImposedValuesIncrements(c, t, t + dt);
    //       bc->setImposedValuesIncrements(this->u1, t, t + dt);
    //     }
    //     this->solver->computeResidual(r, this->u1);
    //     if (!this->solver->computeNewtonCorrection(c, r, this->u1)) {
    //       std::cerr << "Marche pas !\n";
    //       return;
    //     }
    //     c.Print(std::cout);
    //     std::cout << '\n';
    //     this->u1 -= c;
    //     for (const auto& bc : this->dirichlet_boundary_conditions) {
    //       bc->updateImposedValues(this->u1, t + dt);
    //     }
  }  // end of computePrediction

  NonLinearEvolutionProblemImplementationBase::
      ~NonLinearEvolutionProblemImplementationBase() = default;

}  // end of namespace mfem_mgis
