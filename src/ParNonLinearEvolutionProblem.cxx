/*!
 * \file   src/ParNonLinearEvolutionProblem.cxx
 * \brief
 * \author GUillaume Latu
 * \date   03/02/2021
 */

#include "MGIS/Raise.hxx"
#include "MFEMMGIS/MultiMaterialNonLinearIntegrator.hxx"
#include "MFEMMGIS/ParNonLinearEvolutionProblem.hxx"

namespace mfem_mgis {

  ParNonLinearEvolutionProblem::ParNonLinearEvolutionProblem(
      std::shared_ptr<ParFiniteElementDiscretization> fed, const Hypothesis h)
      : mfem::ParNonlinearForm(&(fed->getFiniteElementSpace())),
    mgis_integrator(new MultiMaterialNonLinearIntegrator( // TOCHECK
        fed->getFiniteElementDiscretization(), h)),
        fe_discretization(fed),
        hypothesis(h),
        u0(fed->getFiniteElementSpace().GetTrueVSize()),
        u1(fed->getFiniteElementSpace().GetTrueVSize()) {
    this->u0 = real{0};
    this->u1 = real{0};
    this->solver.SetOperator(*(this));
    this->solver.iterative_mode = true;
    if (this->fe_discretization->getMesh().Dimension() !=
        mgis::behaviour::getSpaceDimension(h)) {
      mgis::raise(
          "NonLinearEvolutionProblem::NonLinearEvolutionProblem: "
          "modelling hypothesis is not consistent with the spatial dimension "
          "of the mesh");
    }
    std::cout << __FILE__ << " LINE: " <<__LINE__ << std::endl;;
    this->AddDomainIntegrator(this->mgis_integrator);
    std::cout << __FILE__ << " LINE: " <<__LINE__ << std::endl;;
  }  // end of NonLinearEvolutionProblem

  mfem::Vector&
  ParNonLinearEvolutionProblem::getUnknownsAtBeginningOfTheTimeStep() {
    return this->u0;
  }  // end of getUnknownsAtBeginningOfTheTimeStep

  const mfem::Vector&
  ParNonLinearEvolutionProblem::getUnknownsAtBeginningOfTheTimeStep() const {
    return this->u0;
  }  // end of getUnknownsAtBeginningOfTheTimeStep

  mfem::Vector& ParNonLinearEvolutionProblem::getUnknownsAtEndOfTheTimeStep() {
    return this->u1;
  }  // end of getUnknownsAtEndOfTheTimeStep

  const mfem::Vector& ParNonLinearEvolutionProblem::getUnknownsAtEndOfTheTimeStep()
      const {
    return this->u1;
  }  // end of getUnknownsAtEndOfTheTimeStep

  const mfem::ParFiniteElementSpace&
  ParNonLinearEvolutionProblem::getFiniteElementSpace() const {
    return this->fe_discretization->getFiniteElementSpace();
  }  // end of NonLinearEvolutionProblem::getFiniteElementSpace

  mfem::ParFiniteElementSpace& ParNonLinearEvolutionProblem::getFiniteElementSpace() {
    return this->fe_discretization->getFiniteElementSpace();
  }  // end of NonLinearEvolutionProblem::getFiniteElementSpace

  mfem::NewtonSolver& ParNonLinearEvolutionProblem::getSolver() {
    return this->solver;
  }  // end of NonLinearEvolutionProblem

  void ParNonLinearEvolutionProblem::revert() {
    this->u1 = this->u0;
    this->mgis_integrator->revert();
  }  // end of revert

  void ParNonLinearEvolutionProblem::update() {
    this->u0 = this->u1;
    this->mgis_integrator->update();
  }  // end of update

  void ParNonLinearEvolutionProblem::addBehaviourIntegrator(const std::string& n,
							    const size_type m,
							    const std::string& l,
							    const std::string& b) {
    this->mgis_integrator->addBehaviourIntegrator(n, m, l, b);
  }  // end of addBehaviourIntegrator

  const Material& ParNonLinearEvolutionProblem::getMaterial(
      const size_type m) const {
    return this->mgis_integrator->getMaterial(m);
  }  // end of ParNonLinearEvolutionProblem::getMaterial

  Material& ParNonLinearEvolutionProblem::getMaterial(const size_type m) {
    return this->mgis_integrator->getMaterial(m);
  }  // end of ParNonLinearEvolutionProblem::getMaterial

  void ParNonLinearEvolutionProblem::solve(const real dt) {
    mfem::Vector zero;
    this->mgis_integrator->setTimeIncrement(dt);
    this->solver.Mult(zero, this->u1);
    if (!this->solver.GetConverged()) {
      mgis::raise("Newton solver did not converge");
    }
  }  // end of solve

  ParNonLinearEvolutionProblem::~ParNonLinearEvolutionProblem() = default;

}  // end of namespace mfem_mgis
