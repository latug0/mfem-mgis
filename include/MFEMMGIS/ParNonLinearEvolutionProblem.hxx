/*!
 * \file   include/MFEMMGIS/NonLinearEvolutionProblem.hxx
 * \brief
 * \author Guillaume Latu
 * \date 02/02/2021
 */

#ifndef LIB_MFEM_MGIS_PAREVOLUTIONPROBLEM_HXX
#define LIB_MFEM_MGIS_PAREVOLUTIONPROBLEM_HXX

#include <memory>
#include "mfem/linalg/vector.hpp"
#include "mfem/linalg/solvers.hpp"
#include "mfem/fem/pnonlinearform.hpp"
#include "MGIS/Behaviour/Hypothesis.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/ParFiniteElementDiscretization.hxx"

#ifdef MFEM_USE_MPI

namespace mfem_mgis {

  // forward declaration
  struct Material;
  // forward declaration
  struct MultiMaterialNonLinearIntegrator;

  /*!
   * \brief class for solving non linear evolution problems
   */
  struct MFEM_MGIS_EXPORT ParNonLinearEvolutionProblem : public mfem::ParNonlinearForm {
    //! \brief a simple alias
    using Hypothesis = mgis::behaviour::Hypothesis;
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     * \param[in] h: modelling hypothesis
     */
    ParNonLinearEvolutionProblem(std::shared_ptr<ParFiniteElementDiscretization>,
                              const Hypothesis);
    //! \return the finite element space
    mfem::ParFiniteElementSpace& getFiniteElementSpace();
    //! \return the finite element space
    const mfem::ParFiniteElementSpace& getFiniteElementSpace() const;
    //! \return the Newton solver
    mfem::NewtonSolver& getSolver();
    //! \return the unknowns at the beginning of the time step
    mfem::Vector& getUnknownsAtBeginningOfTheTimeStep();
    //! \return the unknowns at the beginning of the time step
    const mfem::Vector& getUnknownsAtBeginningOfTheTimeStep() const;
    //! \return the unknowns at the end of the time step
    mfem::Vector& getUnknownsAtEndOfTheTimeStep();
    //! \return the unknowns at the end of the time step
    const mfem::Vector& getUnknownsAtEndOfTheTimeStep() const;
    /*!
     * \return the material with the given id
     * \param[in] m: material id
     */
    const Material& getMaterial(const size_type) const;
    /*!
     * \return the material with the given id
     * \param[in] m: material id
     */
    Material& getMaterial(const size_type);
    /*!
     * \brief add a new material
     * \param[in] n: name of the behaviour integrator
     * \param[in] m: material id
     * \param[in] l: library name
     * \param[in] b: behaviour name
     */
    virtual void addBehaviourIntegrator(const std::string&,
                                        const size_type,
                                        const std::string&,
                                        const std::string&);
    /*!
     * \brief solve the non linear problem over the given time step
     * \param[in] dt: time increment
     */
    void solve(const real);
    /*!
     * \brief revert the internal state variables to the beginning of the time
     * step.
     */
    void revert();
    /*!
     * \brief updat the internal state variables to the end of the time step.
     */
    void update();

    //! \brief destructor
    ~ParNonLinearEvolutionProblem() override;

   private:
    //! \brief pointer to the underlying domain integrator
    MultiMaterialNonLinearIntegrator* const mgis_integrator;
    //! \brief underlying finite element discretization
    const std::shared_ptr<ParFiniteElementDiscretization> fe_discretization;
    //! \brief modelling hypothesis
    const Hypothesis hypothesis;
    //! \brief newton solver
    mfem::NewtonSolver solver;
    //! \brief unknowns at the beginning of the time step
    mfem::Vector u0;
    //! \brief unknowns at the end of the time step
    mfem::Vector u1;
  };  // end of struct ParNonLinearEvolutionProblem

}  // end of namespace mfem_mgis

#endif /* MFEM_USE_MPI */

#endif /* LIB_MFEM_MGIS_PAREVOLUTIONPROBLEM */
