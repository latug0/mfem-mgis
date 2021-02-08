/*!
 * \file   tests/PeriodicTest.cxx
 * \brief
 * This example code solves a simple linear elasticity problem
 * describing a multi-material square.
 * This problem has a 1D analytic solution along x1 dimension,
 * the solution is constant along x2 dimension which is also periodic.
 *
 * The geometry of the domain is assumed to be as
 * follows:
 *
 *             x2=1  +----------+----------+
 *                   | material | material |
 *                   |    1     |    2     |
 *             x2=0  +----------+----------+
 *                 x1=0                   x1=1
 *
 *
 * Specifically, we approximate the weak form of -div(sigma(u))=0
 * where sigma(u)=lambda*div(u)*I+mu*(grad*u+u*grad) is the stress
 * tensor corresponding to displacement field u, and lambda and mu
 * are the material Lame constants. The boundary conditions are
 * periodic.
 *
 * Mechanical strain:
 *                 eps = E + grad_s v
 *
 *           with  E the given macrocoscopic strain
 *                 v the periodic displacement fluctuation
 * Displacement:
 *                   u = U + v
 *
 *           with  U the given displacement associated to E
 *                   E = grad_s U
 *
 * The local microscopic strain is equal, on average, to the macroscopic strain:
 *           <eps> = <E>
 * \author Thomas Helfer, Guillaume Latu
 * \date   14/10/2020
 */

#include <memory>
#include <cstdlib>
#include <iostream>
#include "mfem/general/optparser.hpp"
#include "mfem/linalg/solvers.hpp"
#include "mfem/fem/datacollection.hpp"
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/BehaviourIntegrator.hxx"
#include "MFEMMGIS/BehaviourIntegratorFactory.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

void (*getSolution(const std::size_t i))(const mfem::Vector&, mfem::Vector&)
{
    constexpr const auto xthr = mfem_mgis::real(1) / 2;
    std::array<void (*)(const mfem::Vector&, mfem::Vector&), 6u> solutions =
    {
        +[](const mfem::Vector& x, mfem::Vector& u)
        {
            constexpr const auto gradx = mfem_mgis::real(1) / 3;
            u = mfem_mgis::real{};
            if (x(0) < xthr)
            {
                u(0) = gradx * x(0);
            }
            else
            {
                u(0) = gradx * xthr - gradx * (x(0) - xthr);
            }
        },
        +[](const mfem::Vector& x, mfem::Vector& u)
        {
            constexpr const auto gradx = mfem_mgis::real(4) / 30;
            u = mfem_mgis::real{};
            if (x(0) < xthr)
            {
                u(0) = gradx * x(0);
            }
            else
            {
                u(0) = gradx * xthr - gradx * (x(0) - xthr);
            }
        },
        +[](const mfem::Vector& x, mfem::Vector& u)
        {
            constexpr const auto gradx = mfem_mgis::real(4) / 30;
            u = mfem_mgis::real{};
            if (x(0) < xthr)
            {
                u(0) = gradx * x(0);
            }
            else
            {
                u(0) = gradx * xthr - gradx * (x(0) - xthr);
            }
        },
        +[](const mfem::Vector& x, mfem::Vector& u)
        {
            constexpr const auto gradx = mfem_mgis::real(1) / 3;
            u = mfem_mgis::real{};
            if (x(0) < xthr)
            {
                u(1) = gradx * x(0);
            }
            else
            {
                u(1) = gradx * xthr - gradx * (x(0) - xthr);
            }
        },
        +[](const mfem::Vector& x, mfem::Vector& u)
        {
            constexpr const auto gradx = mfem_mgis::real(1) / 3;
            u = mfem_mgis::real{};
            if (x(0) < xthr)
            {
                u(2) = gradx * x(0);
            }
            else
            {
                u(2) = gradx * xthr - gradx * (x(0) - xthr);
            }
        },
        +[](const mfem::Vector&, mfem::Vector& u)
        {
            u = mfem_mgis::real{};
        }
    };
    return solutions[i];
}

std::shared_ptr<mfem::Solver> getLinearSolver(const std::size_t i)
{
    using generator = std::shared_ptr<mfem::Solver> (*)();
    const generator generators[] =
    {
        []() -> std::shared_ptr<mfem::Solver> {
            std::shared_ptr<mfem::GMRESSolver> pgmres(new mfem::GMRESSolver);
            pgmres->iterative_mode = false;
            pgmres->SetRelTol(1e-12);
            pgmres->SetAbsTol(1e-12);
            pgmres->SetMaxIter(300);
            pgmres->SetPrintLevel(1);
            return pgmres;
        },
        []() -> std::shared_ptr<mfem::Solver> {
            std::shared_ptr<mfem::CGSolver> pcg(new mfem::CGSolver);
            pcg->SetRelTol(1e-12);
            pcg->SetMaxIter(300);
            pcg->SetPrintLevel(1);
            return pcg;
        }
#ifdef MFEM_USE_SUITESPARSE
        ,
        []() -> std::shared_ptr<mfem::Solver> {
            std::shared_ptr<mfem::UMFPackSolver> pumf(new mfem::UMFPackSolver);
            return (pumf);
        }
#endif
    };
    return (generators[i])();
}

void setBoundaryConditions(mfem_mgis::NonLinearEvolutionProblemBase& problem)
{
    // Impose no displacement on the first node
    // which needs to be on x=xmin or x=xmax axis.
    // ux=0, uy=0, uz=0 on this point.
    const auto dim = problem.getFiniteElementSpace().GetMesh()->Dimension();
    const auto nnodes = problem.getFiniteElementSpace().GetTrueVSize() / dim;
    mfem::Array<int> ess_tdof_list;
    ess_tdof_list.SetSize(dim);
    for (int k = 0; k < dim; k++)
    {
        int tgdof = 0 + k * nnodes;
        ess_tdof_list[k] = tgdof;
    }
    problem.SetEssentialTrueDofs(ess_tdof_list);
}

void setSolverParameters(mfem_mgis::NonLinearEvolutionProblemBase& problem,
                         mfem::Solver& lsolver)
{
    auto& solver = problem.getSolver();
    solver.iterative_mode = true;
    solver.SetSolver(lsolver);
    solver.SetPrintLevel(0);
    solver.SetRelTol(1e-12);
    solver.SetAbsTol(1e-12);
    solver.SetMaxIter(10);
}  // end of setSolverParameters

bool checkSolution(mfem_mgis::NonLinearEvolutionProblemBase& problem,
                   const std::size_t i)
{
    constexpr const auto eps = mfem_mgis::real{1e-10};
    const auto dim = problem.getFiniteElementSpace().GetMesh()->Dimension();
    // recover the solution as a grid function
    auto& u1 = problem.getUnknownsAtEndOfTheTimeStep();
    mfem::GridFunction x(&problem.getFiniteElementSpace());
    x.MakeTRef(&problem.getFiniteElementSpace(), u1, 0);
    x.SetFromTrueVector();
    // comparison to analytical solution
    mfem::VectorFunctionCoefficient sol_coef(dim, getSolution(i));
    const auto error = x.ComputeL2Error(sol_coef);
    if (error > eps)
    {
        std::cerr << "Error is greater than threshold (" << error << " > " << eps << ")\n";
        return false;
    }
    std::cerr << "Error is lower than threshold (" << error << " < " << eps
              << ")\n";
    return true;
}

void exportResults(mfem_mgis::NonLinearEvolutionProblemBase& problem,
                   const std::size_t tcase)
{
    auto* const mesh = problem.getFiniteElementSpace().GetMesh();
    auto& u1 = problem.getUnknownsAtEndOfTheTimeStep();
    mfem::GridFunction x(&problem.getFiniteElementSpace());
    x.MakeTRef(&problem.getFiniteElementSpace(), u1, 0);
    x.SetFromTrueVector();

    mfem::ParaViewDataCollection paraview_dc(
        "PeriodicTestOutput-" + std::to_string(tcase), mesh);
    paraview_dc.RegisterField("u", &x);
    paraview_dc.SetCycle(0);
    paraview_dc.SetTime(0.0);
    paraview_dc.Save();
};

struct TestParameters
{
    const char* mesh_file = nullptr;
    const char* library = nullptr;
    int order = 1;
    int tcase = 0;
    int linearsolver = 0;
};

TestParameters parseCommandLineOptions(const int argc, char** const argv)
{
    TestParameters p;
    // options treatment
    mfem::OptionsParser args(argc, argv);
    args.AddOption(&p.mesh_file, "-m", "--mesh", "Mesh file to use.");
    args.AddOption(&p.library, "-l", "--library", "Material library.");
    args.AddOption(&p.order, "-o", "--order",
                   "Finite element order (polynomial degree).");
    args.AddOption(&p.tcase, "-t", "--test-case",
                   "identifier of the case : Exx->0, Eyy->1, Ezz->2, Exy->3, "
                   "Exz->4, Eyz->5");
    args.AddOption(&p.linearsolver, "-ls", "--linearsolver",
                   "identifier of the linear solver: 0 -> GMRES, 1 -> CG, 2 -> UMFPack");
    args.Parse();
    if ((!args.Good()) || (p.mesh_file == nullptr))
    {
        args.PrintUsage(std::cout);
        std::exit(EXIT_FAILURE);
    }
    args.PrintOptions(std::cout);
    if ((p.tcase < 0) || (p.tcase > 5))
    {
        std::cerr << "Invalid test case\n";
        std::exit(EXIT_FAILURE);
    }
    return p;
}

void executeMFEMMGISTest(const TestParameters& p)
{
    constexpr const auto dim = mfem_mgis::size_type{3};
    // creating the finite element workspace
    auto mesh = std::make_shared<mfem::Mesh>(p.mesh_file, 1, 1);
    if (dim != mesh->Dimension())
    {
        std::cerr << "Invalid mesh dimension \n";
        std::exit(EXIT_FAILURE);
    }
    // building the non linear problem
    mfem_mgis::NonLinearEvolutionProblem problem(
        std::make_shared<mfem_mgis::FiniteElementDiscretization>(
            mesh, std::make_shared<mfem::H1_FECollection>(p.order, dim), 3),
        mgis::behaviour::Hypothesis::TRIDIMENSIONAL);
    problem.addBehaviourIntegrator("Mechanics", 1, p.library, "Elasticity");
    problem.addBehaviourIntegrator("Mechanics", 2, p.library, "Elasticity");
    // materials
    auto& m1 = problem.getMaterial(1);
    auto& m2 = problem.getMaterial(2);
    // setting the material properties
    auto set_properties = [](auto& m, const double l, const double mu)
    {
        mgis::behaviour::setMaterialProperty(m.s0, "FirstLameCoefficient", l);
        mgis::behaviour::setMaterialProperty(m.s0, "ShearModulus", mu);
        mgis::behaviour::setMaterialProperty(m.s1, "FirstLameCoefficient", l);
        mgis::behaviour::setMaterialProperty(m.s1, "ShearModulus", mu);
    };
    set_properties(m1, 100, 75);
    set_properties(m2, 200, 150);
    //
    auto set_temperature = [](auto& m)
    {
        mgis::behaviour::setExternalStateVariable(m.s0, "Temperature", 293.15);
        mgis::behaviour::setExternalStateVariable(m.s1, "Temperature", 293.15);
    };
    set_temperature(m1);
    set_temperature(m2);
    // macroscopic strain
    std::vector<mfem_mgis::real> e(6, mfem_mgis::real{});
    if (p.tcase < 3)
    {
        e[p.tcase] = 1;
    }
    else
    {
        e[p.tcase] = 1.41421356237309504880 / 2;
    }
    m1.setMacroscopicGradients(e);
    m2.setMacroscopicGradients(e);
    //
    setBoundaryConditions(problem);
    //
    auto lsolver = getLinearSolver(p.linearsolver);
    setSolverParameters(problem, *(lsolver.get()));
    // solving the problem
    problem.solve(1);
    //
    if (!checkSolution(problem, p.tcase))
    {
        std::exit(EXIT_FAILURE);
    }
    //
    exportResults(problem, p.tcase);
}

namespace mfem_mgis
{
struct ElasticityNonLinearIntegrator final
    : public mfem::NonlinearFormIntegrator
{

    /*!
     * \brief a simple test
     * \param[in] i: pointer to behaviour integrator
     * \param[in] n: name of the calling method
     * \param[in] m: material id
     */
    static void checkIfBehaviourIntegratorIsDefined(
        const BehaviourIntegrator* const i,
        const char* const n,
        const size_type m)
    {
        if (i == nullptr)
        {
            mgis::raise("MultiMaterialNonLinearIntegrator::" + std::string(n) +
                        ": no behaviour integrator associated with material '"+
                        std::to_string(m) + "'");
        }
    }  // end if checkIfBehaviourIntegratorIsDefined

    ElasticityNonLinearIntegrator(std::shared_ptr<FiniteElementDiscretization> &fed, const Hypothesis h)
        : fe_discretization(fed),  hypothesis(h)
    {
        const auto& mesh = this->fe_discretization->getMesh();
        // shifting by one allows to directly use the material id to get
        // the behaviour integrators in all the other methods.
        // However, behaviour_integrators[0] will always be null.
        this->behaviour_integrators.resize(mesh.attributes.Max() + 1);
    }  // end of ElasticityNonLinearIntegrator

    /*!
     * \brief Compute part of the inner forces using a single element.
     *
     * this function is called automatically by NonlinearForm::Mult
     * to performe the assembly of the RHS.
     */
    void AssembleElementVector(const mfem::FiniteElement& e,
                               mfem::ElementTransformation& tr,
                               const mfem::Vector& U,
                               mfem::Vector& F) override
    {
        const auto m = tr.Attribute;
        const auto& bi = this->behaviour_integrators[m];
        checkIfBehaviourIntegratorIsDefined(bi.get(), "AssembleElementVector", m);
        bi->computeInnerForces(F, e, tr, U);
    }  // end of AssembleElementVector


    /*!
     * \brief Compute part of the stiffness matrix using a single element
     *
     * this function is called automatically by NonlinearForm::GetGradient
     * to perform the assembly of Stiffness matrix.
     */
    void AssembleElementGrad(const mfem::FiniteElement& e,
                             mfem::ElementTransformation& tr,
                             const mfem::Vector& U,
                             mfem::DenseMatrix& K) override
    {
        const auto m = tr.Attribute;
        const auto& bi = this->behaviour_integrators[m];
        checkIfBehaviourIntegratorIsDefined(bi.get(), "AssembleElementGrad", m);
        bi->computeStiffnessMatrix(K, e, tr, U);
    }  // end of AssembleElementGrade


    void addBehaviourIntegrator(const std::string& n,
                                const size_type m,
                                const std::string& l,
                                const std::string& b)
    {
        if (this->behaviour_integrators[m] != nullptr)
        {
            mgis::raise(
                "MultiMaterialNonLinearIntegrator::addBehaviourIntegrator: "
                "integrator already defined for material '" +
                std::to_string(m) + "'");
        }
        const auto& f = BehaviourIntegratorFactory::get(this->hypothesis);
        this->behaviour_integrators[m] =
            f.generate(n, *(this->fe_discretization), m,
                       mfem_mgis::load(l, b, this->hypothesis));
    }  // end of addBehaviourIntegrator

    const Material& getMaterial(
        const size_type m) const
    {
        const auto& bi = this->behaviour_integrators[m];
        checkIfBehaviourIntegratorIsDefined(bi.get(), "getMaterial", m);
        return bi->getMaterial();
    }  // end of getMaterial

    Material& getMaterial(const size_type m)
    {
        const auto& bi = this->behaviour_integrators[m];
        checkIfBehaviourIntegratorIsDefined(bi.get(), "getMaterial", m);
        return bi->getMaterial();
    }  // end of getMaterial


private:
    //! \brief Finite Element Discretization
    const std::shared_ptr<FiniteElementDiscretization> fe_discretization;
    //! \brief modelling hypothesis
    const Hypothesis hypothesis;
    /*!
      * \brief mapping between the material integrator and the behaviour
      * integrator.
      */
    std::vector<std::unique_ptr<BehaviourIntegrator>> behaviour_integrators;
};
}

void executeMFEMTest(const TestParameters& p)
{
    constexpr const auto dim = mfem_mgis::size_type{3};
    // creating the finite element workspace
    auto mesh = std::make_shared<mfem::Mesh>(p.mesh_file, 1, 1);
    if (dim != mesh->Dimension())
    {
        std::cerr << "Invalid mesh dimension \n";
        std::exit(EXIT_FAILURE);
    }
    // building the non linear problem
    auto fed = std::make_shared<mfem_mgis::FiniteElementDiscretization>(
                   mesh, std::make_shared<mfem::H1_FECollection>(p.order, dim),
                   3);
    mfem_mgis::NonLinearEvolutionProblemBase problem(fed);
    auto enli = new mfem_mgis::ElasticityNonLinearIntegrator(
        fed,mgis::behaviour::Hypothesis::TRIDIMENSIONAL);

    problem.AddDomainIntegrator(enli);
    enli->addBehaviourIntegrator("Mechanics", 1, p.library, "Elasticity");
    enli->addBehaviourIntegrator("Mechanics", 2, p.library, "Elasticity");
    // materials
    auto& m1 = enli->getMaterial(1);
    auto& m2 = enli->getMaterial(2);
    // setting the material properties
    auto set_properties = [](auto& m, const double l, const double mu)
			  {
			    mgis::behaviour::setMaterialProperty(m.s0, "FirstLameCoefficient", l);
			    mgis::behaviour::setMaterialProperty(m.s0, "ShearModulus", mu);
			    mgis::behaviour::setMaterialProperty(m.s1, "FirstLameCoefficient", l);
			    mgis::behaviour::setMaterialProperty(m.s1, "ShearModulus", mu);
			  };
    set_properties(m1, 100, 75);
    set_properties(m2, 200, 150);
    //
    auto set_temperature = [](auto& m)
			   {
			     mgis::behaviour::setExternalStateVariable(m.s0, "Temperature", 293.15);
			     mgis::behaviour::setExternalStateVariable(m.s1, "Temperature", 293.15);
			   };
    set_temperature(m1);
    set_temperature(m2);
    std::vector<mfem_mgis::real> e(6, mfem_mgis::real{});
    if (p.tcase < 3)
    {
        e[p.tcase] = 1;
    }
    else
    {
        e[p.tcase] = 1.41421356237309504880 / 2;
    }
    m1.setMacroscopicGradients(e);
    m2.setMacroscopicGradients(e);
    
    setBoundaryConditions(problem);
    //
    auto lsolver = getLinearSolver(p.linearsolver);
    setSolverParameters(problem, *(lsolver.get()));
    // solving the problem
    problem.solve(1);
    // check convergence
    if (!checkSolution(problem, p.tcase))
    {
        std::exit(EXIT_FAILURE);
    }
    // export to paraview 
    exportResults(problem, p.tcase);
}

int main(const int argc, char** const argv)
{
    const auto p = parseCommandLineOptions(argc, argv);
    //executeMFEMMGISTest(p);
    executeMFEMTest(p);
    return EXIT_SUCCESS;
}