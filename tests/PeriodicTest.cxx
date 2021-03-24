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
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementationBase.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"

#ifndef MFEM_USE_MPI
#define MPI_COMM_WORLD 0
#define MPI_Finalize(args...) {}
#define MPI_Allreduce(args...) {}
#define MPI_Init(args...) {}
#endif

void (*getSolution(const std::size_t i))(const mfem::Vector&, mfem::Vector&) {
  constexpr const auto xthr = mfem_mgis::real(1) / 2;
  std::array<void (*)(const mfem::Vector&, mfem::Vector&), 6u> solutions = {
      +[](const mfem::Vector& x, mfem::Vector& u) {
        constexpr const auto gradx = mfem_mgis::real(1) / 3;
        u = mfem_mgis::real{};
        if (x(0) < xthr) {
          u(0) = gradx * x(0);
        } else {
          u(0) = gradx * xthr - gradx * (x(0) - xthr);
        }
      },
      +[](const mfem::Vector& x, mfem::Vector& u) {
        constexpr const auto gradx = mfem_mgis::real(4) / 30;
        u = mfem_mgis::real{};
        if (x(0) < xthr) {
          u(0) = gradx * x(0);
        } else {
          u(0) = gradx * xthr - gradx * (x(0) - xthr);
        }
      },
      +[](const mfem::Vector& x, mfem::Vector& u) {
        constexpr const auto gradx = mfem_mgis::real(4) / 30;
        u = mfem_mgis::real{};
        if (x(0) < xthr) {
          u(0) = gradx * x(0);
        } else {
          u(0) = gradx * xthr - gradx * (x(0) - xthr);
        }
      },
      +[](const mfem::Vector& x, mfem::Vector& u) {
        constexpr const auto gradx = mfem_mgis::real(1) / 3;
        u = mfem_mgis::real{};
        if (x(0) < xthr) {
          u(1) = gradx * x(0);
        } else {
          u(1) = gradx * xthr - gradx * (x(0) - xthr);
        }
      },
      +[](const mfem::Vector& x, mfem::Vector& u) {
        constexpr const auto gradx = mfem_mgis::real(1) / 3;
        u = mfem_mgis::real{};
        if (x(0) < xthr) {
          u(2) = gradx * x(0);
        } else {
          u(2) = gradx * xthr - gradx * (x(0) - xthr);
        }
      },
      +[](const mfem::Vector&, mfem::Vector& u) { u = mfem_mgis::real{}; }};
  return solutions[i];
}

template <bool parallel>
std::shared_ptr<mfem::Solver> getLinearSolver(const std::size_t i) {
  using generator = std::shared_ptr<mfem::Solver> (*)();
  const generator generators[] = {
      []() -> std::shared_ptr<mfem::Solver> {
        std::shared_ptr<mfem::GMRESSolver> pgmres(nullptr);
        if constexpr (parallel) 
          pgmres = std::make_shared<mfem::GMRESSolver>(MPI_COMM_WORLD);
        else
          pgmres = std::make_shared<mfem::GMRESSolver>();
        pgmres->iterative_mode = false;
        pgmres->SetRelTol(1e-12);
        pgmres->SetAbsTol(1e-12);
        pgmres->SetMaxIter(300);
        pgmres->SetPrintLevel(1);
        return pgmres;
      },
      []() -> std::shared_ptr<mfem::Solver> {
        std::shared_ptr<mfem::CGSolver> pcg(nullptr);                                
        if constexpr (parallel) 
          pcg = std::make_shared<mfem::CGSolver>(MPI_COMM_WORLD);
        else
          pcg = std::make_shared<mfem::CGSolver>();
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

template <bool parallel>
void setSolverParameters(
    mfem_mgis::NonLinearEvolutionProblemImplementationBase<parallel>& problem,
    mfem::Solver& lsolver) {
  auto& solver = problem.getSolver();
  solver.iterative_mode = true;
  solver.SetSolver(lsolver);
  solver.SetPrintLevel(0);
  solver.SetRelTol(1e-12);
  solver.SetAbsTol(1e-12);
  solver.SetMaxIter(10);
}  // end of setSolverParmeters

template <bool parallel>
bool checkSolution(mfem_mgis::NonLinearEvolutionProblemImplementationBase<parallel>& problem,
                   const std::size_t i) {
  constexpr const auto eps = mfem_mgis::real{1e-10};
  const auto dim = problem.getFiniteElementSpace().GetMesh()->Dimension();
  // recover the solution as a grid function
  auto& u1 = problem.getUnknownsAtEndOfTheTimeStep();
  mfem_mgis::GridFunction<parallel> x(&problem.getFiniteElementSpace());
  x.MakeTRef(&problem.getFiniteElementSpace(), u1, 0);
  x.SetFromTrueVector();
  // comparison to analytical solution
  mfem::VectorFunctionCoefficient sol_coef(dim, getSolution(i));
  const auto error = x.ComputeL2Error(sol_coef);
  if (error > eps) {
    std::cerr << "Error is greater than threshold (" << error << " > " << eps << ")\n";
    return false;
  }
  std::cerr << "Error is lower than threshold (" << error << " < " << eps
            << ")\n";
  return true;
}

template <bool parallel>
void exportResults(mfem_mgis::NonLinearEvolutionProblemImplementationBase<parallel>& problem,
                   const std::size_t tcase) {
  auto* const mesh = problem.getFiniteElementSpace().GetMesh();
  auto& u1 = problem.getUnknownsAtEndOfTheTimeStep();
  mfem_mgis::GridFunction<parallel> x(&problem.getFiniteElementSpace());
  x.MakeTRef(&problem.getFiniteElementSpace(), u1, 0);
  x.SetFromTrueVector();

  mfem::ParaViewDataCollection paraview_dc(
      "PeriodicTestOutput-" + std::to_string(tcase), mesh);
  paraview_dc.RegisterField("u", &x);
  paraview_dc.SetCycle(0);
  paraview_dc.SetTime(0.0);
  paraview_dc.Save();
};

struct TestParameters {
  const char* mesh_file = nullptr;
  const char* library = nullptr;
  int order = 1;
  int tcase = 0;
  int linearsolver = 0;
};

template <bool parallel>
TestParameters parseCommandLineOptions(int &argc, char* argv[]){
  TestParameters p;
  // options treatment
  if constexpr (parallel) MPI_Init(&argc, &argv);
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
  if ((!args.Good()) || (p.mesh_file == nullptr)) {
    args.PrintUsage(std::cout);
    std::exit(EXIT_FAILURE);
  }
  args.PrintOptions(std::cout);
  if ((p.tcase < 0) || (p.tcase > 5)) {
    std::cerr << "Invalid test case\n";
    std::exit(EXIT_FAILURE);
  }
  return p;
}

template <bool parallel>
void exit_on_failure () {
  if constexpr (parallel) MPI_Finalize();
  std::exit(EXIT_FAILURE);
}

template <bool parallel>
void executeMFEMMGISTest(const TestParameters& p) {
  constexpr const auto dim = mfem_mgis::size_type{3};
  // creating the finite element workspace

  std::shared_ptr<mfem_mgis::Mesh<parallel>> mesh;
  if constexpr (parallel) {
      auto smesh = std::make_shared<mfem::Mesh>(p.mesh_file, 1, 1);
      if (dim != smesh->Dimension()) {
        std::cerr << "Invalid mesh dimension \n";
        exit_on_failure<parallel>();
      }
      for (int i = 0 ; i < 2 ; i++)
        smesh->UniformRefinement();
      mesh = std::make_shared<mfem::ParMesh>(MPI_COMM_WORLD, *smesh);
    } else {
      mesh = std::make_shared<mfem::Mesh>(p.mesh_file, 1, 1);
      if (dim != mesh->Dimension()) {
        std::cerr << "Invalid mesh dimension \n";
        exit_on_failure<parallel>();
      }
  }
  // building the non linear problem
  mfem_mgis::PeriodicNonLinearEvolutionProblem problem(
      std::make_shared<mfem_mgis::FiniteElementDiscretization>(
          mesh, std::make_shared<mfem::H1_FECollection>(p.order, dim), 3));
  problem.addBehaviourIntegrator("Mechanics", 1, p.library, "Elasticity");
  problem.addBehaviourIntegrator("Mechanics", 2, p.library, "Elasticity");
  // materials
  auto& m1 = problem.getMaterial(1);
  auto& m2 = problem.getMaterial(2);
  // setting the material properties
  auto set_properties = [](auto& m, const double l, const double mu) {
    mgis::behaviour::setMaterialProperty(m.s0, "FirstLameCoefficient", l);
    mgis::behaviour::setMaterialProperty(m.s0, "ShearModulus", mu);
    mgis::behaviour::setMaterialProperty(m.s1, "FirstLameCoefficient", l);
    mgis::behaviour::setMaterialProperty(m.s1, "ShearModulus", mu);
  };
  set_properties(m1, 100, 75);
  set_properties(m2, 200, 150);
  //
  auto set_temperature = [](auto& m) {
    mgis::behaviour::setExternalStateVariable(m.s0, "Temperature", 293.15);
    mgis::behaviour::setExternalStateVariable(m.s1, "Temperature", 293.15);
  };
  set_temperature(m1);
  set_temperature(m2);
  // macroscopic strain
  std::vector<mfem_mgis::real> e(6, mfem_mgis::real{});
  if (p.tcase < 3) {
    e[p.tcase] = 1;
  } else {
    e[p.tcase] = 1.41421356237309504880 / 2;
  }
  problem.setMacroscopicGradientsEvolution([e](const double) { return e; });
  //
  auto lsolver = getLinearSolver<parallel>(p.linearsolver);
  setSolverParameters(problem, *(lsolver.get()));
  // solving the problem
  problem.solve(0, 1);
  //
  if (!checkSolution(problem, p.tcase)) {
    exit_on_failure<parallel>();
  }
  //
  exportResults(problem, p.tcase);
  if constexpr(parallel) MPI_Finalize();
}

struct ElasticityNonLinearIntegrator final
    : public mfem::NonlinearFormIntegrator {

  void AssembleElementVector(const mfem::FiniteElement&,
                             mfem::ElementTransformation&,
                             const mfem::Vector&,
                             mfem::Vector&) override {}

  void AssembleElementGrad(const mfem::FiniteElement&,
                           mfem::ElementTransformation&,
                           const mfem::Vector&,
                           mfem::DenseMatrix&) override {}

  void setMacroscopicGradients(const std::vector<mfem_mgis::real>& g) {
    this->e = g;
  }

 private:
  //! macroscopic gradients
  std::vector<mfem_mgis::real> e;
};

void executeMFEMMTest(const TestParameters& p) {
  constexpr const auto dim = mfem_mgis::size_type{3};
  // creating the finite element workspace
  auto mesh = std::make_shared<mfem::Mesh>(p.mesh_file, 1, 1);
  if (dim != mesh->Dimension()) {
    std::cerr << "Invalid mesh dimension \n";
    std::exit(EXIT_FAILURE);
  }
  // building the non linear problem
  mfem_mgis::NonLinearEvolutionProblemImplementationBase<false> problem(
      std::make_shared<mfem_mgis::FiniteElementDiscretization>(
          mesh, std::make_shared<mfem::H1_FECollection>(p.order, dim), 3));
  std::vector<mfem_mgis::real> e(6, mfem_mgis::real{});
  if (p.tcase < 3) {
    e[p.tcase] = 1;
  } else {
    e[p.tcase] = 1 / 2;
  }
  //
  problem.AddDomainIntegrator(new ElasticityNonLinearIntegrator());
  //
  mfem_mgis::setPeriodicBoundaryConditions(problem);
  //
  auto lsolver = getLinearSolver<false>(p.linearsolver);
  setSolverParameters(problem, *(lsolver.get()));
  // solving the problem
  problem.solve(0, 1);
  //
  if (!checkSolution(problem, p.tcase)) {
    std::exit(EXIT_FAILURE);
  }
  //
  exportResults(problem, p.tcase);
}

int main(int argc, char* argv[]) {
#ifdef DO_USE_MPI  
  const auto p = parseCommandLineOptions<true>(argc, argv);
  executeMFEMMGISTest<true>(p);
#else
  const auto p = parseCommandLineOptions<false>(argc, argv);
  executeMFEMMGISTest<false>(p);
#endif
  return EXIT_SUCCESS;
}
