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
#include "MFEMMGIS/AnalyticalTests.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"

#ifndef MFEM_USE_MPI
#define MPI_COMM_WORLD 0
#define MPI_Finalize(args...) \
  {}
#define MPI_Init(args...) \
  {}
#endif

void (*getSolution(const std::size_t i))(mfem::Vector&, const mfem::Vector&) {
  constexpr const auto xthr = mfem_mgis::real(1) / 2;
  std::array<void (*)(mfem::Vector&, const mfem::Vector&), 6u> solutions = {
      +[](mfem::Vector& u, const mfem::Vector& x) {
        constexpr const auto gradx = mfem_mgis::real(1) / 3;
        u = mfem_mgis::real{};
        if (x(0) < xthr) {
          u(0) = gradx * x(0);
        } else {
          u(0) = gradx * xthr - gradx * (x(0) - xthr);
        }
      },
      +[](mfem::Vector& u, const mfem::Vector& x) {
        constexpr const auto gradx = mfem_mgis::real(4) / 30;
        u = mfem_mgis::real{};
        if (x(0) < xthr) {
          u(0) = gradx * x(0);
        } else {
          u(0) = gradx * xthr - gradx * (x(0) - xthr);
        }
      },
      +[](mfem::Vector& u, const mfem::Vector& x) {
        constexpr const auto gradx = mfem_mgis::real(4) / 30;
        u = mfem_mgis::real{};
        if (x(0) < xthr) {
          u(0) = gradx * x(0);
        } else {
          u(0) = gradx * xthr - gradx * (x(0) - xthr);
        }
      },
      +[](mfem::Vector& u, const mfem::Vector& x) {
        constexpr const auto gradx = mfem_mgis::real(1) / 3;
        u = mfem_mgis::real{};
        if (x(0) < xthr) {
          u(1) = gradx * x(0);
        } else {
          u(1) = gradx * xthr - gradx * (x(0) - xthr);
        }
      },
      +[](mfem::Vector& u, const mfem::Vector& x) {
        constexpr const auto gradx = mfem_mgis::real(1) / 3;
        u = mfem_mgis::real{};
        if (x(0) < xthr) {
          u(2) = gradx * x(0);
        } else {
          u(2) = gradx * xthr - gradx * (x(0) - xthr);
        }
      },
      +[](mfem::Vector& u, const mfem::Vector&) { u = mfem_mgis::real{}; }};
  return solutions[i];
}

static void setLinearSolver(mfem_mgis::AbstractNonLinearEvolutionProblem& p,
                            const std::size_t i) {
  if (i == 0) {
    p.setLinearSolver("GMRESSolver", {{"VerbosityLevel", 1},
                                      {"AbsoluteTolerance", 1e-12},
                                      {"RelativeTolerance", 1e-12},
                                      {"MaximumNumberOfIterations", 300}});
  } else if (i == 1) {
    p.setLinearSolver("CGSolver", {{"VerbosityLevel", 1},
                                   {"AbsoluteTolerance", 1e-12},
                                   {"RelativeTolerance", 1e-12},
                                   {"MaximumNumberOfIterations", 300}});
#ifdef MFEM_USE_SUITESPARSE
  } else if (i == 2) {
    p.setLinearSolver("UMFPackSolver", {});
#endif
#ifdef MFEM_USE_MUMPS
  } else if (i == 3) {
    p.setLinearSolver("MUMPSSolver", {{"Type", "UNSYMMETRIC"}});
#endif
  } else {
    std::cerr << "unsupported linear solver\n";
    std::exit(EXIT_FAILURE);
  }
}

static void setSolverParameters(
    mfem_mgis::AbstractNonLinearEvolutionProblem& problem) {
  problem.setSolverParameters({{"VerbosityLevel", 0},
                               {"RelativeTolerance", 1e-12},
                               {"AbsoluteTolerance", 1e-12},
                               {"MaximumNumberOfIterations", 10}});
}  // end of setSolverParmeters

bool checkSolution(mfem_mgis::NonLinearEvolutionProblem& problem,
                   const std::size_t i) {
  const auto b = mfem_mgis::compareToAnalyticalSolution(
      problem, getSolution(i), {{"CriterionThreshold", 1e-10}});
  if (!b) {
    std::cerr << "Error is greater than threshold\n";
    return false;
  }
  std::cerr << "Error is lower than threshold\n";
  return true;
}

struct TestParameters {
  const char* mesh_file = nullptr;
  const char* library = nullptr;
  int order = 1;
  int tcase = 0;
  int linearsolver = 0;
};

TestParameters parseCommandLineOptions(int& argc, char* argv[]) {
  TestParameters p;
  // options treatment
  MPI_Init(&argc, &argv);
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&p.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&p.library, "-l", "--library", "Material library.");
  args.AddOption(&p.order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&p.tcase, "-t", "--test-case",
                 "identifier of the case : Exx->0, Eyy->1, Ezz->2, Exy->3, "
                 "Exz->4, Eyz->5");
  args.AddOption(
      &p.linearsolver, "-ls", "--linearsolver",
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

void exit_on_failure() {
  MPI_Finalize();
  std::exit(EXIT_FAILURE);
}

void executeMFEMMGISTest(const TestParameters& p) {
  constexpr const auto dim = mfem_mgis::size_type{3};
  // creating the finite element workspace

  std::shared_ptr<mfem_mgis::FiniteElementDiscretization> fed;
  auto smesh = std::make_shared<mfem::Mesh>(p.mesh_file, 1, 1);
  if (dim != smesh->Dimension()) {
    std::cerr << "Invalid mesh dimension \n";
    exit_on_failure();
  }

  //   auto fed = FiniteElementDiscretization({{"MeshFileName", p.mesh_file},
  //                                           {"FiniteElementCollection",
  //                                           "H1"},
  //                                           {"FiniteElementOrder", p.order},
  //                                           {"parallel", true}});

#ifdef DO_USE_MPI
  auto mesh = std::make_shared<mfem::ParMesh>(MPI_COMM_WORLD, *smesh);
  mesh->UniformRefinement();
  mesh->UniformRefinement();
#else
  auto mesh = smesh;
#endif
  fed = std::make_shared<mfem_mgis::FiniteElementDiscretization>(
      mesh, std::make_shared<mfem::H1_FECollection>(p.order, dim), dim);

  {
    // building the non linear problem
    mfem_mgis::PeriodicNonLinearEvolutionProblem problem(fed);
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
    setLinearSolver(problem, p.linearsolver);
    setSolverParameters(problem);
    //
    problem.addPostProcessing(
        "ParaviewExportResults",
        {{"OutputFileName", "PeriodicTestOutput-" + std::to_string(p.tcase)}});
    // solving the problem
    problem.solve(0, 1);
    problem.executePostProcessings(0, 1);
    //
    if (!checkSolution(problem, p.tcase)) {
      exit_on_failure();
    }
  }
  MPI_Finalize();
}

int main(int argc, char* argv[]) {
  const auto p = parseCommandLineOptions(argc, argv);
  executeMFEMMGISTest(p);
  return EXIT_SUCCESS;
}
