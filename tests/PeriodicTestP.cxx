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
 * \date   03/02/2021
 */

#include <memory>
#include <cstdlib>
#include <iostream>
#include "mfem.hpp"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

int main(int argc, char* argv[]) {
   // 1. Initialize MPI.
  int num_procs, myid;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  std::cerr << "Starting " << myid << std::endl;
  std::cerr << "LINE: " <<__LINE__ << std::endl;;
  constexpr const auto dim = mfem_mgis::size_type{3};
  constexpr const auto eps = mfem_mgis::real{1e-10};
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
  std::shared_ptr<mfem::Solver> (*const solver_list[])() {
    []() ->  std::shared_ptr<mfem::Solver> {
        std::shared_ptr<mfem::CGSolver> pcg(new mfem::CGSolver(MPI_COMM_WORLD));
        pcg->SetRelTol(1e-13);
        pcg->SetMaxIter(300);
        pcg->SetPrintLevel(1);
        return pcg;
    }
  };

  const char* mesh_file = nullptr;
  const char* library = nullptr;
  int order = 1;
  int tcase = 0;
  int linearsolver = 1;
  // options treatment
  mfem::OptionsParser args(argc, argv);

  MPI_Barrier(MPI_COMM_WORLD);
  std::cerr << __FILE__ << " LINE: " <<__LINE__ << std::endl;

  args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&library, "-l", "--library", "Material library.");
  args.AddOption(&order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&tcase, "-t", "--test-case",
                 "identifier of the case : Exx->0, Eyy->1, Ezz->2, Exy->3, "
                 "Exz->4, Eyz->5");
  args.AddOption(&linearsolver, "-ls", "--linearsolver",
                 "identifier of the linear solver: 0 -> GMRES, 1 -> CG, 2 -> UMFPack");
  args.Parse();
  MPI_Barrier(MPI_COMM_WORLD);
  std::cerr << __FILE__ << " LINE: " <<__LINE__ << std::endl;
  if (!args.Good() || (mesh_file == nullptr))
    {
      if (myid == 0)
	{
	  args.PrintUsage(std::cout);
	}
      MPI_Finalize();
      return EXIT_FAILURE;
    }
  if (myid == 0)
    {
      args.PrintOptions(std::cout);
    }
  
  if ((tcase < 0) || (tcase > 5)) {
    if (myid == 0)
      {
	std::cerr << "Invalid test case\n";
      }
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  // creating the finite element workspace
  mfem::Mesh* mesh = new mfem::Mesh(mesh_file, 1, 1);
  if (dim != mesh->Dimension()) {
    std::cerr << "Invalid mesh dimension \n";
    return EXIT_FAILURE;
  }
  for (int i = 0 ; i < 3 ; i++)
    mesh->UniformRefinement();

  mfem::ParMesh *pmesh = new mfem::ParMesh(MPI_COMM_WORLD, *mesh);
  delete mesh;
  MPI_Barrier(MPI_COMM_WORLD);
  std::shared_ptr<mfem::ParMesh> spmesh(pmesh);
  auto pfespace = std::make_shared<mfem_mgis::FiniteElementDiscretization>(
      spmesh, std::make_shared<mfem::H1_FECollection>(order, dim), 3);

  std::cerr << __FILE__ << " LINE: " <<__LINE__ << std::endl;;
  // building the non linear problem
  mfem_mgis::NonLinearEvolutionProblem problem(pfespace,
	       mgis::behaviour::Hypothesis::TRIDIMENSIONAL);

  problem.addBehaviourIntegrator("Mechanics", 1, library, "Elasticity");
  problem.addBehaviourIntegrator("Mechanics", 2, library, "Elasticity");
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
  if (tcase < 3) {
    e[tcase] = 1;
  } else {
    e[tcase] = 1.41421356237309504880 / 2;
  }
  m1.setMacroscopicGradients(e);
  m2.setMacroscopicGradients(e);
  // Impose no displacement on the first node
  // which needs to be on x=xmin or x=xmax axis.
  // ux=0, uy=0, uz=0 on this point.
  const auto nnodes = problem.getFiniteElementSpace().GetTrueVSize() / dim;
  std::cerr << "Number of nodes: " << nnodes << std::endl;
  mfem::Array<int> ess_tdof_list;
  ess_tdof_list.SetSize(0);
  {
    mfem::GridFunction nodes(&problem.getFiniteElementSpace());
    int found = 0;
    int size = nnodes;
    bool reorder_space = true;
    pmesh->GetNodes(nodes);

     // Traversal of all dofs to detect which one is (0,0,0)
     for (int i = 0; i < size; ++i) {
       double coord[dim]; // coordinates of a node
       double dist = 0.;
       for (int j = 0; j < dim; ++j) {
	 if (reorder_space)
	   coord[j] = (nodes)[j * size + i];
	 else
	   coord[j] = (nodes)[i * dim + j];
	 // because of periodic BC, 0. is also 1.
	 if (abs(coord[j] - 1.) < 1e-7)
	   coord[j] = 0.;
	 dist += coord[j] * coord[j];
       }
       // If distance is close to zero, we have our reference point
       if (dist < 1.e-16) {
	 //	 cout << myid << "coord: " <<coord[0] << " " << coord[1] << endl;
	 for (int j = 0; j < dim; ++j) {
	   int id_unk;
	   if (reorder_space)
	     {
	     //id_unk = (j * size + i);
	     id_unk = problem.getFiniteElementSpace().GetLocalTDofNumber(j * size + i);
	     }
	   else
	     {
	       //id_unk = (i * dim + j);
	       id_unk = problem.getFiniteElementSpace().GetLocalTDofNumber(i * dim + j);
	     }
	   if (id_unk >= 0)
	     {
	       found = 1;
	       ess_tdof_list.Append(id_unk);
	       std::cout << myid << "bloqued unknown: " << id_unk << std::endl;
	     }
	 }
       }
     }
     MPI_Allreduce(MPI_IN_PLACE, &found, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
     MFEM_VERIFY(found, "Reference point at (0,0) was not found");
  }
  problem.SetEssentialTrueDofs(ess_tdof_list);
  // solving the problem
  std::shared_ptr<mfem::Solver> lsolver = solver_list[linearsolver]();

  auto& solver = problem.getSolver();
  solver.iterative_mode = true;
  solver.SetSolver(*lsolver);
  solver.SetPrintLevel(0);
  solver.SetRelTol(1e-12);
  solver.SetAbsTol(1e-12);
  solver.SetMaxIter(10);
  problem.solve(1);
  // recover the solution as a grid function
  auto& u1 = problem.getUnknownsAtEndOfTheTimeStep();
  mfem_mgis::mGridFunction x(&problem.getFiniteElementSpace());
  x.MakeTRef(&problem.getFiniteElementSpace(), u1, 0);
  x.SetFromTrueVector();
  // comparison to analytical solution
  mfem::VectorFunctionCoefficient sol_coef(dim, solutions[tcase]);
  const auto error = x.ComputeL2Error(sol_coef);
  if (error > eps) {
    std::cerr << "Error is greater than threshold (" << error << " > " << eps << ")\n";
    return EXIT_FAILURE;
  } else {
    std::cerr << "Error is lower than threshold (" << error << " < " << eps << ")\n";
  }

  // exporting the results
  mfem::ParaViewDataCollection paraview_dc(
      "PeriodicTestOutput-" + std::to_string(tcase), pmesh);
  paraview_dc.RegisterField("u", &x);
  paraview_dc.SetCycle(0);
  paraview_dc.SetTime(0.0);
  paraview_dc.Save();

  MPI_Finalize();
  return EXIT_SUCCESS;
}
