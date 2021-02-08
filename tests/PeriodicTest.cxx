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
#include "mfem/fem/intrules.hpp"
#include "mfem/linalg/solvers.hpp"
#include "mfem/fem/datacollection.hpp"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

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

std::shared_ptr<mfem::Solver> getLinearSolver(const std::size_t i) {
  using generator = std::shared_ptr<mfem::Solver> (*)();
  const generator generators[] = {
      []() -> std::shared_ptr<mfem::Solver> {
        std::shared_ptr<mfem::GMRESSolver> pgmres(new mfem::GMRESSolver);
        pgmres->iterative_mode = true;
        pgmres->SetRelTol(1e-13);
        pgmres->SetAbsTol(1e-13);
        pgmres->SetMaxIter(300);
        pgmres->SetPrintLevel(1);
        return pgmres;
      },
      []() -> std::shared_ptr<mfem::Solver> {
        std::shared_ptr<mfem::CGSolver> pcg(new mfem::CGSolver);
        pcg->SetRelTol(1e-13);
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

void setBoundaryConditions(mfem_mgis::NonLinearEvolutionProblemBase& problem){
  // Impose no displacement on the first node
  // which needs to be on x=xmin or x=xmax axis.
  // ux=0, uy=0, uz=0 on this point.
  const auto dim = problem.getFiniteElementSpace().GetMesh()->Dimension();
  const auto nnodes = problem.getFiniteElementSpace().GetTrueVSize() / dim;
  mfem::Array<int> ess_tdof_list;
  ess_tdof_list.SetSize(dim);
  for (int k = 0; k < dim; k++) {
    int tgdof = 0 + k * nnodes;
    ess_tdof_list[k] = tgdof;
  }
  problem.SetEssentialTrueDofs(ess_tdof_list);
}

void setSolverParameters(mfem_mgis::NonLinearEvolutionProblemBase& problem,
                         mfem::Solver& lsolver) {
  auto& solver = problem.getSolver();
  solver.iterative_mode = false;
  solver.SetSolver(lsolver);
  solver.SetPrintLevel(0);
  solver.SetRelTol(1e-12);
  solver.SetAbsTol(1e-12);
  solver.SetMaxIter(10);
}  // end of setSolverParmeters

bool checkSolution(mfem_mgis::NonLinearEvolutionProblemBase& problem,
                   const std::size_t i) {
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
  if (error > eps) {
    std::cerr << "Error is greater than threshold (" << error << " > " << eps << ")\n";
    return false;
  }
  std::cerr << "Error is lower than threshold (" << error << " < " << eps
            << ")\n";
  return true;
}

void exportResults(mfem_mgis::NonLinearEvolutionProblemBase& problem,
                   const std::size_t tcase) {
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

struct TestParameters {
  const char* mesh_file = nullptr;
  const char* library = nullptr;
  int order = 1;
  int tcase = 0;
  int linearsolver = 0;
};

TestParameters parseCommandLineOptions(const int argc, char** const argv){
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

void executeMFEMMGISTest(const TestParameters& p) {
  constexpr const auto dim = mfem_mgis::size_type{3};
  // creating the finite element workspace
  auto mesh = std::make_shared<mfem::Mesh>(p.mesh_file, 1, 1);
  if (dim != mesh->Dimension()) {
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
  if (!checkSolution(problem, p.tcase)) {
    std::exit(EXIT_FAILURE);
  }
  //
  exportResults(problem, p.tcase);
}

struct ElasticityNonLinearIntegrator final
    : public mfem::NonlinearFormIntegrator {

  ElasticityNonLinearIntegrator(mfem::Coefficient &l,
				mfem::Coefficient &m,
				int _tcase) {
    lambda = &l; mu = &m; tcase = _tcase;
#ifdef MFEM_THREAD_SAFE
    MFEM_VERIFY(0,"MPI_THREAD_SAFE should be off");
#endif
  }
    /*!
     * \brief Compute part of the inner forces using a single element.
     *
     * this function is called automatically by NonlinearForm::Mult
     * to perform the assembly of the RHS.
     */
    void AssembleElementVector(const mfem::FiniteElement& el,
			       mfem::ElementTransformation& Tr,
			       const mfem::Vector& U,
			       mfem::Vector& elvect) override {
    int dof = el.GetDof();
    int spaceDim = Tr.GetSpaceDim();
    
    dshape.SetSize(dof, spaceDim);
    Qvec.SetSize(3);
    Svec.SetSize(3);
    elvect.SetSize(dof*spaceDim);
    elvect = 0.0;
    
    const mfem::IntegrationRule *ir =
      mfem::NonlinearFormIntegrator::IntRule;
    if (ir == NULL)
      {
	int intorder = 2 * el.GetOrder();
	ir = &mfem::IntRules.Get(el.GetGeomType(), intorder);
      }
    
    for (int i = 0; i < ir->GetNPoints(); i++)
      {
	const mfem::IntegrationPoint &ip = ir->IntPoint(i);
	Tr.SetIntPoint(&ip);
	el.CalcPhysDShape(Tr, dshape);

	double M = mu->Eval(Tr, ip);
	double L = lambda->Eval(Tr, ip);
	switch(tcase) {
	case 0:
	  Svec[0]=0;Qvec[0] = -(L+2*M);
	  Svec[1]=1;Qvec[1] = -L;
	  Svec[2]=2;Qvec[2] = -L;
	  break;
	case 1:
	  Svec[0]=0;Qvec[0] = -L;
	  Svec[1]=1;Qvec[1] = -(L+2*M);
	  Svec[2]=2;Qvec[2] = -L;
	  break;
	case 2:
	  Svec[0]=0;Qvec[0] = -L;
	  Svec[1]=1;Qvec[1] = -L;
	  Svec[2]=2;Qvec[2] = -(L+2*M);
	  break;
	case 3:
	  Svec[0]=1;Qvec[0] = -M;
	  Svec[1]=0;Qvec[1] = -M;
	  Svec[2]=2;Qvec[2] = 0;
	  break;
	case 4:
	  Svec[0]=0;Qvec[0] = 0;
	  Svec[1]=2;Qvec[1] = -M;
	  Svec[2]=1;Qvec[2] = -M;
	  break;
	case 5:
	  Svec[0]=2;Qvec[0] = -M;
	  Svec[1]=1;Qvec[1] = 0;
	  Svec[2]=0;Qvec[2] = -M;
	  break;
	}
	Qvec *= ip.weight * Tr.Weight();
	for (int k = 0; k < spaceDim; k++) 
	  {	
	    for (int s = 0; s < dof; s++)
	      {
		elvect(dof*k+s)	+= Qvec[k]*dshape(s,Svec[k]);
	      }
	  }
      }
    }  // end of AssembleElementVector

    
    /*!
     * \brief Compute part of the stiffness matrix using a single element
     *
     * this function is called automatically by NonlinearForm::GetGradient
     * to perform the assembly of Stiffness matrix.
     */
    void AssembleElementGrad(const mfem::FiniteElement& el,
			     mfem::ElementTransformation& Trans,
			     const mfem::Vector& U,
			     mfem::DenseMatrix& elmat) override {
    int dof = el.GetDof();
    int dim = el.GetDim();
    double w, L, M;
     
    MFEM_ASSERT(dim == Trans.GetSpaceDim(), "");
     
    cshape.SetSize(dof);
    gshape.SetSize(dof, dim);
    pelmat.SetSize(dof); // size dof*dof
    divshape.SetSize(dim*dof);
    elmat.SetSize(dof * dim);
    
    const mfem::IntegrationRule *ir =
      mfem::NonlinearFormIntegrator::IntRule;
    if (ir == NULL)
      {
	int order = 2 * Trans.OrderGrad(&el); // correct order?
	ir = &mfem::IntRules.Get(el.GetGeomType(), order);
      }
     
    elmat = 0.0;

    for (int nn = 0; nn < ir -> GetNPoints(); nn++)
      {
	const mfem::IntegrationPoint &ip = ir->IntPoint(nn);
	 
	Trans.SetIntPoint(&ip);
	// Each row of the result dshape contains
	// the derivatives of one shape function at the point ip.
	el.CalcPhysDShape(Trans, gshape);
	el.CalcPhysShape(Trans, cshape);
	// Get the transformation Trans for point ip
	// Get the weights associated to point ip
	w = ip.weight * Trans.Weight();
	// Multiply the derivatives by the inverse jacobian matrix
	// to get the derivatives along x, y and z
	// gshape contains these derivatives gshape(dof,dim)
	M = mu->Eval(Trans, ip);
	L = lambda->Eval(Trans, ip);

	for (int i = 0; i < dim; i++)
	  for (int j = 0; j < dim; j++)
	    {
	      for (int k = 0; k < dof; k++)
		for (int l = 0; l < dof; l++)
		  {
		    elmat(dof*i+k, dof*j+l) += 
		      (L * w) * gshape(k, i) * gshape(l, j) +
		      (M * w) * gshape(k, j) * gshape(l, i) ;
		  }
	      if (j == i) {
		for (int k = 0; k < dof; k++)
		  for (int l = 0; l < dof; l++)
		    {
		      double ftmp = 0;
		      for (int d = 0; d < dim; d++)
			ftmp += gshape(k, d) * gshape(l, d);
		      elmat(dof*i+k, dof*j+l) +=
			(M * w) * ftmp ;
		    }
	      }
	    }
      }
    }  // end of AssembleElementGrad
    

protected:
  // Coefficients related to materials
  mfem::Coefficient *lambda, *mu;
  int tcase;
private:
  // Several temporary buffers
  mfem::Vector shape, Qvec, Svec;
  mfem::DenseMatrix dshape;
  mfem::Vector cshape;
  mfem::DenseMatrix gshape, pelmat;
  mfem::Vector divshape;
};

void executeMFEMTest(const TestParameters& p) {
  constexpr const auto dim = mfem_mgis::size_type{3};
  // creating the finite element workspace
  auto mesh = std::make_shared<mfem::Mesh>(p.mesh_file, 1, 1);
  if (dim != mesh->Dimension()) {
    std::cerr << "Invalid mesh dimension \n";
    std::exit(EXIT_FAILURE);
  }
  // building the non linear problem
  mfem_mgis::NonLinearEvolutionProblemBase problem(
      std::make_shared<mfem_mgis::FiniteElementDiscretization>(
          mesh, std::make_shared<mfem::H1_FECollection>(p.order, dim), 3));
  mfem::Vector lambda(mesh->attributes.Max());
  lambda = 100.0;
  if (mesh->attributes.Max() > 1)
    lambda(1) = lambda(0)*2;
  mfem::PWConstCoefficient lambda_func(lambda); 
  // Question: pourquoi piece wise constant ?
  // lié aux attributs : constant par zone matériau	
  mfem::Vector mu(mesh->attributes.Max());
  mu = 75.0;
  if (mesh->attributes.Max() > 1)
    mu(1) = mu(0)*2;
  mfem::PWConstCoefficient mu_func(mu);
  problem.AddDomainIntegrator(new ElasticityNonLinearIntegrator(lambda_func,mu_func,p.tcase));
  //
  setBoundaryConditions(problem);
  //
  auto lsolver = getLinearSolver(p.linearsolver);
  setSolverParameters(problem, *(lsolver.get()));
  // solving the problem
  problem.solve(1);
  //
  if (!checkSolution(problem, p.tcase)) {
    std::exit(EXIT_FAILURE);
  }
  //
  exportResults(problem, p.tcase);
}

int main(const int argc, char** const argv) {
  const auto p = parseCommandLineOptions(argc, argv);
  //  executeMFEMMGISTest(p);
  executeMFEMTest(p);
  return EXIT_SUCCESS;
}
