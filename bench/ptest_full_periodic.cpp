//                                      
// Compile with: make ptest_full_periodic
//
// Author: Guillaume Latu
//
// Description:  This example code solves a simple linear elasticity problem
//               describing a multi-material square.
//               This problem has a 1D analytic solution along x1 dimension,
//               the solution is constant along x2 dimension which is also periodic.
//
//               The geometry of the domain is assumed to be as
//               follows:
//
//                           x2=1  +----------+----------+
//                                 | material | material |      
//                                 |    1     |    2     |      
//                           x2=0  +----------+----------+      
//                               x1=0                   x1=1
//
//
//               Specifically, we approximate the weak form of -div(sigma(u))=0
//               where sigma(u)=lambda*div(u)*I+mu*(grad*u+u*grad) is the stress
//               tensor corresponding to displacement field u, and lambda and mu
//               are the material Lame constants. The boundary conditions are
//               periodic.
//
//               Mechanical strain:
//                               eps = E + grad_s v
//
//                         with  E the given macrocoscopic strain
//                               v the periodic displacement fluctuation 
//               Displacement:
//                                 u = U + v
//
//                         with  U the given displacement associated to E
//                                 E = grad_s U
//
//               The local microscopic strain is equal, on average, to the macroscopic strain: 
//                         <eps> = <E>
//                    
//               Equation to be solved:   div(C grad_s v)= - div(C grad_s U)
//
//               Thus, we introduce the force term -div(C grad_s U) in addition to the 
//               classical elasticity problem.
//

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <chrono>
#include <ctime>

#define USE_PROFILER 1
#define LIB_PROFILER_PRINTF MpiPrintf
#include "libProfiler.h"
using namespace std;
using namespace mfem;

// central frontier in-between the two materials
const double xthr = 0.5;

/// Class for domain integrator L(v) := (f, grad v)
class DomainLFGrad2Integrator : public LinearFormIntegrator
{
protected:
  Coefficient *lambda, *mu;
private:
  Vector shape, Qvec, Svec;
  DenseMatrix dshape;
  int tcase;

public:
   /// Constructs the domain integrator (Q, grad v)
  DomainLFGrad2Integrator(Coefficient &l, Coefficient &m, int _tcase)
  {
    lambda = &l; mu = &m; tcase = _tcase;
  }

  void AssembleRHSElementVect(
    const FiniteElement &el, ElementTransformation &Tr, Vector &elvect)
  {
    int dof = el.GetDof();
    int spaceDim = Tr.GetSpaceDim();
    
    dshape.SetSize(dof, spaceDim);
    Qvec.SetSize(3);
    Svec.SetSize(3);
    elvect.SetSize(dof*spaceDim);
    elvect = 0.0;
    
    const IntegrationRule *ir = IntRule;
    if (ir == NULL)
      {
	int intorder = 2 * el.GetOrder();
	ir = &IntRules.Get(el.GetGeomType(), intorder);
      }
    
    for (int i = 0; i < ir->GetNPoints(); i++)
      {
	const IntegrationPoint &ip = ir->IntPoint(i);
	Tr.SetIntPoint(&ip);
	el.CalcPhysDShape(Tr, dshape);

	double M = mu->Eval(Tr, ip);
	double L = lambda->Eval(Tr, ip);
	switch(tcase) {
	case 1:
	  Svec[0]=0;Qvec[0] = -(L+2*M);
	  Svec[1]=1;Qvec[1] = -L;
	  Svec[2]=2;Qvec[2] = -L;
	  break;
	case 2:
	  Svec[0]=0;Qvec[0] = -L;
	  Svec[1]=1;Qvec[1] = -(L+2*M);
	  Svec[2]=2;Qvec[2] = -L;
	  break;
	case 3:
	  Svec[0]=0;Qvec[0] = -L;
	  Svec[1]=1;Qvec[1] = -L;
	  Svec[2]=2;Qvec[2] = -(L+2*M);
	  break;
	case 4:
	  Svec[0]=1;Qvec[0] = -M;
	  Svec[1]=0;Qvec[1] = -M;
	  Svec[2]=2;Qvec[2] = 0;
	  break;
	case 5:
	  Svec[0]=0;Qvec[0] = 0;
	  Svec[1]=2;Qvec[1] = -M;
	  Svec[2]=1;Qvec[2] = -M;
	  break;
	case 6:
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
  }
  
  
};

void sol_exact1(const Vector &x, Vector &u);
void sol_exact2(const Vector &x, Vector &u);
void sol_exact3(const Vector &x, Vector &u);
void sol_exact4(const Vector &x, Vector &u);
void sol_exact5(const Vector &x, Vector &u);
void sol_exact6(const Vector &x, Vector &u);

int main(int argc, char *argv[])
{
   PROFILER_ENABLE;
   // Initialize MPI.
   int num_procs, myid;
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);

   PROFILER_START(0_total);
   PROFILER_START(1_initialize);

   // Parse command-line options.
   //   const char *mesh_file = "square_2mat_per.msh";
   const char *mesh_file = "cube_2mat_per.mesh";
   int order = 1;
   int tcase = 1;
   bool visualization = 1;
   bool reorder_space = true;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.AddOption(&tcase, "-t", "--tcase",
                  "identifier of the case : Exx->1, Eyy->2, Ezz->3, Exy->4, Eyz->5, Exz->6");
   args.AddOption(&reorder_space, "-nodes", "--by-nodes", "-vdim", "--by-vdim",
                  "Use byNODES ordering of vector space instead of byVDIM");
   args.Parse();
   if (!args.Good())
   {
      if (myid == 0)
      {
         args.PrintUsage(cout);
      }
      MPI_Finalize();
      return 1;
   }
   if (myid == 0)
   {
      args.PrintOptions(cout);
   }

   PROFILER_END(); PROFILER_START(2_read_mesh);
   //   Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral or hexahedral elements with the same code.
   Mesh *mesh = new Mesh(mesh_file, 1, 0);
   int dim = mesh->Dimension();

   if ( dim == 2 && tcase != 1 && tcase != 2 && tcase != 4 ) {
     if (myid == 0)
       cout << "This test case is undefined in 2D" << endl;
     MPI_Finalize();
     exit(1);
   }
   if ( dim == 3 && ((tcase < 1) || (tcase > 6))) {
     if (myid == 0)
       cout << "This test case is undefined in 3D" << endl;
     MPI_Finalize();
     exit(1);
   }

   PROFILER_END(); PROFILER_START(3_refine_mesh);
   PROFILER_START(3.1_sequential);
   //  Refine the mesh.
   {
      int ref_levels = 0;
      for (int l = 0; l < ref_levels; l++)
      {
         mesh->UniformRefinement();
      }
   }
   
   PROFILER_END(); PROFILER_START(3.2_parallel);
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);

   //  Refine the mesh.
   {
      int par_ref_levels = 4;
      cout << "par_ref_levels " << par_ref_levels << endl;
      for (int l = 0; l < par_ref_levels; l++)
      {
         pmesh->UniformRefinement();
      }
   }
   pmesh->RemoveUnusedVertices();
   delete mesh;

  PROFILER_END(); 
  PROFILER_END(); 
  PROFILER_START(4_initialize_fem); 
  PROFILER_START(4.1_define_system_components);
   //    Define a finite element space on the mesh. Here we use vector finite
   //    elements, i.e. dim copies of a scalar finite element space. The vector
   //    dimension is specified by the last argument of the FiniteElementSpace
   //    constructor. 
   FiniteElementCollection *fec;
   ParFiniteElementSpace *fespace;
   fec = new H1_FECollection(order, dim);
   if (reorder_space)
     {
       fespace = new ParFiniteElementSpace(pmesh, fec, dim, Ordering::byNODES);
     }
   else
     {
       fespace = new ParFiniteElementSpace(pmesh, fec, dim, Ordering::byVDIM);
     }
   
   cout << myid << ": Number of finite element local unknowns: " << fespace->GetTrueVSize() << endl;
   auto gsize = fespace->GlobalVSize();
   if (myid == 0) 
     cout << "Nb global unknowns: "<< gsize << endl;

   //    Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.
   ParGridFunction x(fespace);
   x = 0.; 

   // Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the linear elasticity integrator with piece-wise
   //    constants coefficient lambda and mu.
   Vector lambda(pmesh->attributes.Max());
   lambda = 100.0;
   if (pmesh->attributes.Max() > 1)
     lambda(1) = lambda(0)*2;
   PWConstCoefficient lambda_func(lambda); 
   Vector mu(pmesh->attributes.Max());
   mu = 75.0;
   if (pmesh->attributes.Max() > 1)
     mu(1) = mu(0)*2;
   PWConstCoefficient mu_func(mu);

   // Impose no displacement on the first node
   // which needs to be on x=xmin or x=xmax axis.
   // ux=0, uy=0, uz=0 on this point.
   Array<int> ess_tdof_list;
   ess_tdof_list.SetSize(0);
   {
     GridFunction nodes(fespace);
     int found = 0;
     pmesh->GetNodes(nodes);
     const auto size = nodes.Size()/dim;
     std::cerr << "Number of nodes: " << size << std::endl;

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
	     id_unk = fespace->GetLocalTDofNumber(j * size + i);
	     }
	   else
	     {
	       //id_unk = (i * dim + j);
	       id_unk = fespace->GetLocalTDofNumber(i * dim + j);
	     }
	   if (id_unk >= 0)
	     {
	       found = 1;
	       ess_tdof_list.Append(id_unk);
	       x(id_unk) = 0.;
	       cout << myid << "bloqued unknown: " << id_unk << endl;
	     }
	 }
       }
     }
     MPI_Allreduce(MPI_IN_PLACE, &found, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
     MFEM_VERIFY(found, "Reference point at (0,0) was not found");
   }
   ParBilinearForm *a = new ParBilinearForm(fespace);
   a->AddDomainIntegrator(new ElasticityIntegrator(lambda_func,mu_func));

   // Set up the right-hand side of the FEM linear system.
   ParLinearForm *b = new ParLinearForm(fespace);
   b->AddDomainIntegrator(new DomainLFGrad2Integrator(lambda_func,mu_func,tcase));

   PROFILER_END(); PROFILER_START(4.2_assembly);

   a->Assemble();
   b->Assemble();

   OperatorPtr A;
   Vector B, X;
   a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);
   if (myid == 0)
   {
      cout << "done." << endl;
   }

   PROFILER_END();
   PROFILER_END(); 
   PROFILER_START(5_solve_linear_system);
   CGSolver *pcg = new CGSolver(MPI_COMM_WORLD);
   pcg->SetRelTol(1e-9);
   pcg->SetAbsTol(1e-9);
   pcg->SetMaxIter(1200);
   pcg->SetPrintLevel(5);
   pcg->SetOperator(*A);
   pcg->Mult(B, X);

   // 14. Recover the parallel grid function corresponding to X. This is the
   //     local finite element solution on each processor.
   a->RecoverFEMSolution(X, *b, x);


   PROFILER_END(); 
   PROFILER_START(6_postprocess);
   //  Save the result, the displacement u
   {
     ParaViewDataCollection paraview_dc("per", pmesh);
     paraview_dc.RegisterField("u",&x);
     paraview_dc.SetCycle(0);
     paraview_dc.SetTime(0.0);
     paraview_dc.Save();
   }

   void (*sol_exact)(const Vector &x, Vector &u);
   switch(tcase) {
   case 1:
     sol_exact=sol_exact1;
     break;
   case 2:
     sol_exact=sol_exact2;
     break;
   case 3:
     sol_exact=sol_exact3;
     break;
   case 4:
     sol_exact=sol_exact4;
     break;
   case 5:
     sol_exact=sol_exact5;
     break;
   case 6:
     sol_exact=sol_exact6;
     break;
   default:
     cout << "Test undefined" << endl;
     exit (1);
     break;
   }
   VectorFunctionCoefficient sol_coef (dim, sol_exact);
   double errorL2 = x.ComputeL2Error(sol_coef);
   PROFILER_END(); 
   PROFILER_END(); 
   if (myid == 0) 
     LogProfiler();
   
   PROFILER_DISABLE;
   if (errorL2 < 1e-9)
   {
     if (myid == 0)
     {
       cerr<<"\ntcase " << tcase << " -- L2 norm: " << errorL2 << endl;
       cerr << "OK" << endl;
     }
     MPI_Finalize();
     return 0;
   }
   else
   {
     if (myid == 0)
     {
       cerr<<"\ntcase " << tcase << " -- L2 norm: " << errorL2 << endl;
       cerr << "Fail" << endl;
     }
     MPI_Finalize();
     return 1;
   }
}



void sol_exact1(const Vector &x, Vector &u)
{
  const double gradx = 1./3.;
  u =  0.;
  if (x(0) < xthr)
    {
      u(0) = gradx*x(0);
    }
  else
    {
      u(0) = gradx*xthr - gradx*(x(0)-xthr);
    }
}

void sol_exact2(const Vector &x, Vector &u)
{
  const double gradx = 4./30.;
  u =  0.;
  if (x(0) < xthr)
    {
      u(0) = gradx*x(0);
    }
  else
    {
      u(0) = gradx*xthr - gradx*(x(0)-xthr);
    }
}

void sol_exact3(const Vector &x, Vector &u)
{
  const double gradx = 4./30.;
  u =  0.;
  if (x(0) < xthr)
    {
      u(0) = gradx*x(0);
    }
  else
    {
      u(0) = gradx*xthr - gradx*(x(0)-xthr);
    }
}

void sol_exact4(const Vector &x, Vector &u)
{
  const double gradx = 1./3.;
  u =  0.;
  if (x(0) < xthr)
    {
      u(1) = gradx*x(0);
    }
  else
    {
      u(1) = gradx*xthr - gradx*(x(0)-xthr);
    }
}

void sol_exact5(const Vector &x, Vector &u)
{
  u =  0.;
}

void sol_exact6(const Vector &x, Vector &u)
{
  const double gradx = 1./3.;
  u =  0.;
  if (x(0) < xthr)
    {
      u(2) = gradx*x(0);
    }
  else
    {
      u(2) = gradx*xthr - gradx*(x(0)-xthr);
    }
}
