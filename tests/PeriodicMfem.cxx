//                        Modified version (GL) of MFEM Example 2
//                                      
// Description:  This example code solves a simple linear elasticity problem
//               describing a multi-material square.
//               We want to reproduce the case described here :
//                https://comet-fenics.readthedocs.io/en/latest/demo/periodic_homog_elas/periodic_homog_elas.html
//               but with a simpler setting : case laminé 2d of Marc Josien
//                                            with only a 1D analytic solution along x1 dimension
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
//               Boundary conditions -> not well defined for the moment (Dirichlet

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;
#define USE_PROFILER 1
#define LIB_PROFILER_IMPLEMENTATION
#include "libProfiler.h"

int tcase = 0;

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
  }
  
  
};
///// Class for domain integrator L(v) := (f, grad v)
//class DomainLFGrad2Integrator : public LinearFormIntegrator
//{
//protected:
//  Coefficient *lambda, *mu;
//private:
//  Vector epsM, myvect, myres;
//  DenseMatrix tBmat, gtBmat, Emat, dshape;
//  int width, spaceDim;
//
//public:
//   /// Constructs the domain integrator (Q, grad v)
//  DomainLFGrad2Integrator(Coefficient &l, Coefficient &m, int sdim, int mcase)
//  {
//    lambda = &l; mu = &m;
//    spaceDim = sdim;
//    MFEM_ASSERT((spaceDim == 3)|| (spaceDim == 2), "space dimension has to be 2 or 3");
//    if (spaceDim == 2)
//      width = 3;
//    if (spaceDim == 3)
//      width = 6;
//    myvect.SetSize(width);
//    myres.SetSize(spaceDim);
//    epsM.SetSize(width);
//    epsM = 0.;
//    Emat.SetSize(width);
//    switch(mcase) {
//    case 1:
//      epsM(0)=1.;
//      break;
//    case 2:
//      epsM(1)=1.;
//      break;
//    case 3:
//      epsM(2)=1.;
//      break;
//    case 4:
//      epsM(3)=1.;
//      break;
//    case 5:
//      epsM(4)=1.;
//      break;
//    case 6:
//      epsM(5)=1.;
//      break;
//    }
//  }
//  
//  void AssembleRHSElementVect(
//    const FiniteElement &el, ElementTransformation &Tr, Vector &elvect) override 
//  {
//    int dof = el.GetDof();
//    gtBmat.SetSize(spaceDim,dof*width);
//    tBmat.SetSize(spaceDim,width);
//    dshape.SetSize(dof,spaceDim);
//    elvect.SetSize(dof*spaceDim);
//    elvect = 0.0;
//    MFEM_ASSERT(spaceDim == Tr.GetSpaceDim(), "");
//    
//    const IntegrationRule *ir = IntRule;
//    if (ir == NULL)
//      {
//	int intorder = 2 * el.GetOrder();
//	ir = &IntRules.Get(el.GetGeomType(), intorder);
//      }
//    
//    for (int i = 0; i < ir->GetNPoints(); i++)
//      {
//	const IntegrationPoint &ip = ir->IntPoint(i);
//	Tr.SetIntPoint(&ip);
//	el.CalcPhysDShape(Tr, dshape);
//	double M = mu->Eval(Tr, ip);
//	double L = lambda->Eval(Tr, ip);
//	// Fill Emat
//	Emat = 0.;
//	for (int k = 0; k< spaceDim; k++)
//	  {
//	    // Extra-diagonal terms Emat
//	    for (int l = 0; l< spaceDim; l++)
//	      Emat(k,l) = L;
//	    // Diagonal terms Emat
//	    Emat(k,k) = L+2*M;
//	  }
//	// Diagonal terms
//	for (int k = spaceDim; k < width; k++)
//	  Emat(k,k)= M;
//
//	// Fill gtBmat
//	gtBmat = 0.;
//	for (int s = 0; s < dof; s++)
//	  {
//	    int offset = width*s;
//	    // Fill gtBmat that depends on dof 's'.
//	    // gtBmat(:,offset:offset+width) contains tBmat
//	    // for dof 's'
//	    for (int k = 0; k < spaceDim; k++) 
//	      gtBmat(k,offset+k) = dshape(s,k);
//	    gtBmat(0,offset+spaceDim) = dshape(s,1); // N_y 
//	    gtBmat(1,offset+spaceDim) = dshape(s,0); // N_x 
//	    if (spaceDim == 3) {
//	      gtBmat(1,offset+4) = dshape(s,2); // N_z
//	      gtBmat(2,offset+4) = dshape(s,1); // N_y
//	      gtBmat(0,offset+5) = dshape(s,2); // N_z
//	      gtBmat(2,offset+5) = dshape(s,0); // N_x
//	    }
//	  }
//
//	for (int s = 0; s < dof; s++)
//	  {
//	    // Fill tBmat that depends on 's'
//	    int offset_s = width*s;
//	    tBmat.CopyCols(gtBmat,offset_s,offset_s+width-1);
//	    myvect = 0.;
//	    // Do the dense algebra
//	    Emat.Mult(epsM,myvect);
//	    tBmat.Mult(myvect,myres);
//	    // Multiply by the weight
//	    myres *= -1 * ip.weight * Tr.Weight();
//	    // Add the contribution to the RHS
//	    for (int k = 0; k < spaceDim; k++)
//	      elvect(dof*k+s) += myres[k];
//	  }
//      }
//
//}
//};



class Elasticity2Integrator: public BilinearFormIntegrator
{
protected:
  Coefficient *lambda, *mu;

private:
#ifndef MFEM_THREAD_SAFE
  DenseMatrix dshape, gtBmat, Emat;
  DenseMatrix BmatR, tBmatL, tBmatR, EBmat, tBEBmat;
#endif

public:
  Elasticity2Integrator(Coefficient &l, Coefficient &m)
  {
    lambda = &l; mu = &m;
  }

  /** Given a particular Finite Element
      computes the element stiffness matrix elmat. */
  void AssembleElementMatrix(const FiniteElement &el,
				     ElementTransformation &Trans,
				     DenseMatrix &elmat) override 
  {
    int dof = el.GetDof();
    int dim = el.GetDim();
    double w, L, M;
    int width;
     
    MFEM_VERIFY(dim == Trans.GetSpaceDim(), "");
     
    if (dim == 2)
      width = 3;
    if (dim == 3)
      width = 6;
#ifdef MFEM_THREAD_SAFE
    DenseMatrix gshape(dof, dim);
#else
    Emat.SetSize(width);
    dshape.SetSize(dof, dim);
    tBEBmat.SetSize(dim, dim);
    tBmatL.SetSize(dim, width);
    tBmatR.SetSize(dim, width);
    BmatR.SetSize(width,dim);
    EBmat.SetSize(width,dim);
    gtBmat.SetSize(dim,dof*width);
#endif

    elmat.SetSize(dof * dim);
    
    const IntegrationRule *ir = IntRule;
    if (ir == NULL)
      {
	int order = 2 * Trans.OrderGrad(&el); // correct order?
	ir = &IntRules.Get(el.GetGeomType(), order);
      }
     
    elmat = 0.0;

    for (int i = 0; i < ir -> GetNPoints(); i++)
      {
	const IntegrationPoint &ip = ir->IntPoint(i);

	Trans.SetIntPoint(&ip);
	// Each row of the result dshape contains
	// the derivatives of one shape function at the point ip.
	el.CalcPhysDShape(Trans, dshape);
	// Get the transformation Trans for point ip
	// Get the weights associated to point ip
	w = ip.weight * Trans.Weight();
	// Multiply the derivatives by the inverse jacobian matrix
	// to get the derivatives along x, y and z
	// dshape contains these derivatives dshape(dof,dim)
	M = mu->Eval(Trans, ip);
	L = lambda->Eval(Trans, ip);

	// Fill Emat
	Emat = 0.;
	for (int k = 0; k< dim; k++)
	  {
	    // Extra-diagonal terms Emat
	    for (int l = 0; l< dim; l++)
	      Emat(k,l) = L;
	    // Diagonal terms Emat
	    Emat(k,k) = L+2*M;
	  }
	// Diagonal terms
	for (int k = dim; k < width; k++)
	  Emat(k,k)= M;

	// Fill gtBmat
	gtBmat = 0.;
	for (int s = 0; s < dof; s++)
	  {
	    int offset = width*s;
	    // Fill gtBmat that depends on dof 's'
	    // gtBmat(:,offset:offset+width) contains tBmat for dof 's'
	    for (int k = 0; k < dim; k++) 
	      gtBmat(k,offset+k) = dshape(s,k);
	    gtBmat(0,offset+dim) = dshape(s,1); // N_y 
	    gtBmat(1,offset+dim) = dshape(s,0); // N_x 
	    if (dim == 3) {
	      gtBmat(1,offset+4) = dshape(s,2); // N_z
	      gtBmat(2,offset+4) = dshape(s,1); // N_y
	      gtBmat(0,offset+5) = dshape(s,2); // N_z
	      gtBmat(2,offset+5) = dshape(s,0); // N_x
	    }
	  }

	// Perform the main matrix multiplications
	for (int s = 0; s < dof; s++)
	  {
	    tBmatL.CopyCols(gtBmat,width*s,width*s+width-1);
	    for (int t = 0; t < dof; t++)
	      {
		tBmatR.CopyCols(gtBmat,width*t,width*t+width-1);
		BmatR.Transpose(tBmatR);
		Mult(Emat,BmatR,EBmat);
		Mult(tBmatL,EBmat,tBEBmat);
		// Multiply by the weight
		tBEBmat *= w;
		// Add the contribution to the RHS
		for (int k = 0; k < dim; k++)
		  for (int l = 0; l < dim; l++)
		    elmat(dof*k+s, dof*l+t) += tBEBmat(k,l);
	      }
	  }
      }

  }
};

int main(int argc, char *argv[])
{
   PROFILER_ENABLE;
   // 1. Parse command-line options.
   const char *mesh_file = "cube_2mat_per.mesh";
   int order = 1;

   PROFILER_START(0_total);
   PROFILER_START(1_initialize);
   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.AddOption(&tcase, "-t", "--test-case",
                  "number of the case : Exx->11, Eyy->12, Ezz->13 ...\\");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   PROFILER_END(); PROFILER_START(2_read_mesh);
   // 2. Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral or hexahedral elements with the same code.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();

  PROFILER_END(); PROFILER_START(3_refine_mesh);

  for (int i = 0 ; i < 3 ; i++)
    mesh->UniformRefinement();

  PROFILER_END(); PROFILER_START(4_initialize_fem); 
   // 5. Define a finite element space on the mesh. Here we use vector finite
   //    elements, i.e. dim copies of a scalar finite element space. The vector
   //    dimension is specified by the last argument of the FiniteElementSpace
   //    constructor. 
   FiniteElementCollection *fec;
   FiniteElementSpace *fespace;
   fec = new H1_FECollection(order, dim);
   fespace = new FiniteElementSpace(mesh, fec, dim);

   int ndof = fespace->GetTrueVSize() / dim;
   std::cout << "Number of nodes: " << ndof << std::endl;

   //    Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.
   GridFunction x(fespace);
   x = -1.0; 

   // 6. Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the linear elasticity integrator with piece-wise
   //    constants coefficient lambda and mu.
   Vector lambda(mesh->attributes.Max());
   lambda = 100.0;
   if (mesh->attributes.Max() > 1)
     lambda(1) = lambda(0)*2;
   PWConstCoefficient lambda_func(lambda); 
   // Question: pourquoi piece wise constant ?
   // lié aux attributs : constant par zone matériau	
   Vector mu(mesh->attributes.Max());
   mu = 75.0;
   if (mesh->attributes.Max() > 1)
     mu(1) = mu(0)*2;
   PWConstCoefficient mu_func(mu);

   // Determine the list of true (i.e. conforming) essential boundary dofs.
   Array<int> ess_tdof_list, ess_bdr(mesh->bdr_attributes.Max());
   
   ess_bdr = 0;
   cout << "Full periodic" << endl;
   MFEM_VERIFY( (tcase >= 0) && (tcase <= 5) && (dim == 3),"");

   ess_tdof_list.SetSize(dim);
   for (int k = 0; k < dim; k++) 
     {
       int tgdof = 0+k*ndof;
       ess_tdof_list[k] = tgdof;
       x(tgdof) = 0.0;
     }
  

   BilinearForm *a = new BilinearForm(fespace);
   //   a->AddDomainIntegrator(new Elasticity2Integrator(lambda_func,mu_func));
   a->AddDomainIntegrator(new ElasticityIntegrator(lambda_func,mu_func));
  
   // 7. Set up the right-hand side of the FEM linear system.
   LinearForm rhs(fespace);
   rhs.AddDomainIntegrator(new DomainLFGrad2Integrator(lambda_func,mu_func,tcase));
   PROFILER_END(); PROFILER_START(5_solve);

   a->Assemble();
   rhs.Assemble();

   SparseMatrix A;
   Vector B, X;
   a->FormLinearSystem(ess_tdof_list, x, rhs, A, X, B);
   //   GSSmoother M(A);
   CGSolver pcg;
   pcg.SetRelTol(1e-13);
   pcg.SetMaxIter(300);
   pcg.SetPrintLevel(1);
   //   PCG(A, M, B, X, 1, 500, 1e-20, 0.0);
   pcg.Mult(B, X);

   a->RecoverFEMSolution(X, rhs, x);
   
   PROFILER_END(); PROFILER_START(6_postprocess);
   
   // 14. Save the result
   {
     ParaViewDataCollection paraview_dc("lam", mesh);
     paraview_dc.RegisterField("u",&x);
     //     paraview_dc.SetLevelsOfDetail(4);
     paraview_dc.SetCycle(0);
     paraview_dc.SetTime(0.0);
     //     paraview_dc.SetHighOrderOutput(true);
     paraview_dc.Save();
   }
   
   PROFILER_END(); 
   PROFILER_END(); 
   LogProfiler();
   PROFILER_DISABLE;
   return 0;
}

