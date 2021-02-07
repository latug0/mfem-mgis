/*!
 * \file   include/MFEMMGIS/MFEMForward.hxx
 * \brief    
 * \author Thomas Helfer
 * \date   11/06/2020
 */

#ifndef LIB_MFEM_MGIS_MFEM_FORWARD_HXX
#define LIB_MFEM_MGIS_MFEM_FORWARD_HXX

namespace mfem {

  class Vector;
  class GridFunction;
  class DenseMatrix;
  class Mesh;
  class FiniteElementSpace;
  class FiniteElementCollection;
  class FiniteElement;
  class ElementTransformation;
  class IntegrationRule;
  class NonlinearForm;
  class NonlinearFormIntegrator;
  class ParGridFunction;
  class ParMesh;
  class ParFiniteElementSpace;
  class ParNonlinearForm;

}  // end of namespace mfem


namespace mfem_mgis {
#ifdef LIB_MPI
  using mMesh = mfem::ParMesh;
  using mFiniteElementSpace = mfem::ParFiniteElementSpace;
  using mNonlinearForm = mfem::ParNonlinearForm;
  using mNonlinearFormIntegrator = mfem::NonlinearFormIntegrator;
  using mGridFunction = mfem::ParGridFunction;
#else
  #error "Should not be in sequential mode"
  using mMesh = mfem::Mesh;
  using mFiniteElementSpace = mfem::FiniteElementSpace;
  using mNonlinearForm = mfem::NonlinearForm;
  using mNonlinearFormIntegrator = mfem::NonlinearFormIntegrator;
  using mGridFunction = mfem::GridFunction;
#endif
}
#endif /* LIB_MFEM_MGIS_MFEM_FORWARD_HXX */
