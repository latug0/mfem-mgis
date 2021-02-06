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
  class DenseMatrix;
  class Mesh;
  class FiniteElementSpace;
  class FiniteElementCollection;
  class FiniteElement;
  class ElementTransformation;
  class IntegrationRule;
  class NonlinearForm;
  class NonlinearFormIntegrator;
#ifdef MFEM_USE_MPI
  class ParMesh;
  class ParFiniteElementSpace;
  class ParNonlinearForm;
#endif

}  // end of namespace mfem


namespace mfem_mgis {
#ifdef MFEM_USE_MPI
  using mMesh = mfem::ParMesh;
  using mFiniteElementSpace = mfem::ParFiniteElementSpace;
  using mNonlinearForm = mfem::ParNonlinearForm;
  using mNonlinearFormIntegrator = mfem::NonlinearFormIntegrator;
#else
  using mMesh = mfem::Mesh;
  using mFiniteElementSpace = mfem::FiniteElementSpace;
  using mNonlinearForm = mfem::NonlinearForm;
  using mNonlinearFormIntegrator = mfem::NonlinearFormIntegrator;
#endif
}
#endif /* LIB_MFEM_MGIS_MFEM_FORWARD_HXX */
