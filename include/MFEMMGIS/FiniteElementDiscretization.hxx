/*!
 * \file   include/MFEMMGIS/FiniteElementDiscretization.hxx
 * \brief
 * \author Thomas Helfer
 * \date 16/12/2020
 */

#ifndef LIB_MFEM_MGIS_FINITEELEMENTDISCRETIZATION_HXX
#define LIB_MFEM_MGIS_FINITEELEMENTDISCRETIZATION_HXX

#include <memory>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MFEMForward.hxx"

namespace mfem_mgis {

  /*!
   * \brief a simple class used to:
   * - handle the life time of the mesh and the finite element collection.
   * - create and handle a finite element space
   */
  struct MFEM_MGIS_EXPORT FiniteElementDiscretization {
    /*!
     * \brief constructor
     * \param[in] m: mesh
     * \param[in] c: collection
     * \param[in] d: size of the unknowns
     *
     * \note this methods creates the finite element space.
     */
    FiniteElementDiscretization(
        std::shared_ptr<mMesh>,
        std::shared_ptr<const mfem::FiniteElementCollection>,
        const size_type);
    /*!
     * \brief constructor
     * \param[in] m: mesh
     * \param[in] c: collection
     * \param[in] s: finite element space
     */
    FiniteElementDiscretization(
        std::shared_ptr<mMesh>,
        std::shared_ptr<const mfem::FiniteElementCollection>,
        std::unique_ptr<mFiniteElementSpace>);
    //! \return the mesh
    mMesh& getMesh();
    //! \return the mesh
    const mMesh& getMesh() const;
    //! \return the finite element space
    mFiniteElementSpace& getFiniteElementSpace();
    //! \return the finite element space
    const mFiniteElementSpace& getFiniteElementSpace() const;
    //! \return the finite element collection
    const mfem::FiniteElementCollection& getFiniteElementCollection() const;
    //! \brief destructor
    ~FiniteElementDiscretization();

   private:
    //! \brief mesh
    std::shared_ptr<mMesh> mesh;
    //! \brief finite element collection
    std::shared_ptr<const mfem::FiniteElementCollection> fec;
    //! \brief finite element space
    std::unique_ptr<mFiniteElementSpace> fe_space;
  };  // end of FiniteElementDiscretization

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_FINITEELEMENTDISCRETIZATION_HXX */
