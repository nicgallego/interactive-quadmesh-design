//
// Created by wuschelbueb on 15.09.21.
//

#ifndef OPENFLIPPER_CROSSFIELD_HH
#define OPENFLIPPER_CROSSFIELD_HH

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <ObjectTypes/PolyMesh/PolyMesh.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <MeshTools/MeshSelectionT.hh>
#include <ACG/Geometry/Types/PlaneT.hh>
#include <iostream>
#include <vector>
#include <float.h>

class Crossfield {
public:
    using Point = ACG::Vec3d;

    /*
     * is the same as:
     * Crossfield(TriMesh &trimesh) {
     *      trimesh_ = trimesh;
     * }
     */
    Crossfield(TriMesh &trimesh, std::vector<int> &heInRange)
            : trimesh_{trimesh}, heInRange_{heInRange} {
        trimesh.add_property(inUse, "Edge already associated with a Face");
        trimesh.add_property(barycenter, "Barycenter of each Face");
        trimesh.add_property(associated_edge_to_face, "Edge associated with face, used for local coord sys");
    }

    ~Crossfield() {
        trimesh_.remove_property(inUse);
        trimesh_.remove_property(barycenter);
        trimesh_.remove_property(associated_edge_to_face);
    }

public:

    void getCrossfield();

    void getEnergy();

private:

    void createCrossfields();

    void setlocalCoordFrame();

    void getBaryCenter(std::vector<Point> &barycenters);

    void getConstraints(std::vector<int> &constraints);


    OpenMesh::FPropHandleT<OpenMesh::EdgeHandle> associated_edge_to_face;
    OpenMesh::FPropHandleT<Point> barycenter;
    OpenMesh::EPropHandleT<bool> inUse;

    TriMesh &trimesh_;
    std::vector<int> &heInRange_;
    PolyMesh dualGraph_;
};


#endif //OPENFLIPPER_CROSSFIELD_HH
