//
// Created by wuschelbueb on 15.09.21.
//

#ifndef OPENFLIPPER_CROSSFIELD_HH
#define OPENFLIPPER_CROSSFIELD_HH

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <MeshTools/MeshSelectionT.hh>
#include <ObjectTypes/PolyMesh/PolyMesh.hh>
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
    Crossfield(TriMesh &trimesh)
            : trimesh_{trimesh} {
    }

    ~Crossfield() {}

public:

    void getCrossfield();

    void getEnergy();

private:

    void createCrossfields();

    void setlocalCoordFrame();

    void getCircumCenter();

    std::vector<int> getConstraints();


    OpenMesh::FPropHandleT<int> face_associated_edge;
    OpenMesh::FPropHandleT<ACG::Vec3d> circumcenter;
    OpenMesh::FPropHandleT<ACG::Vec3d> crossfield;

    TriMesh &trimesh_;
    PolyMesh dualGraph_;
};


#endif //OPENFLIPPER_CROSSFIELD_HH
