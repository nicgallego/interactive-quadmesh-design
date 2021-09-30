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
#include <ACG/Utils/ColorCoder.hh>
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
        trimesh.add_property(associatedFace, "Edge already associated with a Face");
        trimesh.add_property(face_color, "halfedge color");
        trimesh.add_property(barycenter, "Barycenter of each Face");
        trimesh.add_property(reference_edge, "Edge associated with face, used for local coord sys");
    }

    ~Crossfield() {
        trimesh_.remove_property(associatedFace);
        trimesh_.remove_property(face_color);
        trimesh_.remove_property(barycenter);
        trimesh_.remove_property(reference_edge);
    }

public:

    void getCrossfield();

    void getEnergy();

private:

    void createCrossfields();

    void getAngleKbetweenTriangles(std::vector<int> &faces);

    void setlocalCoordFrame(std::vector<int> &faces);

    void getBaryCenterAndRefEdge(std::vector<int> &faces);

    void getConstraints(std::vector<int> &asdf);

    OpenMesh::FPropHandleT<TriMesh::Color> face_color;
    OpenMesh::FPropHandleT<std::pair<Point, int>> reference_edge;
    OpenMesh::FPropHandleT<Point> barycenter;
    OpenMesh::EPropHandleT<bool> associatedFace;

    TriMesh &trimesh_;
    std::vector<int> &heInRange_;
    PolyMesh *dualGraph_;
};


#endif //OPENFLIPPER_CROSSFIELD_HH
