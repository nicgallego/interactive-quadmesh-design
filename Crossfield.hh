//
// Created by wuschelbueb on 15.09.21.
//

#ifndef OPENFLIPPER_CROSSFIELD_HH
#define OPENFLIPPER_CROSSFIELD_HH

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <OpenMesh/Core/Utils/Property.hh>
#include <ObjectTypes/PolyMesh/PolyMesh.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <MeshTools/MeshSelectionT.hh>
#include <ACG/Geometry/Types/PlaneT.hh>
#include <ACG/Scenegraph/TransformNode.hh>
#include <ACG/Utils/ColorCoder.hh>
#include <CoMISo/Solver/ConstrainedSolver.hh>
#include <gmm/gmm.h>
#include <functional>
#include <iostream>
#include <vector>
#include <float.h>

class Crossfield {
public:
    using Point = ACG::Vec3d;
public:

    /*
     * is the same as:>
     * Crossfield(TriMesh &trimesh, std::vector<int> &heInRange) {
     *      trimesh_ = trimesh;
     *      heInRange_ = heInRange;
     * }
     */
    Crossfield(TriMesh &trimesh, std::vector<int> &heInRange)
            : trimesh_{trimesh}, heInRange_{heInRange} {
        trimesh_.add_property(theta, "angle to rotate");
        trimesh_.add_property(face_color, "Face");
        trimesh_.add_property(heh_color, "Halfedge");
        trimesh_.add_property(x_axis, "x axis of local coord system");
        trimesh_.add_property(y_axis, "y axis of local coord system");
        trimesh_.add_property(x_axis_r, "x axis of local coord system rotated");
        trimesh_.add_property(y_axis_r, "y axis of local coord system rotated");
        trimesh_.add_property(pos_matrixA, "position matrix A");
        trimesh_.add_property(reference_HEdge, "halfedge, x-axis local coord sys");
    }

    ~Crossfield() {
//        trimesh_.remove_property(face_color);
//        trimesh_.remove_property(pos_matrixA);
//        trimesh_.remove_property(barycenter);
//        trimesh_.remove_property(reference_HEdge);
//        trimesh_.remove_property(x_axis);
//        trimesh_.remove_property(y_axis);
//        trimesh_.remove_property(x_axis_r);
//        trimesh_.remove_property(y_axis_r);
//        trimesh_.remove_property(theta);
    }

    void getCrossfield();

    void getEnergy();

private:

    void createCrossfields();

    double
    getEnergy(const std::map<int, double> &edgeKappa, const std::vector<int> &faces, const std::vector<double> &_rhs,
              const std::vector<double> &_x);

    gmm::row_matrix<gmm::wsvector<double>>
    getConstraints(const std::map<int, double> &edgeKappa, const std::vector<int> &constrainedEdges,
                   const std::vector<int> &faces);

    std::vector<int> getIdxToRound(const std::map<int, double> &edgeKappa, const std::vector<int> &faces);

    std::vector<double> getRHS(const std::map<int, double> &edgeKappa, const std::vector<int> &faces);

    gmm::col_matrix<gmm::wsvector<double>>
    getMatrixA(const std::vector<int> &faces, const std::map<int, double> &edgeKappa);

    std::map<int, double> getKappa(const std::vector<int> &faces);

    void setlocalCoordFrame(const std::vector<int> &faces);

    void rotateLocalCoordFrame(const std::vector<int> &faces, const std::vector<double> _x);

    void colorFaces(const std::vector<int> &faces);

    void colorHEdges(const std::vector<int> &constrainedEdges);

    std::vector<int> getBaryCenterAndRefEdge(const std::vector<int> &constrainedEdges);

    std::vector<int> getConstraints();

    double getTotalArea(const std::vector<int> &faces);

    void removeProperties ();


    TriMesh &trimesh_;
    std::vector<int> &heInRange_;

    OpenMesh::FPropHandleT<int> face_color;
    OpenMesh::FPropHandleT<int> pos_matrixA;
    OpenMesh::FPropHandleT<int> reference_HEdge;
    OpenMesh::FPropHandleT<Point> x_axis;
    OpenMesh::FPropHandleT<Point> y_axis;
    OpenMesh::FPropHandleT<Point> x_axis_r;
    OpenMesh::FPropHandleT<Point> y_axis_r;
    OpenMesh::FPropHandleT<double> theta;
    OpenMesh::HPropHandleT<int> heh_color;

};


#endif //OPENFLIPPER_CROSSFIELD_HH
