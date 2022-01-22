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
        trimesh_.add_property(constraint_angle, "angle constraint");
        trimesh_.add_property(heh_color, "Halfedge");
        trimesh_.add_property(xa, "x axis of local coord system");
        trimesh_.add_property(ya, "y axis of local coord system");
        trimesh_.add_property(x_vec_field, "x vector field");
        trimesh_.add_property(y_vec_field, "y vector field");
        trimesh_.add_property(x_vec_field_r, "x vector field rotated");
        trimesh_.add_property(y_vec_field_r, "y vector field rotated");
        trimesh_.add_property(pos_matrixA, "position matrix A");
        trimesh_.add_property(referenceHeIdx, "ref he, x-axis local coord sys");
    }

    ~Crossfield() {
//        trimesh_.remove_property(face_color);
//        trimesh_.remove_property(pos_matrixA);
//        trimesh_.remove_property(barycenter);
//        trimesh_.remove_property(referenceHeIdx);
//        trimesh_.remove_property(x_vec_field);
//        trimesh_.remove_property(y_vec_field);
//        trimesh_.remove_property(x_vec_field_r);
//        trimesh_.remove_property(y_vec_field_r);
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
    getConstraintMatrix(const std::map<int, double> &edgeKappa, const std::vector<int> &constrainedEdges,
                        const std::vector<int> &faces);

    std::vector<int> getIdxToRound(const std::map<int, double> &edgeKappa, const std::vector<int> &faces);

    std::vector<double> getRHS(const std::map<int, double> &edgeKappa, const std::vector<int> &faces);

    gmm::col_matrix<gmm::wsvector<double>>
    getMatrixA(const std::vector<int> &faces, const std::map<int, double> &edgeKappa);

    std::map<int, double> getMapHeKappa(const std::vector<int> &faces);

    void getStatusNeigh(const OpenMesh::FaceHandle fh, const OpenMesh::FaceHandle fh_neigh, std::map<int, double> &edgeKappa);
    void addKappaHeToMap(const std::pair<int, int> commonEdge, std::map<int, double> &edgeKappa);
    std::pair<int,int> getCommonEdgeBetweenTriangles(const OpenMesh::FaceHandle fh, const OpenMesh::FaceHandle fh_neigh, const int &refEdgeMain, const int &refEdgeNeigh);
    double getKappa(const int refEdgeMain, const int refEdgeNeigh, const std::pair<int,int> commonEdge);


    void setlocalCoordFrame(const std::vector<int> &faces);

    void rotateLocalCoordFrame(const std::vector<int> &faces, const std::vector<double> _x);

    void colorFaces(const std::vector<int> &faces);

    void colorHEdges(const std::vector<int> &constrainedEdges);

    std::vector<int> getReferenceEdge(const std::vector<int> &constrainedHEdges);

    void setFacesVecWithRefHe(const int i, int &temp, std::vector<int> &faces);

    void setFacesVec(const int i, const int temp, std::vector<int> &faces);

    std::vector<int> getConstraints();

    void getSelectedVertices(std::vector<int> &constraints);

    void getSelectedEdges(std::vector<int> &constraints);

    void getSelectedHEdges(std::vector<int> &constraints);

    void getSelectedFaces(std::vector<int> &constraints);

    double getTotalArea(const std::vector<int> &faces);

    std::vector<int> getConstrainedEdges(const std::vector<int> &constrainedHEdges);

    void removeProperties();


    TriMesh &trimesh_;
    std::vector<int> &heInRange_;

    OpenMesh::FPropHandleT<int> face_color;
    OpenMesh::FPropHandleT<int> pos_matrixA;
    OpenMesh::FPropHandleT<int> referenceHeIdx;
    OpenMesh::FPropHandleT<double> constraint_angle;
    OpenMesh::FPropHandleT<Point> xa;
    OpenMesh::FPropHandleT<Point> ya;
    OpenMesh::FPropHandleT<Point> x_vec_field;
    OpenMesh::FPropHandleT<Point> y_vec_field;
    OpenMesh::FPropHandleT<Point> x_vec_field_r;
    OpenMesh::FPropHandleT<Point> y_vec_field_r;
    OpenMesh::FPropHandleT<double> theta;
    OpenMesh::HPropHandleT<int> heh_color;

};


#endif //OPENFLIPPER_CROSSFIELD_HH
