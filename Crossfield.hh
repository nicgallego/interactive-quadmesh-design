//
// Created by wuschelbueb on 15.09.21.
//

#ifndef OPENFLIPPER_CROSSFIELD_HH
#define OPENFLIPPER_CROSSFIELD_HH

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>
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
    }

    ~Crossfield() {
    }

    void getCrossfield();

    void getEnergy();

private:

    void createCrossfields();

    double
    getEnergy(const std::map<int, double> &edgeKappa, const std::vector<int> &faces, const std::vector<double> &_rhs,
              const std::vector<double> &_x);

    gmm::row_matrix<gmm::wsvector<double>>
    getConstraintMatrix(const std::map<int, double> &edgeKappa, const std::vector<int> &constraintHalfEdges,
                        const std::vector<int> &faces);

    std::vector<int> getIdxToRound(const std::map<int, double> &edgeKappa, const std::vector<int> &faces);

    std::vector<double> getRHS(const std::map<int, double> &edgeKappa, const std::vector<int> &faces);

    void getRhsFirstHalf(double const totalArea, const std::vector<int> &faces, std::vector<double> &_rhs,
                         const std::map<int, double> &edgeKappa);

    void getSum(double const totalArea, const int i, std::vector<double> &_rhs, const std::map<int, double> &edgeKappa);

    double setSum(OpenMesh::FaceHandle fh, OpenMesh::HalfedgeHandle fh_it, const std::map<int, double> &edgeKappa,
                  double const totalArea);

    void getRhsSecondHalf(double const totalArea, std::vector<double> &_rhs,
                          const std::map<int, double> &edgeKappa, const int facesPlusOne);

    gmm::col_matrix<gmm::wsvector<double>>
    getMatrixA(const std::vector<int> &faces, const std::map<int, double> &edgeKappa);

    void renumberFace(OpenMesh::FaceHandle fh, int &counter);

    std::map<int, double> getMapHeKappa(const std::vector<int> &faces);

    void getStatusNeigh(const OpenMesh::FaceHandle fh, const OpenMesh::FaceHandle fh_neigh,
                        std::map<int, double> &edgeKappa);

    void addKappaHeToMap(const std::pair<int, int> commonEdge, const double kappa, std::map<int, double> &edgeKappa);

    std::pair<int, int>
    getCommonEdgeBetweenTriangles(const OpenMesh::FaceHandle fh, const OpenMesh::FaceHandle fh_neigh,
                                  const int &refEdgeMain, const int &refEdgeNeigh);

    double getKappa(const int refEdgeMain, const int refEdgeNeigh, const std::pair<int, int> commonEdge);

    void setlocalCoordFrame(const std::vector<int> &faces);

    void rotateLocalCoordFrame(const std::vector<int> &faces, const std::vector<double> _x);

    void colorFaces(const std::vector<int> &faces);

    void colorHEdges(const std::vector<int> &constrainedEdges);

    std::vector<int> getReferenceEdge(const std::vector<int> &constrainedHEdges);

    void setFacesVecWithRefHe(const int i, std::vector<int> &faces);

    void setFacesVec(const int i, std::vector<int> &faces);

    std::vector<int> getConstraints();

    void getSelectedVertices(std::vector<int> &constraints);

    void getSelectedEdges(std::vector<int> &constraints);

    void getSelectedHEdges(std::vector<int> &constraints);

    void getSelectedFaces(std::vector<int> &constraints);

    double getTotalArea(const std::vector<int> &faces);

    std::vector<int> getConstraintEdges(const std::vector<int> &constrainedHEdges);


    TriMesh &trimesh_;
    std::vector<int> &heInRange_;
};


#endif //OPENFLIPPER_CROSSFIELD_HH
