#ifndef OPENFLIPPER_DIJKSTRADISTANCE_HH
#define OPENFLIPPER_DIJKSTRADISTANCE_HH

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <OpenFlipper/BasePlugin/BaseInterface.hh>
#include <OpenFlipper/BasePlugin/PluginFunctions.hh>
#include <MeshTools/MeshSelectionT.hh>
#include <iostream>
#include <vector>
#include <float.h>
#include <algorithm>


class DijkstraDistance {
public:
    using Point = ACG::Vec3d;
public:
    DijkstraDistance(TriMesh &trimesh) : trimesh_{trimesh} {
        // add properties to mesh and give them a name
        trimesh_.add_property(edge_color, "edge_color");
        trimesh_.add_property(distance, "dijkstra distance to starting point");
        trimesh_.add_property(distanceBaryCenter, "used to calc dijkstra of dual graph");
        trimesh_.add_property(origin_constraint, "origin Point of dual spanning tree");
    }

    ~DijkstraDistance() {
    }

public:
    std::map<int, int> getDualGraphDijkstra(std::vector<int> &includedHEdges);

    std::vector<int> calculateDijkstra(const double refDist);

    void includeBoundaryFaces(std::vector<int> &includedVertices, const double refDist);

    std::vector<int>
    getHEinRange(const std::vector<int> &verticesInRange, const double refDist, const bool inclBoundaryF);

    void colorizeEdges(const std::vector<int> &includedHEdges);

private:

    int getSmallestDistPropVertex(const double refDist);

    void initializeDistanceProperty();

    void initializeSelectedVertices(std::vector<int> &selectedVertices, const double zeroDistance);

    void initializeSelectedEdges(std::vector<int> &selectedEdges, const double zeroDistance);

    void initializeSelectedHEdges(std::vector<int> &selectedHEdges, const double zeroDistance);

    void initializeSelectedFaces(std::vector<int> &selectedFaces, const double zeroDistance);

    void checkVOH(const OpenMesh::VertexHandle vh, std::vector<int> &constraintHEdges,
                  const std::vector<int> &includedHEdges);

    void checkOccurrenceInVectors(const OpenMesh::SmartHalfedgeHandle voh_it, std::vector<int> &constraintHEdges,
                                  const std::vector<int> &includedHEdges);

    std::vector<int> createFaceVector(const std::vector<int> constraintHEdges, const std::vector<int> &includedHEdges);

    void setFacesVecWithRefHe(const int i, std::vector<int> &faces);

    void setFacesVec(const int i, std::vector<int> &faces, const std::vector<int> &includedHEdges);

    std::map<int, int> setDualGraphDijkstra(const std::vector<int> &faces);

    int getFaceWithSmallestDist(const std::vector<int> &faces);

    TriMesh &trimesh_;
    // define properties
    OpenMesh::EPropHandleT<TriMesh::Color> edge_color;
    OpenMesh::VPropHandleT<double> distance;
    OpenMesh::FPropHandleT<double> distanceBaryCenter;
    OpenMesh::FPropHandleT<int> origin_constraint;
    std::vector<int> constraintVertices;
};


#endif //OPENFLIPPER_DIJKSTRADISTANCE_HH
