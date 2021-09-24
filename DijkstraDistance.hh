#ifndef OPENFLIPPER_DIJKSTRADISTANCE_HH
#define OPENFLIPPER_DIJKSTRADISTANCE_HH

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
    DijkstraDistance(TriMesh &trimesh) : trimesh_{trimesh} {
        // request color change
        trimesh_.request_edge_colors();
        // add properties to mesh and give them a name
        trimesh_.add_property(edge_color, "edge_color");
        trimesh_.add_property(distance, "distance to origin");
        trimesh_.add_property(visited, "vertex already visited");
    }

    ~DijkstraDistance() {
        trimesh_.remove_property(edge_color);
        trimesh_.remove_property(distance);
        trimesh_.remove_property(visited);
    }

public:


    std::vector<int> calculateDijkstra(const double refDist);

    void includeBoundaryFaces(std::vector<int> &verticesInRange, const double refDist);

    std::vector<int>
    getHEinRange(const std::vector<int> &verticesInRange, const double refDist, const bool inclBoundaryF);

    void colorizeArea(const std::vector<int> &heInRange);

private:
    TriMesh &trimesh_;
    // define properties
    OpenMesh::EPropHandleT<TriMesh::Color> edge_color;
    OpenMesh::VPropHandleT<double> distance;
    OpenMesh::VPropHandleT<bool> visited;

    double getSmallestDistPropVertex(std::vector<int> &allVertices, const double refDist);

    void initializeDistanceProperty(std::vector<int> &allVertices);

    void initializeSelectedVertices(std::vector<int> &selectedVertices, const double zeroDistance);

    void initializeSelectedEdges(std::vector<int> &selectedEdges, const double zeroDistance);

    void initializeSelectedHEdges(std::vector<int> &selectedHEdges, const double zeroDistance);

    void initializeSelectedFaces(std::vector<int> &selectedFaces, const double zeroDistance);
};


#endif //OPENFLIPPER_DIJKSTRADISTANCE_HH
