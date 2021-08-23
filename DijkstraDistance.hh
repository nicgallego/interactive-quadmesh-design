#ifndef OPENFLIPPER_DIJKSTRADISTANCE_HH
#define OPENFLIPPER_DIJKSTRADISTANCE_HH

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <OpenFlipper/BasePlugin/BaseInterface.hh>
#include <OpenFlipper/BasePlugin/PluginFunctions.hh>
#include <MeshTools/MeshSelectionT.hh>
#include <iostream>
#include <vector>
#include <float.h>


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

    void colorizeArea(const double refDist);

    void calculateDijkstra(const double refDist);

private:
    TriMesh &trimesh_;
    // define properties
    OpenMesh::EPropHandleT<TriMesh::Color> edge_color;
    OpenMesh::VPropHandleT<double> distance;
    OpenMesh::VPropHandleT<bool> visited;

    void initialize(std::vector<int> &allVertices);

    double getSmallestDistPropVertex(std::vector<int> &allVertices, const double refDist);
};


#endif //OPENFLIPPER_DIJKSTRADISTANCE_HH
