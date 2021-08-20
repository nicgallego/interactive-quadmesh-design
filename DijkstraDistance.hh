#ifndef OPENFLIPPER_DIJKSTRADISTANCE_HH
#define OPENFLIPPER_DIJKSTRADISTANCE_HH

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <OpenFlipper/BasePlugin/BaseInterface.hh>


class DijkstraDistance {
public:
    DijkstraDistance(TriMesh &trimesh) : trimesh_{trimesh} {
        // request color change
        trimesh_.request_edge_colors();
        // add properties to mesh and give them a name
        trimesh_.add_property(edge_color, "edge_color");
        trimesh_.add_property(edge_length, "edge_length");
        trimesh_.add_property(distance, "distance to origin");
        trimesh_.add_property(pred, "predecessor");
        trimesh_.add_property(visited, "vertex already visited");
    }

    ~DijkstraDistance() {
        trimesh_.remove_property(edge_color);
        trimesh_.remove_property(edge_length);
        trimesh_.remove_property(distance);
        trimesh_.remove_property(pred);
        trimesh_.remove_property(visited);
    }

public:

    void colorizeEdgeSelection();

    void calculateDijkstra(const double refDist);

private:
    TriMesh &trimesh_;
    // define properties
    OpenMesh::VPropHandleT<double> edge_length;
    OpenMesh::EPropHandleT<TriMesh::Color> edge_color;
    OpenMesh::VPropHandleT<double> distance;
    OpenMesh::VPropHandleT<double> pred;
    OpenMesh::VPropHandleT<bool> visited;

    void initialize(std::vector<int> &allVertices);

    double smallestDistVertex(std::vector<int> &allVertices);
};


#endif //OPENFLIPPER_DIJKSTRADISTANCE_HH
