#ifndef OPENFLIPPER_DIJKSTRADISTANCE_HH
#define OPENFLIPPER_DIJKSTRADISTANCE_HH

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <OpenFlipper/BasePlugin/BaseInterface.hh>



class DijkstraDistance {
public:
    DijkstraDistance(TriMesh &trimesh): trimesh_{trimesh} {
        // request color change
        trimesh_.request_edge_colors();
        // add properties to mesh and give them a name
        trimesh_.add_property(edge_color, "edge_color");
        trimesh_.add_property(edge_length, "edge_length");
    }

    ~DijkstraDistance() {
        trimesh_.remove_property(edge_color);
        trimesh_.remove_property(edge_length);
    }

public:

    void colorizeEdgeSelection();


private:
    TriMesh &trimesh_;
    // define properties
    OpenMesh::EPropHandleT<float> edge_length;
    OpenMesh::EPropHandleT<TriMesh::Color> edge_color;
};


#endif //OPENFLIPPER_DIJKSTRADISTANCE_HH
