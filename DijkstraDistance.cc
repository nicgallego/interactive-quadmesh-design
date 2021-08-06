#include <iostream>
#include <vector>
#include "DijkstraDistance.hh"
#include <MeshTools/MeshSelectionT.hh>
#include "OpenFlipper/BasePlugin/PluginFunctions.hh"


void DijkstraDistance::colorizeEdgeSelection() {
    BaseObjectData* object;
    // define colors
    TriMesh::Color green = {0, 255, 0, 255};
    TriMesh::Color white = {255, 255, 255, 255};
    // colorize all edges white
    for (OpenMesh::EdgeHandle eh: trimesh_.edges()) {
        // write to the property
        trimesh_.property(edge_color, eh);
        trimesh_.set_color(eh, white);
    }
    // get selected vertices
    std::vector<int> selectedVertices = MeshSelection::getVertexSelection(&trimesh_);
    for (int i: selectedVertices) {
        // add handle to vertices so you can work with them
        OpenMesh::VertexHandle vh = trimesh_.vertex_handle(i);
        // find adjacent edges
        for (TriMesh::VertexEdgeIter ve_it = trimesh_.ve_iter(vh); ve_it.is_valid(); ++ve_it) {
            // write to the property
//            trimesh_.property(edge_length, *ve_it) = trimesh_.calc_edge_length(*ve_it);
            trimesh_.property(edge_color, *ve_it);
            trimesh_.set_color(*ve_it, green);
        }
        std::cout << "\n=============\n";
    }
}

void DijkstraDistance::calculateDijkstra(const double ref_distance) {
    double infinite_distance = 1000;
    double start_distance = 0;
    for (OpenMesh::VertexHandle vh: trimesh_.vertices()) {
        trimesh_.property(visited, vh) = false;
        trimesh_.property(distance, vh) = infinite_distance;
    };
    std::vector<int> selectedVertices = MeshSelection::getVertexSelection(&trimesh_);
    for (int i: selectedVertices) {
        OpenMesh::VertexHandle vh = trimesh_.vertex_handle(i);
        trimesh_.property(distance, vh) = start_distance;
//        trimesh_.property(visited, vh) = true;
    }

    for(int i: selectedVertices) {
        OpenMesh::VertexHandle vh = trimesh_.vertex_handle(i);
        std::cout << "Our initial vertex is: " << vh << " and has the following neighbouring vertices:" << std::endl;
        for (TriMesh::VertexEdgeIter ve_it = trimesh_.ve_iter(vh); ve_it.is_valid(); ++ve_it) {
            // write to the property
            OpenMesh::HalfedgeHandle voh = trimesh_.halfedge_handle(ve_it,1);
            OpenMesh::VertexHandle vh_next = trimesh_.to_vertex_handle(voh);
            trimesh_.property(edge_length, vh_next) = trimesh_.calc_edge_length(*ve_it);

            std::cout << vh_next << " with the length " << trimesh_.property(edge_length, vh_next) << "\t";

        }
        if(trimesh_.property(distance, vh) != 0) {};
    }
    // get selection
    // cycle through selection; set property of selected one zero
    // create property bool visited = false
    // create property double distance = infinite
    // set property of all the other vertices infinite and put them in a list
    // get edge length to next vertices; update properties to distance (sum up distance)
    // go to next vertex and repeat for his neighbouring vertices
    // check if distance is shorter than ref_distance
    // remove the ones with visited status from list
    std::cout << "this is the value: " << ref_distance << std::endl;

}
