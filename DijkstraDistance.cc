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
    // get selected vertices in vector
    std::vector<int> selectedVertices = MeshSelection::getVertexSelection(&trimesh_);
    for (int i: selectedVertices) {
        // add handle to vertices so you can work with them
        OpenMesh::VertexHandle vh = trimesh_.vertex_handle(i);
        std::cout << "neighbouring vertices of " << vh << ": " << trimesh_.point(vh) << " are:\n";
        // find adjacent edges
        for (TriMesh::VertexEdgeIter ve_it = trimesh_.ve_iter(vh); ve_it.is_valid(); ++ve_it) {
            // write to the property
            trimesh_.property(edge_length, *ve_it) = trimesh_.calc_edge_length(*ve_it);
            trimesh_.property(edge_color, *ve_it);
            trimesh_.set_color(*ve_it, green);
            float LL = trimesh_.property(edge_length, *ve_it);
            std::cout << *ve_it << " length is: " << LL << " ||\t";
        }
        std::cout << "\n=============\n";
    }
}
