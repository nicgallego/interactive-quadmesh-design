#include <iostream>
#include <vector>
#include <float.h>
#include "DijkstraDistance.hh"
#include <MeshTools/MeshSelectionT.hh>
#include "OpenFlipper/BasePlugin/PluginFunctions.hh"

void DijkstraDistance::colorizeEdgeSelection() {
    BaseObjectData *object;
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
    }
}

void DijkstraDistance::calculateDijkstra(const double refDist) {
    const double infiniteDistance = DBL_MAX;
    const double nullDistance = 0;
    bool rmPrevElem {false};
    double totEdgeLen;
    // set visited status of all edges false and distance to a big number
    for (OpenMesh::VertexHandle vh: trimesh_.vertices()) {
        trimesh_.property(visited, vh) = false;
        trimesh_.property(distance, vh) = infiniteDistance;
    };
    // get vertex selection
    std::vector<int> selectedVertices = MeshSelection::getVertexSelection(&trimesh_);
    std::vector<int> visVertices = selectedVertices;

    // set the property of the selection to vistied and distance zero
    for (int i: selectedVertices) {
        OpenMesh::VertexHandle vh = trimesh_.vertex_handle(i);
        trimesh_.property(distance, vh) = nullDistance;
        trimesh_.property(visited, vh) = true;
    }
    while (true) {
        for (int i: selectedVertices) {
            OpenMesh::VertexHandle vh = trimesh_.vertex_handle(i);

            //TODO: is auto replaceable with a VertexOHalfedgeIter?
            for (auto voh_it = trimesh_.voh_iter(vh); voh_it.is_valid(); ++voh_it) {
                OpenMesh::VertexHandle vh_next = trimesh_.to_vertex_handle(voh_it);
                totEdgeLen = trimesh_.property(distance, vh) +
                             trimesh_.calc_edge_length(voh_it);
                // if not smaller, skip distance
                if (totEdgeLen < refDist) {
                    // if totEdgeLen is shorter than property and visited = false -> set distance = totEdgeLen and visited = true
                    if (totEdgeLen < trimesh_.property(distance, vh_next) && !trimesh_.property(visited, vh_next)) {
                        trimesh_.property(distance, vh_next) = totEdgeLen;
                        trimesh_.property(visited, vh_next) = true;
                        // not a good idea to make infinite loops!!!!!!!!!!
                        // find workaround
                        visVertices.push_back(vh_next.idx());
                    }
                    // if totEdgeLen is shorter than property and visited = true -> replace totEdgeLen
                    if (totEdgeLen < trimesh_.property(distance, vh_next) && trimesh_.property(visited, vh_next))
                        trimesh_.property(distance, vh_next) = totEdgeLen;
                        visVertices.push_back(vh_next.idx());
                }
            }
        }
        // temporary
        if (totEdgeLen > refDist)
            break;


    }
    // get selection
    // cycle through selection; set property of selected one zero
    // create property bool visited = false
    // create property double distance = infinite
    // set property of all the other vertices infinite and put them in a list
    // get edge length to next vertices; update properties to distance (sum up distance)
    // go to next vertex and repeat for his neighbouring vertices
    // check if distance is shorter than refDist
    // remove the ones with visited status from list
    std::cout << "this is the value: " << refDist << std::endl;

}
