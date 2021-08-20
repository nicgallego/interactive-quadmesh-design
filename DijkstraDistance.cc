#include <iostream>
#include <vector>
#include <float.h>
#include "DijkstraDistance.hh"
#include <MeshTools/MeshSelectionT.hh>
#include "OpenFlipper/BasePlugin/PluginFunctions.hh"
#include "OpenMesh/Core/Utils/PropertyManager.hh"

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
    std::vector<int> allVertices;
    initialize(allVertices);

    while (true) {
        double totEdgeLen = 0;
        double vertexIndex = smallestDistVertex(allVertices);
        // if all vertices are visited the algo stops
        if (vertexIndex == DBL_MAX)
            break;
        OpenMesh::VertexHandle vh = trimesh_.vertex_handle(vertexIndex);
        trimesh_.property(visited, vh) = true;

        for (auto voh_it = trimesh_.voh_iter(vh); voh_it.is_valid(); ++voh_it) {
            OpenMesh::VertexHandle vh_neighbour = trimesh_.to_vertex_handle(*voh_it);
            totEdgeLen = trimesh_.property(distance, vh) + trimesh_.calc_edge_length(*voh_it);
            if (!trimesh_.property(visited, vh_neighbour) && totEdgeLen < trimesh_.property(distance, vh_neighbour))
                trimesh_.property(distance, vh_neighbour) = totEdgeLen;
        }
    }
    //testing vector to see if everything worked
    std::vector<std::vector<double>> distVector;
    for (OpenMesh::VertexHandle vh: trimesh_.vertices()) {
        double dist = trimesh_.property(distance, vh);
        distVector.push_back({(double) vh.idx(), dist});
    };
}

double DijkstraDistance::smallestDistVertex(std::vector<int> &allVertices) {
    double minDistance = DBL_MAX;
    double vertex = DBL_MAX;
    for (int i: allVertices) {
        OpenMesh::VertexHandle vh = trimesh_.vertex_handle(i);
        if (!trimesh_.property(visited, vh) && trimesh_.property(distance, vh) < minDistance) {
            minDistance = trimesh_.property(distance, vh);
            vertex = vh.idx();
        }
    }
    return vertex;
}

void DijkstraDistance::initialize(std::vector<int> &allVertices) {
    std::vector<int> selectedVertices = MeshSelection::getVertexSelection(&trimesh_);
    const double infiniteDistance = DBL_MAX;
    const double zeroDistance = 0;

    for (OpenMesh::VertexHandle vh: trimesh_.vertices()) {
        trimesh_.property(distance, vh) = infiniteDistance;
        trimesh_.property(visited, vh) = false;
    };

    // set the property of the distance zero
    for (int i: selectedVertices) {
        OpenMesh::VertexHandle vh = trimesh_.vertex_handle(i);
        trimesh_.property(distance, vh) = zeroDistance;
    }

    // create list with all vertices
    for (auto vh: trimesh_.vertices())
        allVertices.push_back(vh.idx());
}

