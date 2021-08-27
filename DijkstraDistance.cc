#include "DijkstraDistance.hh"

void DijkstraDistance::colorizeArea(const double refDist) {
    BaseObjectData *object;
    // define colors
    TriMesh::Color babyblue = {0, 123, 123, 255};
    TriMesh::Color white = {255, 255, 255, 255};
    // colorize all edges white
    for (OpenMesh::EdgeHandle eh: trimesh_.edges()) {
        // write to the property
        trimesh_.property(edge_color, eh);
        trimesh_.set_color(eh, white);
    }
    // colorize edges where vertices have a smaller distance than the refDist
    for (OpenMesh::VertexHandle vh: trimesh_.vertices()) {
        if (trimesh_.property(distance, vh) < refDist) {
            for (auto voh_it = trimesh_.voh_iter(vh); voh_it.is_valid(); ++voh_it) {
                OpenMesh::VertexHandle vh_neighbour = trimesh_.to_vertex_handle(*voh_it);
                if (trimesh_.property(distance, vh_neighbour) < refDist) {
                    trimesh_.set_color(trimesh_.edge_handle(*voh_it), babyblue);
                }
            }
        }
    }
}

std::vector<int> DijkstraDistance::calculateDijkstra(const double refDist) {
    std::vector<int> allVertices;
    std::vector<int> verticesInRange;

    initializeDistanceProperty(allVertices);

    while (true) {
        double totEdgeLen = 0;
        double vertexIndex = getSmallestDistPropVertex(allVertices, refDist);
        // if all vertices are visited the algo stops
        if (vertexIndex == DBL_MAX)
            break;
        OpenMesh::VertexHandle vh = trimesh_.vertex_handle(vertexIndex);
        trimesh_.property(visited, vh) = true;
        verticesInRange.push_back(vh.idx());
        for (auto voh_it = trimesh_.voh_iter(vh); voh_it.is_valid(); ++voh_it) {
            OpenMesh::VertexHandle vh_neighbour = trimesh_.to_vertex_handle(*voh_it);
            totEdgeLen = trimesh_.property(distance, vh) + trimesh_.calc_edge_length(*voh_it);
            if (!trimesh_.property(visited, vh_neighbour) && totEdgeLen < trimesh_.property(distance, vh_neighbour))
                trimesh_.property(distance, vh_neighbour) = totEdgeLen;
        }
    }
    return verticesInRange;
}

// check every vertex and return vertex with the smallest "distance" property which is still unvisited
double DijkstraDistance::getSmallestDistPropVertex(std::vector<int> &allVertices, const double refDist) {
    double minDistance = DBL_MAX;
    double anyVertex = DBL_MAX;
    for (int i: allVertices) {
        OpenMesh::VertexHandle vh = trimesh_.vertex_handle(i);
        if (!trimesh_.property(visited, vh) &&
            trimesh_.property(distance, vh) < minDistance &&
            trimesh_.property(distance, vh) < refDist) {
            minDistance = trimesh_.property(distance, vh);
            anyVertex = vh.idx();
        }
    }
    return anyVertex;
}

void DijkstraDistance::initializeDistanceProperty(std::vector<int> &allVertices) {
    std::vector<int> selectedVertices = MeshSelection::getVertexSelection(&trimesh_);
    std::vector<int> selectedEdges = MeshSelection::getEdgeSelection(&trimesh_);
    std::vector<int> selectedHEdges = MeshSelection::getHalfedgeSelection(&trimesh_);
    std::vector<int> selectedFaces = MeshSelection::getFaceSelection(&trimesh_);
    const double infiniteDistance = DBL_MAX;
    const double zeroDistance = 0;

    // give all vertices the property "distance" with infiniteDistance
    for (OpenMesh::VertexHandle vh: trimesh_.vertices()) {
        trimesh_.property(distance, vh) = infiniteDistance;
        trimesh_.property(visited, vh) = false;
    };
    if (!selectedVertices.empty())
        initializeSelectedVertices(selectedVertices, zeroDistance);

    if (!selectedEdges.empty())
        initializeSelectedEdges(selectedEdges, zeroDistance);

    if (!selectedHEdges.empty())
        initializeSelectedHEdges(selectedHEdges, zeroDistance);

    if (!selectedFaces.empty())
        initializeSelectedFaces(selectedFaces, zeroDistance);


    // create list with all vertices
    for (auto vh: trimesh_.vertices())
        allVertices.push_back(vh.idx());
}

// sets vertex property distance to zero
void DijkstraDistance::initializeSelectedVertices(std::vector<int> &selectedVertices, const double zeroDistance) {
    // set the property of the distance zero
    for (int i: selectedVertices) {
        OpenMesh::VertexHandle vh = trimesh_.vertex_handle(i);
        trimesh_.property(distance, vh) = zeroDistance;
    }
}

// sets the property "distance" on both adjacent vertices of the selected edges to zero
void DijkstraDistance::initializeSelectedEdges(std::vector<int> &selectedEdges, const double zeroDistance) {
    for (int i: selectedEdges) {
        OpenMesh::EdgeHandle eh = trimesh_.edge_handle(i);
        OpenMesh::VertexHandle vh1 = trimesh_.to_vertex_handle(trimesh_.halfedge_handle(eh, 0));
        OpenMesh::VertexHandle vh2 = trimesh_.from_vertex_handle(trimesh_.halfedge_handle(eh, 0));
        trimesh_.property(distance, vh1) = zeroDistance;
        trimesh_.property(distance, vh2) = zeroDistance;
    }
}

// sets the property "distance" on both adjacent vertices of the selected halfedges to zero
void DijkstraDistance::initializeSelectedHEdges(std::vector<int> &selectedHEdges, const double zeroDistance) {
    for (int i: selectedHEdges) {
        OpenMesh::VertexHandle vh1 = trimesh_.to_vertex_handle(trimesh_.halfedge_handle(i));
        OpenMesh::VertexHandle vh2 = trimesh_.from_vertex_handle(trimesh_.halfedge_handle(i));
        trimesh_.property(distance, vh1) = zeroDistance;
        trimesh_.property(distance, vh2) = zeroDistance;
    }
}

// sets the property "distance" on all adjacent vertices of the selected faces to zero
void DijkstraDistance::initializeSelectedFaces(std::vector<int> &selectedFaces, const double zeroDistance) {
    for (int i: selectedFaces) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        for (auto fv_it = trimesh_.fv_iter(fh); fv_it.is_valid(); fv_it++) {
            trimesh_.property(distance, *fv_it) = zeroDistance;
        }
    }
}