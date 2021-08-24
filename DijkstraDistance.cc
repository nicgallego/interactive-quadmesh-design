#include "DijkstraDistance.hh"

//bool sortcol(const std::vector<double> &v1,
//             const std::vector<double> &v2) {
//    return v1[1] < v2[1];
//}

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

    initialize(allVertices);

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
    //testing vector to see if everything worked
//    std::vector<std::vector<double>> distVector;
//    for (OpenMesh::VertexHandle vh: trimesh_.vertices()) {
//        double dist = trimesh_.property(distance, vh);
//        distVector.push_back({(double) vh.idx(), dist});
//    };
//    std::vector<std::vector<double>> vertexProperties;
//
//    for (OpenMesh::VertexHandle vh: trimesh_.vertices()) {
//        vertexProperties.push_back({(double) vh.idx(), trimesh_.property(distance, vh)});
//    };
//
//    std::sort(vertexProperties.begin(), vertexProperties.end(), sortcol);
//
//    for (int i = 0; i < (int) vertexProperties.size(); ++i) {
//        for (int j = 0; j < 2; ++j)
//            std::cout << vertexProperties[i][j] << " ";
//        std::cout << std::endl;
//    }
//    std::cout << "==============================\n";

}

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

