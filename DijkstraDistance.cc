#include "DijkstraDistance.hh"

std::vector<int> DijkstraDistance::calculateDijkstra(const double refDist) {
    trimesh_.release_vertex_status();
    trimesh_.request_vertex_status();
    std::vector<int> includedVertices;

    initializeDistanceProperty();

    while (true) {
        double totEdgeLen = 0;
        double vertexIndex = getSmallestDistPropVertex(constrainedVertices, refDist);
        // if all vertices are visited the algo stops
        if (vertexIndex == DBL_MAX)
            break;
        OpenMesh::VertexHandle vh = trimesh_.vertex_handle(vertexIndex);
        trimesh_.status(vh).set_tagged(true);
        includedVertices.push_back(vh.idx());
        for (auto voh_it = trimesh_.voh_iter(vh); voh_it.is_valid(); ++voh_it) {
            OpenMesh::VertexHandle vh_neighbour = trimesh_.to_vertex_handle(*voh_it);
            totEdgeLen = trimesh_.property(distance, vh) + trimesh_.calc_edge_length(*voh_it);
            if (!trimesh_.status(vh_neighbour).tagged()
                && totEdgeLen < trimesh_.property(distance, vh_neighbour)) {
                trimesh_.property(distance, vh_neighbour) = totEdgeLen;
            }
        }
    }
    return includedVertices;
}


/*
 * if some face has a vertex that is smaller than the refDist, face gets added to the vertices in Range
 */
void DijkstraDistance::includeBoundaryFaces(std::vector<int> &includedVertices, const double refDist) {
    std::vector<int> tempVertices;
    for (int i: includedVertices) {
        OpenMesh::VertexHandle vh = trimesh_.vertex_handle(i);
        for (auto voh_it = trimesh_.voh_iter(vh); voh_it.is_valid(); ++voh_it) {
            OpenMesh::VertexHandle vh_neighbour = trimesh_.to_vertex_handle(*voh_it);
            if (trimesh_.property(distance, vh_neighbour) != DBL_MAX &&
                trimesh_.property(distance, vh_neighbour) >= refDist &&
                !voh_it->is_boundary()) {
                auto fh = trimesh_.face_handle(*voh_it);
                for (auto fv_it = trimesh_.fv_iter(fh); fv_it.is_valid(); ++fv_it)
                    tempVertices.push_back(fv_it->idx());
            }
        }
    }
    //erase duplicates
    std::set<int> s;
    unsigned size = tempVertices.size();
    for (unsigned i = 0; i < size; ++i) s.insert(tempVertices[i]);
    tempVertices.assign(s.begin(), s.end());
    //add to vertices in range
    for (int i: tempVertices) {
        if (std::find(includedVertices.begin(), includedVertices.end(), i) ==
            includedVertices.end())
            includedVertices.push_back(i);
    }
}

/*
 * gets all halfedges in the range. can be used to create a crossfield
 */
std::vector<int> DijkstraDistance::getHEinRange(const std::vector<int> &includedVertices, const double refDist,
                                                const bool inclBoundaryF) {
    for (int i: includedVertices) {
        OpenMesh::VertexHandle vh = trimesh_.vertex_handle(i);
        for (auto voh_it = trimesh_.voh_iter(vh); voh_it.is_valid(); ++voh_it) {
            OpenMesh::VertexHandle vh_neighbour = trimesh_.to_vertex_handle(*voh_it);
            if (trimesh_.property(distance, vh) < refDist &&
                trimesh_.property(distance, vh_neighbour) < refDist) {
                includedHEdges.push_back(voh_it->idx());
            } else if (trimesh_.property(distance, vh_neighbour) < DBL_MAX && inclBoundaryF)
                includedHEdges.push_back(voh_it->idx());
        }
    }
    return includedHEdges;
}

void DijkstraDistance::colorizeEdges(const std::vector<int> &includedHEdges) {
    // define colors
    TriMesh::Color green = {0, 1, 0, 1};
    TriMesh::Color white = {1, 1, 1, 1};
    // colorize all edges white
    for (OpenMesh::EdgeHandle eh: trimesh_.edges()) {
        trimesh_.set_color(eh, white);
    }

    // colorize edges where vertices have a smaller distance than the refDist blue
    // and edges where the refDist is bigger but some vertices of the face are smaller than refDist green
    for (int i: includedHEdges) {
        OpenMesh::HalfedgeHandle ehh = trimesh_.halfedge_handle(i);
        trimesh_.set_color(trimesh_.edge_handle(ehh), green);
    }
}

// check every vertex and return vertex with the smallest "distance" property which is still unvisited
double DijkstraDistance::getSmallestDistPropVertex(std::vector<int> &allVertices, const double refDist) {
    double minDistance = DBL_MAX;
    double anyVertex = DBL_MAX;
    for (int i: allVertices) {
        OpenMesh::VertexHandle vh = trimesh_.vertex_handle(i);
        if (!trimesh_.status(vh).tagged() &&
            trimesh_.property(distance, vh) < minDistance &&
            trimesh_.property(distance, vh) < refDist) {
            minDistance = trimesh_.property(distance, vh);
            anyVertex = vh.idx();
        }
    }
    return anyVertex;
}


void DijkstraDistance::initializeDistanceProperty() {
    std::vector<int> selectedVertices = MeshSelection::getVertexSelection(&trimesh_);
    std::vector<int> selectedEdges = MeshSelection::getEdgeSelection(&trimesh_);
    std::vector<int> selectedHEdges = MeshSelection::getHalfedgeSelection(&trimesh_);
    std::vector<int> selectedFaces = MeshSelection::getFaceSelection(&trimesh_);
    const double infiniteDistance = DBL_MAX;
    const double zeroDistance = 0;

    // give all vertices the property "distance" with infiniteDistance
    for (OpenMesh::VertexHandle vh: trimesh_.vertices()) {
        trimesh_.property(distance, vh) = infiniteDistance;
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
        constrainedVertices.push_back(vh.idx());
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

void DijkstraDistance::getDualGraphDijkstra(std::vector<int> &includedHEdges) {
    trimesh_.release_face_status();
    trimesh_.request_face_status();
    std::vector<int> constraintHEdges;
    for (int i: constrainedVertices) {
        checkVOH(i, constraintHEdges);
    }
    if (constraintHEdges.empty()) {
        throw std::runtime_error("getDualGraphDijkstra: constraintHEdges can't be empty!\n");
    }
    std::vector<int> faces = createFaceVector(constraintHEdges);
    setDualGraphDijkstra(faces);
}

void DijkstraDistance::checkVOH(const int i, std::vector<int> &constraintHEdges) {
    OpenMesh::VertexHandle vh = trimesh_.vertex_handle(i);
    for (auto voh_it = trimesh_.voh_iter(vh); voh_it.is_valid(); ++voh_it) {
        if (std::find(includedHEdges.begin(), includedHEdges.end(), voh_it->idx()) != includedHEdges.end()
            &&
            std::find(constraintHEdges.begin(), constraintHEdges.end(), voh_it->idx()) == constraintHEdges.end()) {
            constraintHEdges.push_back(voh_it->idx());
        }
    }
}

std::vector<int> DijkstraDistance::createFaceVector(const std::vector<int> constraintHEdges) {
    std::vector<int> faces;
    trimesh_.release_halfedge_status();
    trimesh_.release_face_status();
    // request to change the status
    trimesh_.request_halfedge_status();
    trimesh_.request_face_status();

    // assign constraint edges to their faces
    for (int i: constraintHEdges) {
        setFacesVecWithRefHe(i, faces);
    }

    // assign reference edges in range to faces
    for (int i: includedHEdges) {
        setFacesVec(i, faces);

    }
    return faces;
}

void DijkstraDistance::setFacesVecWithRefHe(const int i, std::vector<int> &faces) {
    OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i);
    if (trimesh_.is_boundary(heh)) {
        heh = trimesh_.opposite_halfedge_handle(heh);
    }
    OpenMesh::FaceHandle fh = trimesh_.face_handle(heh);
    if (!trimesh_.status(heh).tagged() && !trimesh_.status(fh).tagged()) {
        trimesh_.property(distanceBaryCenter, fh) = 0.0;
        trimesh_.property(origin_constraint, fh) = fh.idx();
        // set face and edge as used
        trimesh_.status(heh).set_tagged(true);
        trimesh_.status(fh).set_tagged(true);
        faces.push_back(fh.idx());
    }
}

void DijkstraDistance::setFacesVec(const int i, std::vector<int> &faces) {
    OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i);
    if (trimesh_.is_boundary(heh))
        heh = trimesh_.opposite_halfedge_handle(heh);
    OpenMesh::FaceHandle fh = trimesh_.face_handle(heh);
    // avoids including faces which only have one halfedge in heInRange
    OpenMesh::HalfedgeHandle nheh = trimesh_.next_halfedge_handle(heh);
    // halfedges or faces which are already tagged, won't be used anymore, no duplicates are allowed
    if ((std::find(includedHEdges.begin(), includedHEdges.end(), nheh.idx()) != includedHEdges.end()) &&
        (!trimesh_.status(heh).tagged() && !trimesh_.status(fh).tagged())) {
        // both halfedges need to be tagged, else it is possible opposite faces share an edge
        trimesh_.property(distanceBaryCenter, fh) = DBL_MAX;
        trimesh_.property(origin_constraint, fh) = fh.idx();
        trimesh_.status(heh).set_tagged(true);
        trimesh_.status(trimesh_.opposite_halfedge_handle(heh)).set_tagged(true);
        trimesh_.status(fh).set_tagged(true);
        faces.push_back(fh.idx());
    }
}

void DijkstraDistance::setDualGraphDijkstra(const std::vector<int> &faces) {
    trimesh_.release_face_status();
    trimesh_.request_face_status();
    while (true) {
        double distance = 0;
        int faceIdx = getFaceWithSmallestDist(faces);
        // if all vertices are visited the algo stops
        if (faceIdx == DBL_MAX)
            break;
        OpenMesh::FaceHandle fh = trimesh_.face_handle(faceIdx);
        trimesh_.status(fh).set_tagged(true);
        for (auto ff_it = trimesh_.ff_iter(fh); ff_it.is_valid(); ++ff_it) {
            if (!trimesh_.property(distanceBaryCenter, *ff_it) == 0) {
                Point origin = trimesh_.calc_face_centroid(fh);
                Point neighbour = trimesh_.calc_face_centroid(*ff_it);
                Point temp = neighbour - origin;
                double distance_norm = temp.norm();
                distance = trimesh_.property(distanceBaryCenter, fh) + distance_norm;
                if (!trimesh_.status(*ff_it).tagged()
                    && distance < trimesh_.property(distanceBaryCenter, *ff_it)) {
                    trimesh_.property(distanceBaryCenter, *ff_it) = distance;
                    trimesh_.property(origin_constraint, *ff_it) = fh.idx();
                }
            }
        }
    }
}

int DijkstraDistance::getFaceWithSmallestDist(const std::vector<int> &faces) {
    double minDistance = DBL_MAX;
    double faceIdx = DBL_MAX;
    for (int i: faces) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        if (!trimesh_.status(fh).tagged() &&
            trimesh_.property(distanceBaryCenter, fh) < minDistance) {
            minDistance = trimesh_.property(distanceBaryCenter, fh);
            faceIdx = fh.idx();
        }
    }
    return faceIdx;
}