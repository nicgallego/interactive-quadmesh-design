#include "DijkstraDistance.hh"

std::vector<int> DijkstraDistance::calculateDijkstra(const double refDist) {
    trimesh_.release_vertex_status();
    trimesh_.request_vertex_status();
    auto distance = OpenMesh::VProp<double>(trimesh_, "distance");
    std::vector<int> includedVertices;
    initializeDistanceProperty();

    while (true) {
        double totEdgeLen = 0;
        int vertexIndex = getSmallestDistPropVertex(refDist);
        // if all vertices are visited the algo stops
        if (vertexIndex == INT_MAX)
            break;
        OpenMesh::VertexHandle vh = trimesh_.vertex_handle(vertexIndex);
        trimesh_.status(vh).set_tagged(true);
        includedVertices.push_back(vh.idx());
        for (TriMesh::VertexOHalfedgeIter voh_it = trimesh_.voh_iter(vh); voh_it.is_valid(); ++voh_it) {
            OpenMesh::VertexHandle vh_neighbour = trimesh_.to_vertex_handle(*voh_it);
            totEdgeLen = distance[vh] + trimesh_.calc_edge_length(*voh_it);
//            totEdgeLen = trimesh_.property(distance, vh) + trimesh_.calc_edge_length(*voh_it);
            if (!trimesh_.status(vh_neighbour).tagged() && totEdgeLen < distance[vh_neighbour]) {
                distance[vh_neighbour] = totEdgeLen;
            }
        }
    }
    return includedVertices;
}


/*
 * if some face has a vertex that is smaller than the refDist, face gets added to the vertices in Range
 */
void DijkstraDistance::includeBoundaryFaces(std::vector<int> &includedVertices, const double refDist) {
    auto distance = OpenMesh::VProp<double>(trimesh_, "distance");
    std::vector<int> tempVertices;
    for (int i: includedVertices) {
        OpenMesh::VertexHandle vh = trimesh_.vertex_handle(i);
        for (auto voh_it = trimesh_.voh_iter(vh); voh_it.is_valid(); ++voh_it) {
            OpenMesh::VertexHandle vh_neighbour = trimesh_.to_vertex_handle(*voh_it);
            if (distance[vh_neighbour] != DBL_MAX && distance[vh_neighbour] >= refDist && !voh_it->is_boundary()) {
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
    auto distance = OpenMesh::VProp<double>(trimesh_, "distance");
    std::vector<int> includedHEdges;
    for (int i: includedVertices) {
        OpenMesh::VertexHandle vh = trimesh_.vertex_handle(i);
        for (auto voh_it = trimesh_.voh_iter(vh); voh_it.is_valid(); ++voh_it) {
            OpenMesh::VertexHandle vh_neighbour = trimesh_.to_vertex_handle(*voh_it);
            if (distance[vh] < refDist && distance[vh_neighbour] < refDist) {
                includedHEdges.push_back(voh_it->idx());
            } else if (distance[vh_neighbour] < DBL_MAX && inclBoundaryF)
                includedHEdges.push_back(voh_it->idx());
        }
    }
    return includedHEdges;
}

void DijkstraDistance::colorizeEdges(const std::vector<int> &includedHEdges) {
    // request color change
    trimesh_.request_edge_colors();
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
int DijkstraDistance::getSmallestDistPropVertex(const double refDist) {
    auto distance = OpenMesh::VProp<double>(trimesh_, "distance");
    double minDistance = DBL_MAX;
    int vertexIdx = INT_MAX;
    for (auto vh: trimesh_.vertices()) {
        if (!trimesh_.status(vh).tagged() && distance[vh] < minDistance && distance[vh] < refDist) {
            minDistance = distance[vh];
            vertexIdx = vh.idx();
        }
    }
    return vertexIdx;
}


void DijkstraDistance::initializeDistanceProperty() {
    std::vector<int> selectedVertices = MeshSelection::getVertexSelection(&trimesh_);
    std::vector<int> selectedEdges = MeshSelection::getEdgeSelection(&trimesh_);
    std::vector<int> selectedHEdges = MeshSelection::getHalfedgeSelection(&trimesh_);
    std::vector<int> selectedFaces = MeshSelection::getFaceSelection(&trimesh_);
    auto distance = OpenMesh::VProp<double>(trimesh_, "distance");
    const double infiniteDistance = DBL_MAX;
    const double zeroDistance = 0.0;

    // give all vertices the property "distance" with infiniteDistance
    for (auto vh: trimesh_.vertices()) {
        distance[vh] = infiniteDistance;
    };
    if (!selectedVertices.empty())
        initializeSelectedVertices(selectedVertices, zeroDistance);

    if (!selectedEdges.empty())
        initializeSelectedEdges(selectedEdges, zeroDistance);

    if (!selectedHEdges.empty())
        initializeSelectedHEdges(selectedHEdges, zeroDistance);

    if (!selectedFaces.empty())
        initializeSelectedFaces(selectedFaces, zeroDistance);
}

// sets vertex property distance to zero
void DijkstraDistance::initializeSelectedVertices(std::vector<int> &selectedVertices, const double zeroDistance) {
    auto distance = OpenMesh::VProp<double>(trimesh_, "distance");
    // set the property of the distance zero
    for (int i: selectedVertices) {
        OpenMesh::VertexHandle vh = trimesh_.vertex_handle(i);
        distance[vh] = zeroDistance;
        constraintVertices.push_back(vh.idx());
    }
}

// sets the property "distance" on both adjacent vertices of the selected edges to zero
void DijkstraDistance::initializeSelectedEdges(std::vector<int> &selectedEdges, const double zeroDistance) {
    auto distance = OpenMesh::VProp<double>(trimesh_, "distance");
    for (int i: selectedEdges) {
        OpenMesh::EdgeHandle eh = trimesh_.edge_handle(i);
        OpenMesh::VertexHandle vh1 = trimesh_.to_vertex_handle(trimesh_.halfedge_handle(eh, 0));
        OpenMesh::VertexHandle vh2 = trimesh_.from_vertex_handle(trimesh_.halfedge_handle(eh, 0));
        distance[vh1] = zeroDistance;
        distance[vh2] = zeroDistance;
        constraintVertices.push_back(vh1.idx());
        constraintVertices.push_back(vh2.idx());
    }
}

// sets the property "distance" on both adjacent vertices of the selected halfedges to zero
void DijkstraDistance::initializeSelectedHEdges(std::vector<int> &selectedHEdges, const double zeroDistance) {
    auto distance = OpenMesh::VProp<double>(trimesh_, "distance");
    for (int i: selectedHEdges) {
        OpenMesh::VertexHandle vh1 = trimesh_.to_vertex_handle(trimesh_.halfedge_handle(i));
        OpenMesh::VertexHandle vh2 = trimesh_.from_vertex_handle(trimesh_.halfedge_handle(i));
        distance[vh1] = zeroDistance;
        distance[vh2] = zeroDistance;
        constraintVertices.push_back(vh1.idx());
        constraintVertices.push_back(vh2.idx());
    }
}

// sets the property "distance" on all adjacent vertices of the selected faces to zero
void DijkstraDistance::initializeSelectedFaces(std::vector<int> &selectedFaces, const double zeroDistance) {
    auto distance = OpenMesh::VProp<double>(trimesh_, "distance");
    for (int i: selectedFaces) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        for (auto fv_it = trimesh_.fv_iter(fh); fv_it.is_valid(); fv_it++) {
            distance[*fv_it] = zeroDistance;
            constraintVertices.push_back(fv_it->idx());
        }
    }
}

void DijkstraDistance::getDualGraphDijkstra(std::vector<int> &includedHEdges) {
    trimesh_.release_face_status();
    trimesh_.request_face_status();
    std::vector<int> constraintHEdges;
    for (int i: constraintVertices) {
        auto vh = trimesh_.vertex_handle(i);
        checkVOH(vh, constraintHEdges, includedHEdges);
    }
    if (constraintHEdges.empty()) {
        throw std::runtime_error("getDualGraphDijkstra: constraintHEdges can't be empty!\n");
    }
    std::vector<int> faces = createFaceVector(constraintHEdges, includedHEdges);
    setDualGraphDijkstra(faces);
}

void DijkstraDistance::checkVOH(const OpenMesh::VertexHandle vh, std::vector<int> &constraintHEdges,
                                const std::vector<int> &includedHEdges) {
    for (auto voh_it = trimesh_.voh_iter(vh); voh_it.is_valid(); ++voh_it) {
        checkOccurrenceInVectors(*voh_it, constraintHEdges, includedHEdges);
    }
}

void DijkstraDistance::checkOccurrenceInVectors(const OpenMesh::SmartHalfedgeHandle voh_it,
                                                std::vector<int> &constraintHEdges,
                                                const std::vector<int> &includedHEdges) {
    // can't be in constraintHEgdes and has to be in includedHEdges and to_vertex_handle has to be in constraintVertices
    OpenMesh::VertexHandle vh1 = trimesh_.to_vertex_handle(voh_it);
    if (std::find(includedHEdges.begin(), includedHEdges.end(), voh_it.idx()) != includedHEdges.end()
        && std::find(constraintVertices.begin(), constraintVertices.end(), vh1.idx()) != constraintVertices.end() &&
        std::find(constraintHEdges.begin(), constraintHEdges.end(), voh_it.idx()) == constraintHEdges.end()) {
        constraintHEdges.push_back(voh_it.idx());
    }
}


std::vector<int>
DijkstraDistance::createFaceVector(const std::vector<int> constraintHEdges, const std::vector<int> &includedHEdges) {
    auto origin_constraint = OpenMesh::FProp<int>(trimesh_, "origin_constraint");
    std::vector<int> faces;
    trimesh_.release_halfedge_status();
    trimesh_.release_face_status();
    // request to change the status
    trimesh_.request_halfedge_status();
    trimesh_.request_face_status();
    for (TriMesh::FaceIter f_it = trimesh_.faces_begin(); f_it != trimesh_.faces_end(); ++f_it) {
        origin_constraint[*f_it] = 0;
    }

    // assign constraint edges to their faces
    for (int i: constraintHEdges) {
        setFacesVecWithRefHe(i, faces, constraintHEdges);
    }

    // assign reference edges in range to faces
    for (int i: includedHEdges) {
        setFacesVec(i, faces, includedHEdges);

    }
    return faces;
}

void
DijkstraDistance::setFacesVecWithRefHe(const int i, std::vector<int> &faces, const std::vector<int> constraintHEdges) {
    auto distanceBaryCenter = OpenMesh::FProp<double>(trimesh_, "distanceBaryCenter");
    auto origin_constraint = OpenMesh::FProp<int>(trimesh_, "origin_constraint");
    auto heh = trimesh_.halfedge_handle(i);
    auto oheh = trimesh_.opposite_halfedge_handle(heh);
    auto nheh = trimesh_.next_halfedge_handle(heh);
    auto pheh = trimesh_.prev_halfedge_handle(heh);
    auto onheh = trimesh_.next_halfedge_handle(oheh);
    auto opheh = trimesh_.prev_halfedge_handle(oheh);
    auto fh = trimesh_.face_handle(heh);
    auto ofh = trimesh_.face_handle(oheh);
    if (!trimesh_.is_boundary(heh) && !trimesh_.is_boundary(oheh) && !trimesh_.status(fh).tagged()) {
        if ((std::find(constraintHEdges.begin(), constraintHEdges.end(), nheh.idx()) != constraintHEdges.end()) ||
            (std::find(constraintHEdges.begin(), constraintHEdges.end(), pheh.idx()) != constraintHEdges.end())) {
            addFaceToVector(fh, faces);
        } else if (
                (std::find(constraintHEdges.begin(), constraintHEdges.end(), onheh.idx()) != constraintHEdges.end()) ||
                (std::find(constraintHEdges.begin(), constraintHEdges.end(), opheh.idx()) !=
                 constraintHEdges.end())) {
            addFaceToVector(ofh, faces);
        } else {
            addFaceToVector(fh, faces);
            distanceBaryCenter[ofh] = DBL_MAX;
            origin_constraint[ofh] = fh.idx();
            trimesh_.status(ofh).set_tagged(true);
        }
    } else if (trimesh_.is_boundary(heh) && !trimesh_.status(ofh).tagged()) {
        addFaceToVector(ofh, faces);
    }
}

void DijkstraDistance::addFaceToVector(const OpenMesh::FaceHandle fh, std::vector<int> &faces) {
    auto distanceBaryCenter = OpenMesh::FProp<double>(trimesh_, "distanceBaryCenter");
    auto origin_constraint = OpenMesh::FProp<int>(trimesh_, "origin_constraint");
    distanceBaryCenter[fh] = 0.0;
    origin_constraint[fh] = fh.idx();
    std::cout << "Face " << origin_constraint[fh] << " with Ref HEdge\n";
    trimesh_.status(fh).set_tagged(true);
    faces.push_back(fh.idx());
}

void DijkstraDistance::setFacesVec(const int i, std::vector<int> &faces, const std::vector<int> &includedHEdges) {
    auto origin_constraint = OpenMesh::FProp<int>(trimesh_, "origin_constraint");
    auto distanceBaryCenter = OpenMesh::FProp<double>(trimesh_, "distanceBaryCenter");
    auto heh_color = OpenMesh::HProp<int>(trimesh_, "heh_color");
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
        distanceBaryCenter[fh] = DBL_MAX;
        origin_constraint[fh] = fh.idx();
        std::cout << "Face " << origin_constraint[fh] << " without Ref\n";
        heh_color[heh] = 2;
        trimesh_.status(heh).set_tagged(true);
        trimesh_.status(trimesh_.opposite_halfedge_handle(heh)).set_tagged(true);
        trimesh_.status(fh).set_tagged(true);
        faces.push_back(fh.idx());
    }
}

void DijkstraDistance::setDualGraphDijkstra(const std::vector<int> &faces) {
    trimesh_.release_face_status();
    trimesh_.request_face_status();
    auto distanceBaryCenter = OpenMesh::FProp<double>(trimesh_, "distanceBaryCenter");
    auto origin_constraint = OpenMesh::FProp<int>(trimesh_, "origin_constraint");
    std::map<int, int> dualSpanningTree;
    while (true) {
        double distance = 0.0;
        int faceIdx = getFaceWithSmallestDist(faces);
        // if all vertices are visited the algo stops
        if (faceIdx == INT_MAX)
            break;
        OpenMesh::FaceHandle fh = trimesh_.face_handle(faceIdx);
        trimesh_.status(fh).set_tagged(true);
        for (TriMesh::FaceFaceIter ff_it = trimesh_.ff_iter(fh); ff_it.is_valid(); ++ff_it) {
            if (!distanceBaryCenter[*ff_it] == 0) {
                Point origin = trimesh_.calc_face_centroid(fh);
                Point neighbour = trimesh_.calc_face_centroid(*ff_it);
                Point temp = neighbour - origin;
                double distance_norm = temp.norm();
                distance = distanceBaryCenter[fh] + distance_norm;
                if (!trimesh_.status(*ff_it).tagged()
                    && distance < distanceBaryCenter[*ff_it]) {
                    distanceBaryCenter[*ff_it] = distance;
                    int prevFaceIdx = origin_constraint[fh];
                    origin_constraint[*ff_it] = prevFaceIdx;
                    std::cout << "Face " << *ff_it << " with property " << origin_constraint[*ff_it]
                              << " and distance " << distanceBaryCenter[*ff_it] << std::endl;
                }
            }
        }
    }
}

int DijkstraDistance::getFaceWithSmallestDist(const std::vector<int> &faces) {
    auto distanceBaryCenter = OpenMesh::FProp<double>(trimesh_, "distanceBaryCenter");
    double minDistance = DBL_MAX;
    int faceIdx = INT_MAX;
    for (int i: faces) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        if (!trimesh_.status(fh).tagged() &&
            distanceBaryCenter[fh] < minDistance) {
            minDistance = distanceBaryCenter[fh];
            faceIdx = fh.idx();
        }
    }
    return faceIdx;
}