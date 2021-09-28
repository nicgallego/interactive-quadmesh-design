//
// Created by wuschelbueb on 15.09.21.
//

#include "Crossfield.hh"

void Crossfield::getCrossfield() {
    createCrossfields();
}

void Crossfield::createCrossfields() {
    setlocalCoordFrame();
}

void Crossfield::setlocalCoordFrame() {
    std::vector<Point> barycenters;
    getBaryCenterAndRefEdge(barycenters);
}

void Crossfield::getBaryCenterAndRefEdge(std::vector<Point> &barycenters) {
    // status needs to be released before using, in case it still has saved some stuff from before
    trimesh_.release_halfedge_status();
    trimesh_.release_face_status();
    // request to change the status
    trimesh_.request_halfedge_status();
    trimesh_.request_face_status();
    std::vector<int> constraints;
    getConstraints(constraints);

    // assign constraint edges to their faces
    for (int i: constraints) {
        Point bCenter;
        int valence = 0;
        OpenMesh::EdgeHandle eh = trimesh_.edge_handle(i);
        // get random halfedge
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(eh, 1);
        if (trimesh_.is_boundary(heh))
            heh = trimesh_.opposite_halfedge_handle(heh);
        OpenMesh::FaceHandle fh = trimesh_.face_handle(heh);
        for (auto fv_it = trimesh_.fv_iter(fh); fv_it.is_valid(); ++fv_it) {
            bCenter += trimesh_.point(*fv_it);
            ++valence;
        }
        barycenters.push_back(bCenter / valence);
        // add barycenter to face
        trimesh_.property(barycenter, fh) = (bCenter / valence);
        // add reference edge to face
        trimesh_.property(reference_edge, fh) = trimesh_.calc_edge_vector(heh);
        // set face and edge as used
        trimesh_.status(heh).set_tagged(true);
        trimesh_.status(fh).set_tagged(true);
    }

    // assign reference edges in range to faces
    for (int i: heInRange_) {
        Point bCenter;
        int valence = 0;
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i);
        if (trimesh_.is_boundary(heh))
            heh = trimesh_.opposite_halfedge_handle(heh);
        OpenMesh::FaceHandle fh = trimesh_.face_handle(heh);
        OpenMesh::HalfedgeHandle nheh = trimesh_.next_halfedge_handle(heh);
        // halfedges or faces which are already tagged, won't be used anymore, no duplicates are allowed
        if ((std::find(heInRange_.begin(), heInRange_.end(), nheh.idx()) != heInRange_.end()) &&
            !trimesh_.status(heh).tagged() && !trimesh_.status(fh).tagged()) {

            for (auto fv_it = trimesh_.fv_iter(fh); fv_it.is_valid(); ++fv_it) {
                bCenter += trimesh_.point(*fv_it);
                ++valence;
            }
            barycenters.push_back(bCenter / valence);
            trimesh_.property(barycenter, fh) = (bCenter / valence);
            trimesh_.property(reference_edge, fh) = trimesh_.calc_edge_vector(heh);
            trimesh_.status(heh).set_tagged(true);
            trimesh_.status(fh).set_tagged(true);
        }
    }
}

// convert array of halfedges to faces and fill constraints with them
void Crossfield::getConstraints(std::vector<int> &constraints) {
    std::vector<int> selectedVertices = MeshSelection::getVertexSelection(&trimesh_);
    std::vector<int> selectedEdges = MeshSelection::getEdgeSelection(&trimesh_);
    std::vector<int> selectedHEdges = MeshSelection::getHalfedgeSelection(&trimesh_);
    std::vector<int> selectedFaces = MeshSelection::getFaceSelection(&trimesh_);

    // if vertex not empty checks neighbouring vertices. if they are also selected, add edge between
    if (!selectedVertices.empty()) {
        for (int i: selectedVertices) {
            OpenMesh::VertexHandle vh = trimesh_.vertex_handle(i);
            for (int j: selectedVertices) {
                OpenMesh::VertexHandle vhj = trimesh_.vertex_handle(j);
                OpenMesh::HalfedgeHandle heh = trimesh_.find_halfedge(vhj, vh);
                // avoids duplicates with std::find
                if (heh.is_valid() &&
                    !(std::find(constraints.begin(), constraints.end(), trimesh_.edge_handle(heh).idx()) !=
                      constraints.end())) {
                    constraints.push_back(trimesh_.edge_handle(heh).idx());
                }
            }
        }
    }

    if (!selectedEdges.empty()) {
        for (int i: selectedEdges) {
            // avoids duplicates with std::find
            if (!(std::find(constraints.begin(), constraints.end(), i) != constraints.end()))
                constraints.push_back(i);
        }
    }

    if (!selectedHEdges.empty()) {
        for (int i: selectedHEdges) {
            OpenMesh::EdgeHandle eh = trimesh_.edge_handle(trimesh_.halfedge_handle(i));
            // avoids duplicates with std::find
            if (!(std::find(constraints.begin(), constraints.end(), eh.idx()) != constraints.end()))
                constraints.push_back(eh.idx());
        }
    }

    if (!selectedFaces.empty()) {
        for (int i: selectedFaces) {
            OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
            for (auto fe_it = trimesh_.fe_iter(fh); fe_it.is_valid(); ++fe_it) {
                // avoids duplicates with std::find
                if (!(std::find(constraints.begin(), constraints.end(), fe_it->idx()) != constraints.end()))
                    constraints.push_back(fe_it->idx());
            }
        }
    }
}



