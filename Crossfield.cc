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
    trimesh_.release_halfedge_status();
    trimesh_.release_face_status();
    trimesh_.request_halfedge_status();
    trimesh_.request_face_status();
    std::vector<int> constraints;
    std::vector<int> faces;
    TriMesh::Color white = {1, 1, 1, 1};
    TriMesh::Color green = {0, 1, 0, 1};

    int counter = 0;
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
        trimesh_.property(barycenter, fh) = (bCenter / valence);
        trimesh_.property(reference_edge, fh) = trimesh_.calc_edge_vector(heh);
        trimesh_.status(heh).set_tagged(true);
        trimesh_.status(fh).set_tagged(true);
        faces.push_back(fh.idx());
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
        if ((std::find(heInRange_.begin(), heInRange_.end(), nheh.idx()) != heInRange_.end()) &&
            !trimesh_.status(heh).tagged() && !trimesh_.status(fh).tagged()) {
            ++counter;
            for (auto fv_it = trimesh_.fv_iter(fh); fv_it.is_valid(); ++fv_it) {
                bCenter += trimesh_.point(*fv_it);
                ++valence;
            }
            barycenters.push_back(bCenter / valence);
            trimesh_.property(barycenter, fh) = (bCenter / valence);
            trimesh_.property(reference_edge, fh) = trimesh_.calc_edge_vector(heh);
            trimesh_.status(heh).set_tagged(true);
            trimesh_.status(fh).set_tagged(true);
            faces.push_back(fh.idx());
        }
    }

    for (OpenMesh::FaceHandle fh: trimesh_.faces()) {
        trimesh_.property(face_color, fh);
        trimesh_.set_color(fh, white);
    }

    for (int i: faces) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        trimesh_.property(face_color, fh);
        trimesh_.set_color(fh, green);
    }

    std::cout << "he in range: " << heInRange_.size()
              << "\nconstraints: " << constraints.size()
              << "\nand we counted: " << counter
              << "\nwhich gives us the sum: " << counter + constraints.size()
              << "\nwith faces: " << faces.size() << std::endl;
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
                if (heh.is_valid() &&
                    !(std::find(constraints.begin(), constraints.end(), trimesh_.edge_handle(heh).idx()) !=
                      constraints.end())) {
                    constraints.push_back(trimesh_.edge_handle(heh).idx());
                }
            }

        }
    }

    if (!selectedEdges.empty())
        for (int i: selectedEdges)
            constraints.push_back(i);

    if (!selectedHEdges.empty()) {
        for (int i: selectedHEdges) {
            OpenMesh::EdgeHandle eh = trimesh_.edge_handle(trimesh_.halfedge_handle(i));
            constraints.push_back(eh.idx());
        }
    }

    if (!selectedFaces.empty()) {
        for (int i: selectedFaces) {
            OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
            for (auto fe_it = trimesh_.fe_iter(fh); fe_it.is_valid(); ++fe_it) {
                if (!(std::find(constraints.begin(), constraints.end(), fe_it->idx()) != constraints.end()))
                    constraints.push_back(fe_it->idx());
            }
        }
    }
}



