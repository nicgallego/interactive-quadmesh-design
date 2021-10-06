//
// Created by wuschelbueb on 15.09.21.
//

#include "Crossfield.hh"


void Crossfield::getCrossfield() {
    createCrossfields();
}

void Crossfield::createCrossfields() {
    // first column has the constraints, second column contains the rest of the faces in range; do i need that?
    std::vector<std::vector<int>> importantFaces;
    std::vector<int> faces;
    setlocalCoordFrame(faces);
    getKappa(faces);
}

void Crossfield::getKappa(std::vector<int> &faces) {
    for (int i: faces) {
        std::cout << "ith triangle: " << i << std::endl;
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        for (auto ff_it = trimesh_.ff_iter(fh); ff_it.is_valid(); ++ff_it) {
            std::cout << "new neighbour: " << ff_it->idx() << std::endl;
            // neighbour needs to be in faces
            if (trimesh_.status(*ff_it).tagged()) {
                double alpha, beta, kappa = DBL_MAX;
                std::pair<int, int> commonEdge = {INT_MAX, INT_MAX};
                // get index of ref edges
                int refEdgeMain = trimesh_.property(reference_edge, fh).second;
                int refEdgeNeigh = trimesh_.property(reference_edge, *ff_it).second;
                // get the common edge between the triangles
                for (auto fhe_it_main = trimesh_.fh_iter(fh); fhe_it_main.is_valid(); ++fhe_it_main) {
                    for (auto fhe_it_neigh = trimesh_.fh_iter(*ff_it); fhe_it_neigh.is_valid(); ++fhe_it_neigh) {
                        if (trimesh_.opposite_halfedge_handle(*fhe_it_main).idx() == fhe_it_neigh->idx()) {
                            commonEdge.first = fhe_it_main->idx();
                            commonEdge.second = fhe_it_neigh->idx();
                            // as soon as the common edge is found, the loop is obsolete
                            goto END_LOOP;
                        }
                    }
                }
                END_LOOP:
                // there are 9 different scenarios on how reference edges can be placed in two adjacent triangles
                // refEdge shares common edge
                if (refEdgeMain == commonEdge.first) {
                    alpha = 0;
                    if (refEdgeNeigh == commonEdge.second)
                        beta = 0;
                    if (refEdgeNeigh ==
                        trimesh_.prev_halfedge_handle(trimesh_.halfedge_handle(commonEdge.second)).idx())
                        beta = trimesh_.calc_sector_angle(trimesh_.halfedge_handle(refEdgeNeigh));
                    if (refEdgeNeigh ==
                        trimesh_.next_halfedge_handle(trimesh_.halfedge_handle(commonEdge.second)).idx())
                        beta = trimesh_.calc_sector_angle(trimesh_.halfedge_handle(commonEdge.second));

                    // refEdge is the previous halfedge of the common edge
                } else if (refEdgeMain ==
                           trimesh_.prev_halfedge_handle(trimesh_.halfedge_handle(commonEdge.first)).idx()) {
                    alpha = trimesh_.calc_sector_angle(trimesh_.halfedge_handle(refEdgeMain));
                    if (refEdgeNeigh == commonEdge.second)
                        beta = 0;
                    if (refEdgeNeigh ==
                        trimesh_.prev_halfedge_handle(trimesh_.halfedge_handle(commonEdge.second)).idx())
                        beta = M_PI - trimesh_.calc_sector_angle(trimesh_.halfedge_handle(refEdgeNeigh));
                    if (refEdgeNeigh ==
                        trimesh_.next_halfedge_handle(trimesh_.halfedge_handle(commonEdge.second)).idx())
                        beta = trimesh_.calc_sector_angle(trimesh_.halfedge_handle(commonEdge.second));

                    // refEdge is the next halfedge of the common edge
                } else if (refEdgeMain ==
                           trimesh_.next_halfedge_handle(trimesh_.halfedge_handle(commonEdge.first)).idx()) {
                    alpha = trimesh_.calc_sector_angle(trimesh_.halfedge_handle(commonEdge.first));
                    if (refEdgeNeigh == commonEdge.second)
                        beta = 0;
                    if (refEdgeNeigh ==
                        trimesh_.prev_halfedge_handle(trimesh_.halfedge_handle(commonEdge.second)).idx())
                        beta = trimesh_.calc_sector_angle(trimesh_.halfedge_handle(refEdgeNeigh));
                    if (refEdgeNeigh ==
                        trimesh_.next_halfedge_handle(trimesh_.halfedge_handle(commonEdge.second)).idx())
                        beta = M_PI - trimesh_.calc_sector_angle(trimesh_.halfedge_handle(commonEdge.second));
                } else {
                    std::cerr << "Something went wrong, there needs to be a ref edge in the main triangle!\n";
                }
                // do something with kappa
                kappa = alpha + beta;
                std::cout << "Kappa is: " << kappa << " (" << kappa * 180 / M_PI << ") with alpha: " << alpha
                          << " and beta: " << beta << std::endl;
            }
        }
    }
}


void Crossfield::setlocalCoordFrame(std::vector<int> &faces) {
    // 4 works for vlr 55 works for hr file
    int shrinkingFactor = 4;
    getBaryCenterAndRefEdge(faces);

    for (int i: faces) {
        OpenMesh::VertexHandle vhandle[5];
        std::vector<OpenMesh::VertexHandle> face_handles;
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        Point x_scale = trimesh_.property(reference_edge, fh).first;
        Point x_scale_n = x_scale.normalize();
        Point x_scale_from_bary_p = trimesh_.property(barycenter, fh) + x_scale_n / shrinkingFactor;
        Point x_scale_from_bary_m = trimesh_.property(barycenter, fh) - x_scale_n / shrinkingFactor;

        Point face_n = trimesh_.calc_face_normal(fh).normalize();
        Point y_scale = OpenMesh::cross(x_scale_n, face_n);
        Point y_scale_n = y_scale.normalize();
        Point y_scale_from_bary_p = trimesh_.property(barycenter, fh) + y_scale_n / shrinkingFactor;
        Point y_scale_from_bary_n = trimesh_.property(barycenter, fh) - y_scale_n / shrinkingFactor;

        vhandle[0] = trimesh_.add_vertex(trimesh_.property(barycenter, fh));
        vhandle[1] = trimesh_.add_vertex(x_scale_from_bary_p);
        vhandle[2] = trimesh_.add_vertex((x_scale_from_bary_m));
        vhandle[3] = trimesh_.add_vertex(y_scale_from_bary_p);
        vhandle[4] = trimesh_.add_vertex((y_scale_from_bary_n));

        // create x axis
        face_handles.clear();
        face_handles.push_back(vhandle[1]);
        face_handles.push_back(vhandle[0]);
        face_handles.push_back(vhandle[2]);
        trimesh_.add_face(face_handles);
        //create y axis
        face_handles.clear();
        face_handles.push_back(vhandle[3]);
        face_handles.push_back(vhandle[0]);
        face_handles.push_back(vhandle[4]);
        trimesh_.add_face(face_handles);
    }
}

void Crossfield::getBaryCenterAndRefEdge(std::vector<int> &faces) {
    // status needs to be released before using, in case it still has saved some stuff from before
    trimesh_.release_halfedge_status();
    trimesh_.release_face_status();
    // request to change the status
    trimesh_.request_halfedge_status();
    trimesh_.request_face_status();
    std::vector<int> constrainedHEdges;
    getConstraints(constrainedHEdges);

    // assign constraint edges to their faces
    for (int i: constrainedHEdges) {
        Point bCenter = {0, 0, 0};
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
        // add barycenter to face
        trimesh_.property(barycenter, fh) = (bCenter / valence);
        // add reference edge to face
        trimesh_.property(reference_edge, fh) = {trimesh_.calc_edge_vector(heh), heh.idx()};
        // set face and edge as used
        trimesh_.status(heh).set_tagged(true);
        trimesh_.status(trimesh_.opposite_halfedge_handle(heh)).set_tagged(true);
        trimesh_.status(fh).set_tagged(true);
        faces.push_back(fh.idx());
    }

    // assign reference edges in range to faces
    for (int i: heInRange_) {
        Point bCenter = {0, 0, 0};
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
            trimesh_.property(barycenter, fh) = (bCenter / valence);
            trimesh_.property(reference_edge, fh) = {trimesh_.calc_edge_vector(heh), heh.idx()};
            // both halfedges need to be tagged, else it is possible opposite faces share an edge
            trimesh_.status(heh).set_tagged(true);
            trimesh_.status(trimesh_.opposite_halfedge_handle(heh)).set_tagged(true);
            trimesh_.status(fh).set_tagged(true);
            faces.push_back(fh.idx());
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



