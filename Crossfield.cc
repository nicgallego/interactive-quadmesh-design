//
// Created by wuschelbueb on 15.09.21.
//

#include "Crossfield.hh"

void Crossfield::getCrossfield() {
    createCrossfields();
}

void Crossfield::createCrossfields() {
    std::vector<int> constrainedEdges = getConstraints();
    std::vector<int> faces = getBaryCenterAndRefEdge(constrainedEdges);
    setlocalCoordFrame(faces);
    std::map<int, double> edgeKappa = getKappa(faces);
    gmm::col_matrix<gmm::wsvector<double>> _A = getMatrixA(faces, edgeKappa);
    gmm::row_matrix<gmm::wsvector<double>> _constraints = getConstraints(edgeKappa, constrainedEdges, faces);
    std::vector<double> _x(edgeKappa.size() + faces.size(), 0.0);
    std::vector<double> _rhs = getRHS(edgeKappa, faces);
    std::vector<int> _idx_to_round = getIdxToRound(edgeKappa, faces);

    std::vector<double> _rhsOld = _rhs;
    std::vector<int> _idx_to_roundOld = _idx_to_round;

    COMISO::ConstrainedSolver csolver;
    csolver.misolver().set_iter_full(false);
    csolver.misolver().set_local_iters(50000);
    csolver.misolver().set_cg_iters(20);
    csolver.misolver().set_local_error(1e-3);
    csolver.misolver().set_cg_error(1e-3);
    csolver.misolver().set_multiple_rounding();
    std::cout << "Dimensions of parameters:\n" << "Constraints Rows:\t" << gmm::mat_nrows(_constraints)
              << " and columns: " << gmm::mat_ncols(_constraints) << std::endl
              << "A Rows:\t\t\t\t" << gmm::mat_nrows(_A) << " and columns: " << gmm::mat_ncols(_A) << std::endl
              << "Size of _x:\t\t\t" << _x.size() << std::endl << "Size of _rhs:\t\t" << _rhs.size() << std::endl
              << "Size of edgeKappa:\t" << edgeKappa.size() << std::endl << "Size of idx:\t\t" << _idx_to_round.size()
              << std::endl;
    csolver.solve(_constraints, _A, _x, _rhs, _idx_to_round);
//    double energy_before = getEnergy(edgeKappa, faces, _rhsOld, _x);
//    double energy_after = getEnergy(edgeKappa, faces, _rhs, _x);
//    std::cout << "The energy before the smoothing is:\t" << energy_before << std::endl
//              << "The energy after the smoothing is:\t" << energy_after << std::endl;
    std::cout << "hello world\n";

}

double Crossfield::getEnergy(const std::map<int, double> &edgeKappa, const std::vector<int> &faces,
                             const std::vector<double> &_rhs, const std::vector<double> &_x) {
    int pj_start = faces.size();
    double sum = 0.0;
    for (auto i: edgeKappa) {

        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i.first);
        OpenMesh::HalfedgeHandle oheh = trimesh_.opposite_halfedge_handle(heh);
        OpenMesh::FaceHandle fh = trimesh_.face_handle(heh);
        OpenMesh::FaceHandle ofh = trimesh_.face_handle(oheh);
        int position = trimesh_.property(pos_matrixA, fh);
        int opposite_position = trimesh_.property(pos_matrixA, ofh);
//        std::cout << "theta_i: " << _x[position] << std::endl << "kappa: " << _rhs[pj_start] << std::endl << "periodjump: "
//                  << _x[pj_start] << std::endl << "theta_j: " << _x[opposite_position] << std::endl;
        sum += pow((_x[position] + _rhs[pj_start] + (M_PI / 2 * _x[pj_start]) - _x[opposite_position]), 2);
        pj_start++;
    }
    return sum;
}

gmm::row_matrix<gmm::wsvector<double>>
Crossfield::getConstraints(const std::map<int, double> &edgeKappa, const std::vector<int> &constrainedEdges,
                           const std::vector<int> &faces) {
    //cNplusOne = angle between local coordinate x-axis and constraint
    int cNplusOne = 1, counter = 0;
    int n_row = constrainedEdges.size(), n_col = edgeKappa.size() + faces.size();
    gmm::row_matrix<gmm::wsvector<double>> _constraints(n_row, n_col + cNplusOne);
    for (int i: constrainedEdges) {
        OpenMesh::EdgeHandle ehCons = trimesh_.edge_handle(i);
        OpenMesh::HalfedgeHandle hehCons = trimesh_.halfedge_handle(ehCons, 1);
        if (trimesh_.is_boundary(hehCons))
            hehCons = trimesh_.opposite_halfedge_handle(hehCons);
        OpenMesh::FaceHandle fh = trimesh_.face_handle(hehCons);
        int position = trimesh_.property(pos_matrixA, fh);
        int idx = trimesh_.property(reference_HEdge, fh).second;
        OpenMesh::HalfedgeHandle hehRef = trimesh_.halfedge_handle(idx);
        _constraints(counter, position) = 1.0;
        if (hehCons.idx() == hehRef.idx())
            _constraints(counter, n_col) = -1 * 0;
        if (hehCons.idx() == trimesh_.prev_halfedge_handle(hehRef).idx())
            _constraints(counter, n_col) = -1 * trimesh_.calc_sector_angle(hehCons);
        if (hehCons.idx() == trimesh_.next_halfedge_handle(hehRef).idx())
            _constraints(counter, n_col) = -1 * trimesh_.calc_sector_angle(hehRef);
        counter++;
    }
    gmm::clean(_constraints, 1E-10);
    std::cout << "We have a matrix with " << gmm::mat_nrows(_constraints) <<
              " rows and " << gmm::mat_ncols(_constraints) << " columns:\n" << _constraints << std::endl;
    return _constraints;
}

std::vector<int> Crossfield::getIdxToRound(const std::map<int, double> &edgeKappa, const std::vector<int> &faces) {
    std::vector<int> _idx_to_round;
    int pj_start = faces.size();
    for (const auto &i: edgeKappa)
        _idx_to_round.push_back(pj_start++);
    return _idx_to_round;
}

std::vector<double> Crossfield::getRHS(const std::map<int, double> &edgeKappa, const std::vector<int> &faces) {
    std::vector<double> _rhs(faces.size() + edgeKappa.size());
    double totalArea = getTotalArea(faces), edge_weight = 0.0;
    int facesPlusOne = faces.size();
    for (int i: faces) {
        double sum = 0;
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        for (auto fh_it = trimesh_.fh_iter(fh); fh_it.is_valid(); ++fh_it) {
            OpenMesh::HalfedgeHandle opposite_heh = trimesh_.opposite_halfedge_handle(*fh_it);
            OpenMesh::SmartFaceHandle opposite_fh;
            // this condition is necessary because an opposite face doesn't exist if the opposite heh is boundary
            if (!trimesh_.is_boundary(opposite_heh)) {
                opposite_fh = trimesh_.face_handle(trimesh_.opposite_halfedge_handle(*fh_it));
                edge_weight =
                        ((trimesh_.calc_face_area(fh) / 3) + (trimesh_.calc_face_area(opposite_fh) / 3)) / totalArea;
                auto it = edgeKappa.find(fh_it->idx());
                auto it2 = edgeKappa.find(opposite_heh.idx());
                // it doesn't matter that edge_weight doesn't have a value, since this condition is never true, if edge_weight is valueless
                if (it != edgeKappa.end())
                    sum += (it->second * edge_weight);
                else if (it2 != edgeKappa.end())
                    sum -= (it2->second * edge_weight);
                else if (it != edgeKappa.end() && it2 != edgeKappa.end()) {
                    throw std::runtime_error(
                            "Opposite Halfedges can't both be in edgeKappa!\n");
                }
            }
        }
        int position = trimesh_.property(pos_matrixA, fh);
        _rhs[position] = sum;
    }
    for (const auto &i: edgeKappa) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i.first);
        OpenMesh::FaceHandle fh = trimesh_.face_handle(heh);
        OpenMesh::FaceHandle opposite_fh = trimesh_.face_handle(trimesh_.opposite_halfedge_handle(heh));
        edge_weight = ((trimesh_.calc_face_area(fh) / 3) + (trimesh_.calc_face_area(opposite_fh) / 3)) / totalArea;
        _rhs[facesPlusOne++] = (i.second * edge_weight);
    }
    // scalar multiplication with vector, b*-1 = -b
    gmm::scale(_rhs, -1.0);
    return _rhs;
}

gmm::col_matrix<gmm::wsvector<double>>
Crossfield::getMatrixA(const std::vector<int> &faces, const std::map<int, double> &edgeKappa) {
    trimesh_.release_halfedge_status();
    trimesh_.release_face_status();
    // request to change the status
    trimesh_.request_halfedge_status();
    trimesh_.request_face_status();
    int counter = 0, iteration = 0, pj_start = faces.size(), n = edgeKappa.size() + faces.size();
    double edge_weight = 0, totalArea = getTotalArea(faces);
    gmm::col_matrix<gmm::wsvector<double>> _A(n, n);
    // indexes the faces from 0 to n
    // this is needed in case a face index is higher than the matrix size
    for (const auto &i: edgeKappa) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i.first);
        OpenMesh::FaceHandle fh1 = trimesh_.face_handle(heh);
        OpenMesh::FaceHandle fh2 = trimesh_.face_handle(trimesh_.opposite_halfedge_handle(heh));
        if (!trimesh_.status(fh1).tagged()) {
            trimesh_.property(pos_matrixA, fh1) = counter;
            trimesh_.status(fh1).set_tagged(true);
            counter++;
        }
        if (!trimesh_.status(fh2).tagged()) {
            trimesh_.property(pos_matrixA, fh2) = counter;
            trimesh_.status(fh2).set_tagged(true);
            counter++;
        }
    }
    // fills up sparse column matrix _A
    for (const auto &i: edgeKappa) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i.first);
        OpenMesh::FaceHandle fh1 = trimesh_.face_handle(heh);
        OpenMesh::FaceHandle fh2 = trimesh_.face_handle(trimesh_.opposite_halfedge_handle(heh));
        edge_weight = ((trimesh_.calc_face_area(fh1) / 3) + (trimesh_.calc_face_area(fh2) / 3)) / totalArea;
        int pos_i = trimesh_.property(pos_matrixA, fh1);
        int pos_j = trimesh_.property(pos_matrixA, fh2);
        // |  2*w_ij   -2*w_ij     pi*w_ij |
        // | -2*w_ij    2*w_ij    -pi*w_ij |
        // | pi*w_ij  -pi*w_ij pi^2/2*w_ij |
        _A(pos_i, pos_i) = 2.0 * edge_weight;
        _A(pos_j, pos_j) = 2.0 * edge_weight;
        _A(pos_i, pos_j) = -2.0 * edge_weight;
        _A(pos_j, pos_i) = -2.0 * edge_weight;
        _A(pos_i, pj_start + iteration) = M_PI * edge_weight;
        _A(pj_start + iteration, pos_i) = M_PI * edge_weight;
        _A(pos_j, pj_start + iteration) = -M_PI * edge_weight;
        _A(pj_start + iteration, pos_j) = -M_PI * edge_weight;
        _A(pj_start + iteration, pj_start + iteration) = (pow(M_PI, 2) / 2) * edge_weight;
        iteration++;
    }
    gmm::clean(_A, 1E-10);
    return _A;
}

std::map<int, double> Crossfield::getKappa(const std::vector<int> &faces) {
    std::map<int, double> edgeKappa;
    for (int i: faces) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        for (auto ff_it = trimesh_.ff_iter(fh); ff_it.is_valid(); ++ff_it) {
            // neighbour needs to be in faces
            if (trimesh_.status(*ff_it).tagged()) {
                double alpha = DBL_MAX, beta = alpha, kappa = alpha;
                std::pair<int, int> commonEdge = {INT_MAX, INT_MAX};
                // get index of ref edges
                int refEdgeMain = trimesh_.property(reference_HEdge, fh).second;
                int refEdgeNeigh = trimesh_.property(reference_HEdge, *ff_it).second;
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
                    throw std::runtime_error(
                            "Something went wrong, there needs to be a ref edge in the main triangle!\n");
                }
                // if opposite halfedge isn't in edgeKappa add it else do nothing
                if (edgeKappa.find(
                        trimesh_.opposite_halfedge_handle(trimesh_.halfedge_handle(commonEdge.first)).idx()) ==
                    edgeKappa.end()) {
                    kappa = (alpha + beta);
                    int edgeIndex = commonEdge.first;
                    edgeKappa[edgeIndex] = kappa;
                }
            }
        }
    }
    return edgeKappa;
}


void Crossfield::setlocalCoordFrame(const std::vector<int> &faces) {
    // 4 works for vlr 55 works for hr file
    int shrinkingFactor = 4;

    for (int i: faces) {
        OpenMesh::VertexHandle vhandle[5];
        std::vector<OpenMesh::VertexHandle> face_handles;
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        Point x_scale = trimesh_.property(reference_HEdge, fh).first;
        Point x_scale_n = x_scale.normalize();
        Point x_scale_from_bary_p = trimesh_.property(barycenter, fh) + x_scale_n / shrinkingFactor;
        Point x_scale_from_bary_m = trimesh_.property(barycenter, fh) - x_scale_n / shrinkingFactor;

        Point face_n = trimesh_.calc_face_normal(fh).normalize();
        Point y_scale = OpenMesh::cross(x_scale_n, face_n);
        Point y_scale_n = y_scale.normalize();
        Point y_scale_from_bary_p = trimesh_.property(barycenter, fh) + y_scale_n / shrinkingFactor;
        Point y_scale_from_bary_n = trimesh_.property(barycenter, fh) - y_scale_n / shrinkingFactor;


//        vhandle[0] = trimesh_.add_vertex(trimesh_.property(barycenter, fh));
//        vhandle[1] = trimesh_.add_vertex(x_scale_from_bary_p);
//        vhandle[2] = trimesh_.add_vertex((x_scale_from_bary_m));
//        vhandle[3] = trimesh_.add_vertex(y_scale_from_bary_p);
//        vhandle[4] = trimesh_.add_vertex((y_scale_from_bary_n));
//
//        // create x axis
//        face_handles.clear();
//        face_handles.push_back(vhandle[1]);
//        face_handles.push_back(vhandle[0]);
//        face_handles.push_back(vhandle[2]);
//        trimesh_.add_face(face_handles);
//        //create y axis
//        face_handles.clear();
//        face_handles.push_back(vhandle[3]);
//        face_handles.push_back(vhandle[0]);
//        face_handles.push_back(vhandle[4]);
//        trimesh_.add_face(face_handles);
    }
}

std::vector<int> Crossfield::getBaryCenterAndRefEdge(const std::vector<int> &constrainedEdges) {
    // status needs to be released before using, in case it still has saved some stuff from before
    trimesh_.release_halfedge_status();
    trimesh_.release_face_status();
    // request to change the status
    trimesh_.request_halfedge_status();
    trimesh_.request_face_status();
    std::vector<int> faces;

    // assign constraint edges to their faces
    for (int i: constrainedEdges) {
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
        if (!trimesh_.status(heh).tagged() && !trimesh_.status(fh).tagged()) {
            // add barycenter to face
            trimesh_.property(barycenter, fh) = (bCenter / valence);
            // add reference edge to face
            trimesh_.property(reference_HEdge, fh) = {trimesh_.calc_edge_vector(heh), heh.idx()};
            // set face and edge as used
            trimesh_.status(heh).set_tagged(true);
            trimesh_.status(trimesh_.opposite_halfedge_handle(heh)).set_tagged(true);
            trimesh_.status(fh).set_tagged(true);
            faces.push_back(fh.idx());
        }
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
            trimesh_.property(reference_HEdge, fh) = {trimesh_.calc_edge_vector(heh), heh.idx()};
            // both halfedges need to be tagged, else it is possible opposite faces share an edge
            trimesh_.status(heh).set_tagged(true);
            trimesh_.status(trimesh_.opposite_halfedge_handle(heh)).set_tagged(true);
            trimesh_.status(fh).set_tagged(true);
            faces.push_back(fh.idx());
        }
    }
    return faces;
}

// convert array of halfedges to faces and fill constraints with them
std::vector<int> Crossfield::getConstraints() {
    std::vector<int> constraints;
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
    return constraints;
}

double Crossfield::getTotalArea(const std::vector<int> &faces) {
    double totalarea = 0.0;
    for (int i: faces)
        totalarea += trimesh_.calc_face_area(trimesh_.face_handle(i));
    return totalarea;
}

