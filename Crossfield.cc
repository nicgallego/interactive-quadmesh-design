//
// Created by wuschelbueb on 15.09.21.
//

#include "Crossfield.hh"

void Crossfield::getCrossfield() {
    createCrossfields();
}

void Crossfield::createCrossfields() {
//    removeProperties();
    std::vector<int> constrainedHalfEdges = getConstraints();
    std::vector<int> faces = getReferenceEdge(constrainedHalfEdges);
    std::map<int, double> edgeKappa = getMapHeKappa(faces);
    gmm::col_matrix<gmm::wsvector<double>> _A = getMatrixA(faces, edgeKappa);
    gmm::row_matrix<gmm::wsvector<double>> _constraints = getConstraintMatrix(edgeKappa, constrainedHalfEdges, faces);
    std::vector<double> _x(edgeKappa.size() + faces.size(), 0.0);
    std::vector<double> _rhs = getRHS(edgeKappa, faces);
    std::vector<int> _idx_to_round = getIdxToRound(edgeKappa, faces);
    setlocalCoordFrame(faces);
    std::vector<double> _rhsOld = _rhs;

    std::cout << "A (before): " << _A << std::endl;
    for (std::size_t i = 0, max = _rhs.size(); i != max; ++i) {
        std::cout << "rhs position (" << i << ") with value (before calc): " << _rhs[i] << std::endl;
    }
    std::cout << "constraints: " << _constraints << std::endl;

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

    rotateLocalCoordFrame(faces, _x);

    for (std::size_t i = 0, max = _x.size(); i != max; ++i) {
        std::cout << "x position " << i << " with value: " << _x[i] << std::endl;
    }

    std::cout << "A (after): " << _A << std::endl;
    for (std::size_t i = 0, max = _rhs.size(); i != max; ++i) {
        std::cout << "rhs position (" << i << ") with value (after calc): " << _rhs[i] << std::endl;
    }

    std::cout << "constraints after: " << _constraints << std::endl;
    double energy_before = getEnergy(edgeKappa, faces, _rhsOld, _x);
    double energy_after = getEnergy(edgeKappa, faces, _rhs, _x);
    std::cout << "The energy before the smoothing is:\t" << energy_before << std::endl
              << "The energy after the smoothing is:\t" << energy_after << std::endl;
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
//        std::cout << "theta_i: " << _x[position] << std::endl << "kappa: " << _rhs[pj_start] << std::endl
//                  << "periodjump: "
//                  << _x[pj_start] << std::endl << "theta_j: " << _x[opposite_position] << std::endl;
        sum += pow((_x[position] + _rhs[pj_start] + (M_PI / 2 * _x[pj_start]) - _x[opposite_position]), 2);
        pj_start++;
    }
    return sum;
}

gmm::row_matrix<gmm::wsvector<double>>
Crossfield::getConstraintMatrix(const std::map<int, double> &edgeKappa, const std::vector<int> &constrainedHalfEdges,
                                const std::vector<int> &faces) {
    //cNplusOne = angle between local coordinate x-axis and constraint
    int cNplusOne = 1, counter = 0;
    // this results in have the rows of the constrainedhalfedges
    std::vector<int> constrainedEdges = getConstrainedEdges(constrainedHalfEdges);
    int n_row = constrainedEdges.size(), n_col = edgeKappa.size() + faces.size();
    gmm::row_matrix<gmm::wsvector<double>> _constraints(n_row, n_col + cNplusOne);
    for (int i: constrainedEdges) {
        OpenMesh::HalfedgeHandle hehCons = trimesh_.halfedge_handle(i);
        if (trimesh_.is_boundary(hehCons)) {

        } else {

            OpenMesh::FaceHandle fh = trimesh_.face_handle(hehCons);
            double constraint = trimesh_.property(constraint_angle, fh);
            int position = trimesh_.property(pos_matrixA, fh);
            int idx = trimesh_.property(referenceHeIdx, fh);
            OpenMesh::HalfedgeHandle hehRef = trimesh_.halfedge_handle(idx);
            _constraints(counter, position++) = 1.0;
            _constraints(counter, n_col) = constraint;
//        if (hehCons.idx() == hehRef.idx())
//            _constraints(counter, n_col) = -1.0 * 0;
//        if (hehCons.idx() == trimesh_.prev_halfedge_handle(hehRef).idx())
//            _constraints(counter, n_col) = -1.0 * trimesh_.calc_sector_angle(hehCons);
//        if (hehCons.idx() == trimesh_.next_halfedge_handle(hehRef).idx())
//            _constraints(counter, n_col) = -1.0 * trimesh_.calc_sector_angle(hehRef);
            counter++;
        }
    }
    gmm::clean(_constraints, 1E-10);
    std::cout << "We have a matrix with " << gmm::mat_nrows(_constraints) <<
              " rows and " << gmm::mat_ncols(_constraints) << " columns:" << std::endl;
    return _constraints;
}

std::vector<int> Crossfield::getConstrainedEdges(const std::vector<int> &constrainedHEdges) {
    std::vector<int> edges;
    for (int i: constrainedHEdges) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i);
        OpenMesh::EdgeHandle eh = trimesh_.edge_handle(heh);
        if (std::find(edges.begin(), edges.end(), eh.idx()) == edges.end()) {
            edges.push_back(eh.idx());
        }
    };
    std::cout << "The constrained Edge vector has the size of:\t" << edges.size()
              << "\nwhile the Halfedge vector has the size of:\t" << constrainedHEdges.size() << std::endl;
    return edges;
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
    double totalArea = getTotalArea(faces), edge_weight = DBL_MAX, sum = DBL_MAX;
    int facesPlusOne = faces.size();
    for (int i: faces) {
        sum = 0;
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        for (auto fh_it = trimesh_.fh_iter(fh); fh_it.is_valid(); ++fh_it) {
            OpenMesh::HalfedgeHandle opposite_heh = trimesh_.opposite_halfedge_handle(*fh_it);
            OpenMesh::SmartFaceHandle opposite_fh;
            // this condition is necessary because an opposite face doesn't exist if the opposite heh is boundary
            if (!trimesh_.is_boundary(opposite_heh)) {
                opposite_fh = trimesh_.face_handle(trimesh_.opposite_halfedge_handle(*fh_it));
                edge_weight = 1;
//                edge_weight =((trimesh_.calc_face_area(fh) / 3) + (trimesh_.calc_face_area(opposite_fh) / 3)) / totalArea;
                auto it = edgeKappa.find(fh_it->idx());
                auto it2 = edgeKappa.find(opposite_heh.idx());
                // TODO check calculations here. not sure if they are correct
                // it doesn't matter that edge_weight doesn't have a value, since this condition is never true,
                // if edge_weight is valueless
                if (it != edgeKappa.end() && it2 != edgeKappa.end()) {
                    sum += (it->second * edge_weight);
                } else if (it2 != edgeKappa.end()) {
                    sum -= (it2->second * edge_weight);
                }

//                else if (it != edgeKappa.end() && it2 != edgeKappa.end()) {
//                    throw std::runtime_error(
//                            "Opposite Halfedges can't both be in edgeKappa!\n");
//                }
            }
        }
        int position = trimesh_.property(pos_matrixA, fh);
        _rhs[position] = sum;
    }
    for (const auto &i: edgeKappa) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i.first);
        OpenMesh::FaceHandle fh = trimesh_.face_handle(heh);
        OpenMesh::FaceHandle opposite_fh = trimesh_.face_handle(trimesh_.opposite_halfedge_handle(heh));
        edge_weight = 1;
//        edge_weight = ((trimesh_.calc_face_area(fh) / 3) + (trimesh_.calc_face_area(opposite_fh) / 3)) / totalArea;
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
    double edge_weight = DBL_MAX, totalArea = getTotalArea(faces);
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
        edge_weight = 1;
//        edge_weight = ((trimesh_.calc_face_area(fh1) / 3) + (trimesh_.calc_face_area(fh2) / 3)) / totalArea;
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

std::map<int, double> Crossfield::getMapHeKappa(const std::vector<int> &faces) {
    std::map<int, double> edgeKappa;
    for (int i: faces) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        for (auto ff_it = trimesh_.ff_iter(fh); ff_it.is_valid(); ++ff_it) {
            getStatusNeigh(fh, *ff_it, edgeKappa);
        }
    }
    if (edgeKappa.empty()) {
        throw std::runtime_error("Something went wrong, edgeKappa Map can't be empty!\n");
    }
    return edgeKappa;
}

void Crossfield::getStatusNeigh(const OpenMesh::FaceHandle fh, const OpenMesh::FaceHandle fh_neigh,
                                std::map<int, double> &edgeKappa) {
    double kappa = DBL_MAX;
    std::pair<int, int> commonEdge = {INT_MAX, INT_MAX};
    // neighbour needs to be in faces
    if (trimesh_.status(fh_neigh).tagged()) {
        // get index of ref edges
        int refEdgeMain = trimesh_.property(referenceHeIdx, fh);
        int refEdgeNeigh = trimesh_.property(referenceHeIdx, fh_neigh);
        // get the common edge between the triangles
        commonEdge = getCommonEdgeBetweenTriangles(fh, fh_neigh, refEdgeMain, refEdgeNeigh);
        kappa = getKappa(refEdgeMain, refEdgeNeigh, commonEdge);
        addKappaHeToMap(commonEdge, edgeKappa);
    }
}

std::pair<int, int>
Crossfield::getCommonEdgeBetweenTriangles(const OpenMesh::FaceHandle fh, const OpenMesh::FaceHandle fh_neigh,
                                          const int &refEdgeMain, const int &refEdgeNeigh) {
    std::pair<int, int> commonEdge = {INT_MAX, INT_MAX};
    for (auto fhe_it_main = trimesh_.fh_iter(fh); fhe_it_main.is_valid(); ++fhe_it_main) {
        for (auto fhe_it_neigh = trimesh_.fh_iter(fh_neigh); fhe_it_neigh.is_valid(); ++fhe_it_neigh) {
            if (trimesh_.opposite_halfedge_handle(*fhe_it_main).idx() == fhe_it_neigh->idx()) {
                commonEdge.first = fhe_it_main->idx();
                commonEdge.second = fhe_it_neigh->idx();
                // as soon as the common edge is found, the loop is obsolete
                goto END_LOOP;
            }
        }
    }
    END_LOOP:
    if (commonEdge.first == INT_MAX || commonEdge.second == INT_MAX) {
        throw std::runtime_error("Something went wrong, commonEdge can't have the value INT_MAX!\n");
    }
    return commonEdge;
}

double
Crossfield::getKappa(const int refEdgeMain, const int refEdgeNeigh, const std::pair<int, int> commonEdge) {
    // there are 9 different scenarios on how reference edges can be placed in two adjacent triangles
    // refEdge shares common edge
    double alpha = DBL_MAX, beta = alpha, kappa = alpha;
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
//    } else {
//        throw std::runtime_error("Something went wrong, there needs to be a ref edge in the main triangle!\n");
    }
    return kappa = alpha + beta;
}

void Crossfield::addKappaHeToMap(const std::pair<int, int> commonEdge, std::map<int, double> &edgeKappa) {
    if (edgeKappa.find(commonEdge.second) == edgeKappa.end()) {
        int edgeIndex = commonEdge.first;
        edgeKappa[edgeIndex] = 0;
    }
}

std::vector<int> Crossfield::getReferenceEdge(const std::vector<int> &constrainedHEdges) {
    // status needs to be released before using, in case it still has saved some stuff from before
    trimesh_.release_halfedge_status();
    trimesh_.release_face_status();
    // request to change the status
    trimesh_.request_halfedge_status();
    trimesh_.request_face_status();
    std::vector<int> faces;
    // is introduced until we fix the rotation problem
    int temp = INT_MAX;

    // assign constraint edges to their faces
    for (int i: constrainedHEdges) {
        setFacesVecWithRefHe(i, temp, faces);
    }

    // assign reference edges in range to faces
    for (int i: heInRange_) {
        setFacesVec(i, temp, faces);

    }
    colorHEdges(constrainedHEdges);
    colorFaces(faces);
    return faces;
}

void Crossfield::setFacesVecWithRefHe(const int i, int &temp, std::vector<int> &faces) {
    OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i);
    if (trimesh_.is_boundary(heh)) {
        heh = trimesh_.opposite_halfedge_handle(heh);
    }
    OpenMesh::FaceHandle fh = trimesh_.face_handle(heh);
    if (!trimesh_.status(heh).tagged() && !trimesh_.status(fh).tagged()) {
        // add reference edge to face
        trimesh_.property(referenceHeIdx, fh) = heh.idx();
        temp = heh.idx();
        // set face and edge as used
        trimesh_.status(heh).set_tagged(true);
        trimesh_.status(fh).set_tagged(true);
        faces.push_back(fh.idx());
    }
}

void Crossfield::setFacesVec(const int i, const int temp, std::vector<int> &faces) {
    OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i);
    if (temp == INT_MAX) {
        throw std::runtime_error("the temp value can't be INT MAX !\n");
    }
    if (trimesh_.is_boundary(heh))
        heh = trimesh_.opposite_halfedge_handle(heh);
    trimesh_.property(heh_color, heh) = 2;
    OpenMesh::FaceHandle fh = trimesh_.face_handle(heh);
    // avoids including faces which only have one halfedge in heInRange
    OpenMesh::HalfedgeHandle nheh = trimesh_.next_halfedge_handle(heh);
    // halfedges or faces which are already tagged, won't be used anymore, no duplicates are allowed
    if ((std::find(heInRange_.begin(), heInRange_.end(), nheh.idx()) != heInRange_.end()) &&
        (!trimesh_.status(heh).tagged() && !trimesh_.status(fh).tagged())) {
        trimesh_.property(referenceHeIdx, fh) = temp;
        // both halfedges need to be tagged, else it is possible opposite faces share an edge
        trimesh_.status(heh).set_tagged(true);
        trimesh_.status(trimesh_.opposite_halfedge_handle(heh)).set_tagged(true);
        trimesh_.status(fh).set_tagged(true);
        faces.push_back(fh.idx());
    }
}


// convert array of edges to faces and fill constraints with them
std::vector<int> Crossfield::getConstraints() {
    std::vector<int> constraints;
    getSelectedVertices(constraints);
    getSelectedEdges(constraints);
    getSelectedHEdges(constraints);
    getSelectedFaces(constraints);
    return constraints;
}

void Crossfield::getSelectedVertices(std::vector<int> &constraints) {
    std::vector<int> selection = MeshSelection::getVertexSelection(&trimesh_);
    // if vertex not empty checks neighbouring vertices. if they are also selected, add edge between
    for (int i: selection) {
        OpenMesh::VertexHandle vh = trimesh_.vertex_handle(i);
        for (int j: selection) {
            OpenMesh::VertexHandle vhj = trimesh_.vertex_handle(j);
            OpenMesh::HalfedgeHandle heh1 = trimesh_.find_halfedge(vhj, vh);
            OpenMesh::HalfedgeHandle heh2 = trimesh_.opposite_halfedge_handle(heh1);
            // avoids duplicates with std::find
            if (std::find(constraints.begin(), constraints.end(), heh1.idx()) ==
                constraints.end()) {
                constraints.push_back(heh1.idx());
            }
            if (std::find(constraints.begin(), constraints.end(), heh2.idx()) ==
                constraints.end()) {
                constraints.push_back(heh2.idx());
            }
        }
    }
}

void Crossfield::getSelectedEdges(std::vector<int> &constraints) {
    std::vector<int> selection = MeshSelection::getEdgeSelection(&trimesh_);
    for (int i: selection) {
        // avoids duplicates with std::find
        OpenMesh::EdgeHandle eh = trimesh_.edge_handle(i);
        OpenMesh::HalfedgeHandle heh1 = trimesh_.halfedge_handle(eh, 0);
        OpenMesh::HalfedgeHandle heh2 = trimesh_.halfedge_handle(eh, 1);
        if (std::find(constraints.begin(), constraints.end(), heh1.idx()) == constraints.end()) {
            constraints.push_back(heh1.idx());
        }
        if (std::find(constraints.begin(), constraints.end(), heh2.idx()) == constraints.end()) {
            constraints.push_back(heh2.idx());
        }
    }
}

void Crossfield::getSelectedHEdges(std::vector<int> &constraints) {
    std::vector<int> selection = MeshSelection::getHalfedgeSelection(&trimesh_);
    for (int i: selection) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i);
        // avoids duplicates with std::find
        if (std::find(constraints.begin(), constraints.end(), heh.idx()) == constraints.end())
            constraints.push_back(heh.idx());
    }
}

void Crossfield::getSelectedFaces(std::vector<int> &constraints) {
    std::vector<int> selection = MeshSelection::getFaceSelection(&trimesh_);
    for (int i: selection) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        for (auto fh_it = trimesh_.fh_iter(fh); fh_it.is_valid(); ++fh_it) {
            // avoids duplicates with std::find
            if (std::find(constraints.begin(), constraints.end(), fh_it->idx()) == constraints.end())
                constraints.push_back(fh_it->idx());
        }
    }
}

void Crossfield::setlocalCoordFrame(const std::vector<int> &faces) {
    for (int i: faces) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(trimesh_.property(referenceHeIdx, fh));
        Point a = trimesh_.calc_edge_vector(heh).normalize();
        double temp1 = OpenMesh::dot(a, {1, 0, 0});
        double temp2 = OpenMesh::dot(a, {0, 1, 0});
        double alpha2 = std::atan2(temp2, temp1);
        trimesh_.property(constraint_angle, fh) = alpha2;
        std::cout << "theta_c (" << i << ") is: " << alpha2 * 180 / M_PI << "degrees." << std::endl;
        Point b = a % trimesh_.calc_face_normal(fh);
        trimesh_.property(xa, fh) = {1, 0, 0};
        trimesh_.property(ya, fh) = {0, 1, 0};
        trimesh_.property(x_vec_field, fh) = a;
        trimesh_.property(y_vec_field, fh) = b;
    }
}

void Crossfield::rotateLocalCoordFrame(const std::vector<int> &faces, const std::vector<double> _x) {
    for (int i: faces) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        int position = trimesh_.property(pos_matrixA, fh);
        Point p_x = trimesh_.property(x_vec_field, fh);
        Point p_y = trimesh_.property(y_vec_field, fh);
        double degrees = std::fmod(_x[position] * 180 / M_PI, 360);
        std::cout << "face " << i << ": with angle: " << degrees << std::endl;
        trimesh_.property(x_vec_field_r, fh) = p_x * cos(_x[position]) + p_y * sin(_x[position]);
        trimesh_.property(y_vec_field_r, fh) = p_x * -sin(_x[position]) + p_y * cos(_x[position]);
    }
}

double Crossfield::getTotalArea(const std::vector<int> &faces) {
    double totalarea = 0.0;
    for (int i: faces)
        totalarea += trimesh_.calc_face_area(trimesh_.face_handle(i));
    return totalarea;
}


void Crossfield::colorFaces(const std::vector<int> &faces) {
    for (auto f_it = trimesh_.faces_begin(); f_it != trimesh_.faces_end(); ++f_it) {
        trimesh_.property(face_color, *f_it) = 0;
    }
    for (int i: faces) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        trimesh_.property(face_color, fh) = 1;
    }
}

void Crossfield::colorHEdges(const std::vector<int> &constrainedEdges) {
    for (int i: heInRange_) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i);
        trimesh_.property(heh_color, heh) = 1;
    }
    for (int i: constrainedEdges) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i);
        trimesh_.property(heh_color, heh) = 2;
    }
}


void Crossfield::removeProperties() {
    trimesh_.remove_property(face_color);
    trimesh_.remove_property(pos_matrixA);
    trimesh_.remove_property(referenceHeIdx);
    trimesh_.remove_property(x_vec_field);
    trimesh_.remove_property(y_vec_field);
    trimesh_.remove_property(x_vec_field_r);
    trimesh_.remove_property(y_vec_field_r);
    trimesh_.remove_property(theta);
    trimesh_.remove_property(heh_color);
}
