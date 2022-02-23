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
    setlocalCoordFrame(faces);
    std::map<int, double> edgeKappa = getMapHeKappa(faces);
    gmm::col_matrix<gmm::wsvector<double>> _A = getMatrixA(faces, edgeKappa);
    gmm::row_matrix<gmm::wsvector<double>> _constraints = getConstraintMatrix(edgeKappa, constrainedHalfEdges, faces);
    std::vector<double> _x(edgeKappa.size() + faces.size(), 0.0);
    std::vector<double> _rhs = getRHS(edgeKappa, faces);
    std::vector<int> _idx_to_round = getIdxToRound(edgeKappa, faces);
    std::vector<double> _rhsOld = _rhs;
    std::vector<int> _rhsNew;
    std::vector<int> _c_elim;

    std::cout << "A (before): " << _A << std::endl;

    for (std::size_t i = 0, max = _rhs.size(); i != max; ++i) {
        std::cout << "rhs position (" << i << ") with value (before calc): " << _rhs[i] << std::endl;
    }
    std::cout << std::endl;
    for (std::size_t i = 0, max = _idx_to_round.size(); i != max; ++i) {
        std::cout << "idx to round position (" << i << ") with value (before calc): " << _idx_to_round[i] << std::endl;
    }
    std::cout << "Constraint matrix: " << gmm::mat_nrows(_constraints) <<
              " rows and " << gmm::mat_ncols(_constraints) << " columns:" << std::endl;
    std::cout << "constraints: " << _constraints << std::endl << std::endl;

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
    csolver.solve_const(_constraints, _A, _x, _rhs, _idx_to_round);

    rotateLocalCoordFrame(faces, _x);

//    for (std::size_t i = 0, max = _x.size(); i != max; ++i) {
//        std::cout << "_x position (" << i << ") with value: " << _x[i] << std::endl;
//    }
//
//    std::cout << "A (after): " << _A << std::endl;
//    for (std::size_t i = 0, max = _rhs.size(); i != max; ++i) {
//        std::cout << "_rhs position (" << i << ") with value (after calc): " << _rhs[i] << std::endl;
//    }
//
//    std::cout << "constraints after: " << _constraints << std::endl;
    double energy_before = getEnergy(edgeKappa, faces, _rhsOld, _x);
    double energy_after = getEnergy(edgeKappa, faces, _rhs, _x);
    std::cout << "The energy before the smoothing is:\t" << energy_before << std::endl
              << "The energy after the smoothing is:\t" << energy_after << std::endl;

    // test examples
    const double tol = 1e-6;
    int n = gmm::mat_nrows(_A);
    CVectorType y(n), b(n);
    gmm::copy(_x, y);
    gmm::copy(_rhs, b);
    double ea = computeEnergy(_A, y, b);
    bool eb_ok = std::abs(ea - energy_after) < tol;
//    assert(eb_ok);
    if (!eb_ok) {
        std::cerr << "the energies do not coincide: " << energy_after
                  << " x^T A x + b^t x = " << ea << std::endl;
    }
}

double Crossfield::getEnergy(const std::map<int, double> &edgeKappa, const std::vector<int> &faces,
                             const std::vector<double> &_rhs, const std::vector<double> &_x) {
    auto pos_matrixA = OpenMesh::FProp<int>(trimesh_, "pos_matrixA");
    int pj_start = faces.size();
    double m_pi = M_PI;
    double sum = 0.0;
    for (auto i: edgeKappa) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i.first);
        OpenMesh::FaceHandle fh = trimesh_.face_handle(heh);
        OpenMesh::HalfedgeHandle oheh = trimesh_.opposite_halfedge_handle(heh);
        OpenMesh::FaceHandle ofh = trimesh_.face_handle(oheh);
        int position = pos_matrixA[fh];
        int opposite_position = pos_matrixA[ofh];
//        std::cout << "theta_i: " << _x[position] << std::endl << "kappa: " << _rhs[pj_start] << std::endl
//                  << "periodjump: "
//                  << _x[pj_start] << std::endl << "theta_j: " << _x[opposite_position] << std::endl;
        sum += pow((_x[position] + _rhs[pj_start] + (M_PI / 2. * _x[pj_start]) - _x[opposite_position]), 2);
        pj_start++;
    }
    return sum;
}

gmm::row_matrix<gmm::wsvector<double>>
Crossfield::getConstraintMatrix(const std::map<int, double> &edgeKappa, const std::vector<int> &constraintHalfEdges,
                                const std::vector<int> &faces) {
    auto constraint_angle = OpenMesh::FProp<double>(trimesh_, "constraint_angle");
    auto pos_matrixA = OpenMesh::FProp<int>(trimesh_, "pos_matrixA");
    auto referenceHeIdx = OpenMesh::FProp<int>(trimesh_, "referenceHeIdx");
    //cNplusOne = angle between local coordinate x-axis and constraint
    int cNplusOne = 1, counter = 0, pj_start = faces.size();
    // this results in have the rows of the constrainedhalfedges
    std::vector<int> constraintEdges = getConstraintEdges(
            constraintHalfEdges), noOriginConst = getNumberOriginConstraints(faces);
    int n_row = constraintEdges.size() + noOriginConst.size(), n_col = edgeKappa.size() + faces.size();
    gmm::row_matrix<gmm::wsvector<double>> _constraints(n_row, n_col + cNplusOne);
    for (int i: constraintEdges) {
        // what if edge is between faces, do both faces get a 1 in their respective position in the row of the matrix
        OpenMesh::EdgeHandle eh = trimesh_.edge_handle(i);
        OpenMesh::HalfedgeHandle hehConst = trimesh_.halfedge_handle(eh, 0);
        if (trimesh_.is_boundary(hehConst)) {
            hehConst = trimesh_.opposite_halfedge_handle(hehConst);
        }
        OpenMesh::FaceHandle fh = trimesh_.face_handle(hehConst);
        double constraint = constraint_angle[fh];
        int position = pos_matrixA[fh];
        int idx = referenceHeIdx[fh];
        OpenMesh::HalfedgeHandle hehRef = trimesh_.halfedge_handle(idx);
        _constraints(counter, position) = 1.0;
        if (hehConst.idx() == hehRef.idx())
            _constraints(counter, n_col) = 1.0 * constraint;
        if (hehConst.idx() == trimesh_.prev_halfedge_handle(hehRef).idx())
            _constraints(counter, n_col) = 1.0 * trimesh_.calc_sector_angle(hehConst) + constraint;
        if (hehConst.idx() == trimesh_.next_halfedge_handle(hehRef).idx())
            _constraints(counter, n_col) = 1.0 * trimesh_.calc_sector_angle(hehRef) + constraint;
        counter++;
    }
    dualSpanningTreeConstraint(edgeKappa, counter, pj_start, noOriginConst, _constraints);
    if (counter != n_row) {
        std::ostringstream oss;
        oss << "getConstraintMatrix: amount of rows has to be the: " << n_row << " and not: " << counter << "\n";
        throw std::runtime_error(oss.str());
    }
    gmm::clean(_constraints, 1E-10);
    return _constraints;
}

void Crossfield::dualSpanningTreeConstraint(const std::map<int, double> &edgeKappa, int &counter, const int pj_start,
                                            const std::vector<int> &noOriginConst,
                                            gmm::row_matrix<gmm::wsvector<double>> &_constraints) {
    auto origin_constraint = OpenMesh::FProp<int>(trimesh_, "origin_constraint");
    for (int i: noOriginConst) {
        int iteration = 0;
        for (auto j: edgeKappa) {
            OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(j.first);
            OpenMesh::FaceHandle fh1 = trimesh_.face_handle(heh);
            OpenMesh::FaceHandle fh2 = trimesh_.face_handle(trimesh_.opposite_halfedge_handle(heh));
            int fh1_origin = origin_constraint[fh1];
            int fh2_origin = origin_constraint[fh2];
//            std::cout << "Halfedge: " << heh.idx() << " with faces: " << fh1.idx() << " (originconst: "
//                      << origin_constraint[fh1] << ") and " << fh2.idx() << " (originconst: " << origin_constraint[fh2]
//                      << ")\n";
            if (fh1_origin == fh2_origin && fh1_origin == i) {
                _constraints(counter, pj_start + iteration++) = 1;
            } else {
                _constraints(counter, pj_start + iteration++) = 0;
            }
        }
        counter++;
    }
}

std::vector<int> Crossfield::getNumberOriginConstraints(const std::vector<int> &faces) {
    auto origin_constraint = OpenMesh::FProp<int>(trimesh_, "origin_constraint");
    std::vector<int> noOrigins;
    for (int i: faces) {
        auto fh = trimesh_.face_handle(i);
        if (std::find(noOrigins.begin(), noOrigins.end(), origin_constraint[fh]) == noOrigins.end()) {
            noOrigins.push_back(origin_constraint[fh]);
        }
    }
    return noOrigins;
}

std::vector<int> Crossfield::getConstraintEdges(const std::vector<int> &constrainedHEdges) {
    std::vector<int> edges;
    for (int i: constrainedHEdges) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i);
        OpenMesh::EdgeHandle eh = trimesh_.edge_handle(heh);
        if (std::find(edges.begin(), edges.end(), eh.idx()) == edges.end()) {
            edges.push_back(eh.idx());
        }
    };
    return edges;
}

std::vector<int> Crossfield::getIdxToRound(const std::map<int, double> &edgeKappa, const std::vector<int> &faces) {
    std::vector<int> _idx_to_round;
    int pj_start = faces.size();
    for (const auto &i: edgeKappa) {
        _idx_to_round.push_back(pj_start++);
    }
    return _idx_to_round;
}

std::vector<double> Crossfield::getRHS(const std::map<int, double> &edgeKappa, const std::vector<int> &faces) {
    std::vector<double> _rhs(faces.size() + edgeKappa.size());
    double totalArea = getTotalArea(faces);
    int facesPlusOne = faces.size();

    getRhsFirstHalf(totalArea, faces, _rhs, edgeKappa);
    getRhsSecondHalf(totalArea, _rhs, edgeKappa, facesPlusOne);

    // scalar multiplication with vector, b*-1 = -b
    gmm::scale(_rhs, -1.0);
    return _rhs;
}

void Crossfield::getRhsFirstHalf(double const totalArea, const std::vector<int> &faces, std::vector<double> &_rhs,
                                 const std::map<int, double> &edgeKappa) {
    for (int i: faces) {
        getSum(totalArea, i, _rhs, edgeKappa);
    }
}

void Crossfield::getSum(double const totalArea, const int i, std::vector<double> &_rhs,
                        const std::map<int, double> &edgeKappa) {
    auto pos_matrixA = OpenMesh::FProp<int>(trimesh_, "pos_matrixA");
    double sum = 0;
    OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
    for (TriMesh::FaceHalfedgeIter fh_it = trimesh_.fh_iter(fh); fh_it.is_valid(); ++fh_it) {
        // this condition is necessary because an opposite face doesn't exist if the opposite heh is boundary
        sum += setSum(fh, *fh_it, edgeKappa, totalArea);
    }
    int position = pos_matrixA[fh];
    _rhs[position] = sum;
}

double
Crossfield::setSum(OpenMesh::FaceHandle fh, OpenMesh::HalfedgeHandle fh_it, const std::map<int, double> &edgeKappa,
                   double const totalArea) {
    double sum = 0;
    OpenMesh::HalfedgeHandle opposite_heh = trimesh_.opposite_halfedge_handle(fh_it);
    if (!trimesh_.is_boundary(opposite_heh)) {
        OpenMesh::FaceHandle opposite_fh = trimesh_.face_handle(trimesh_.opposite_halfedge_handle(fh_it));
        double edge_weight = 1;
        edge_weight = ((trimesh_.calc_face_area(fh) / 3) + (trimesh_.calc_face_area(opposite_fh) / 3)) / totalArea;
        auto it = edgeKappa.find(fh_it.idx());
        auto it2 = edgeKappa.find(opposite_heh.idx());
        // it doesn't matter that edge_weight doesn't have a value, since this condition is never true,
        // if edge_weight is valueless
        if (it != edgeKappa.end()) {
            sum += (it->second * edge_weight);
//            std::cout << "Kappa value main triangle: " << it->second << "\nsum main:\t" << sum << std::endl;
        }
        if (it2 != edgeKappa.end()) {
            sum -= (it2->second * edge_weight);
//            std::cout << "Kappa value neighbour triangle: " << it2->second << "\nsum neigh:\t" << sum << std::endl;
//        } else if (it != edgeKappa.end() && it2 != edgeKappa.end()) {
//            throw std::runtime_error(
//                    "Opposite Halfedges can't both be in edgeKappa!\n");
        }
    }
    return sum;
}

void Crossfield::getRhsSecondHalf(double const totalArea, std::vector<double> &_rhs,
                                  const std::map<int, double> &edgeKappa, const int facesPlusOne) {
    int counter = facesPlusOne;
    for (const auto &i: edgeKappa) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i.first);
        OpenMesh::FaceHandle fh = trimesh_.face_handle(heh);
        OpenMesh::FaceHandle opposite_fh = trimesh_.face_handle(trimesh_.opposite_halfedge_handle(heh));
        double edge_weight = 1;
        edge_weight = ((trimesh_.calc_face_area(fh) / 3) + (trimesh_.calc_face_area(opposite_fh) / 3)) / totalArea;
        _rhs[counter++] = (i.second * edge_weight);
    }
}

gmm::col_matrix<gmm::wsvector<double>>
Crossfield::getMatrixA(const std::vector<int> &faces, const std::map<int, double> &edgeKappa) {
    trimesh_.release_halfedge_status();
    trimesh_.release_face_status();
    // request to change the status
    trimesh_.request_halfedge_status();
    trimesh_.request_face_status();
    auto pos_matrixA = OpenMesh::FProp<int>(trimesh_, "pos_matrixA");
    int counter = 0, iteration = 0, pj_start = faces.size(), n = edgeKappa.size() + faces.size();
    double edge_weight = DBL_MAX, totalArea = getTotalArea(faces);
    gmm::col_matrix<gmm::wsvector<double>> _A(n, n);
    // indexes the faces from 0 to n
    // this is needed in case a face index is higher than the matrix size
    for (const auto &i: edgeKappa) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i.first);
        OpenMesh::FaceHandle fh1 = trimesh_.face_handle(heh);
        OpenMesh::FaceHandle fh2 = trimesh_.face_handle(trimesh_.opposite_halfedge_handle(heh));
        setPositionForFace(fh1, counter);
        setPositionForFace(fh2, counter);
    }
    // fills up sparse column matrix _A
    for (auto i: edgeKappa) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i.first);
        OpenMesh::FaceHandle fh1 = trimesh_.face_handle(heh);
        OpenMesh::FaceHandle fh2 = trimesh_.face_handle(trimesh_.opposite_halfedge_handle(heh));
        edge_weight = ((trimesh_.calc_face_area(fh1) / 3) + (trimesh_.calc_face_area(fh2) / 3)) / totalArea;
        int pos_i = pos_matrixA[fh1];
        int pos_j = pos_matrixA[fh2];
        // |  2*w_ij   -2*w_ij     pi*w_ij |
        // | -2*w_ij    2*w_ij    -pi*w_ij |
        // | pi*w_ij  -pi*w_ij pi^2/2*w_ij |
        // added epsilon 0.1 to make the matrix positive definite, careful!
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

void Crossfield::setPositionForFace(OpenMesh::FaceHandle fh, int &counter) {
    auto pos_matrixA = OpenMesh::FProp<int>(trimesh_, "pos_matrixA");
    auto referenceHeIdx = OpenMesh::FProp<int>(trimesh_, "referenceHeIdx");
    if (!trimesh_.status(fh).tagged()) {
        pos_matrixA[fh] = counter;
        trimesh_.status(fh).set_tagged(true);
        counter++;
    }
}

std::map<int, double> Crossfield::getMapHeKappa(const std::vector<int> &faces) {
    std::map<int, double> edgeKappa;
    for (int i: faces) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        for (TriMesh::FaceFaceIter ff_it = trimesh_.ff_iter(fh); ff_it.is_valid(); ++ff_it) {
            getStatusNeigh(fh, *ff_it, edgeKappa);
        }
    }
    if (edgeKappa.empty()) {
        throw std::runtime_error("getMapHeKappa: edgeKappa Map can't be empty!\n");
    }
    return edgeKappa;
}

void Crossfield::getStatusNeigh(const OpenMesh::FaceHandle fh, const OpenMesh::FaceHandle fh_neigh,
                                std::map<int, double> &edgeKappa) {
    auto referenceHeIdx = OpenMesh::FProp<int>(trimesh_, "referenceHeIdx");
    auto constraint_angle = OpenMesh::FProp<double>(trimesh_, "constraint_angle");
    double kappa = DBL_MAX;
    std::pair<int, int> commonEdge = {INT_MAX, INT_MAX};
    // neighbour needs to be in faces
    if (trimesh_.status(fh_neigh).tagged()) {
        // get index of ref edges
        int refEdgeMain = referenceHeIdx[fh];
        int refEdgeNeigh = referenceHeIdx[fh_neigh];
        // get the common edge between the triangles
        commonEdge = getCommonEdgeBetweenTriangles(fh, fh_neigh, refEdgeMain, refEdgeNeigh);
        kappa = getKappa(refEdgeMain, refEdgeNeigh, commonEdge);

        // example test kappa
        bool kappa_ok = testKappa(refEdgeMain, refEdgeNeigh, commonEdge);
//        assert(kappa_ok);
        addKappaHeToMap(commonEdge, kappa, edgeKappa);
    }
}

std::pair<int, int>
Crossfield::getCommonEdgeBetweenTriangles(const OpenMesh::FaceHandle fh, const OpenMesh::FaceHandle fh_neigh,
                                          const int &refEdgeMain, const int &refEdgeNeigh) {
    std::pair<int, int> commonEdge = {INT_MAX, INT_MAX};
    for (TriMesh::FaceHalfedgeIter fhe_it_main = trimesh_.fh_iter(fh); fhe_it_main.is_valid(); ++fhe_it_main) {
        for (TriMesh::FaceHalfedgeIter fhe_it_neigh = trimesh_.fh_iter(
                fh_neigh); fhe_it_neigh.is_valid(); ++fhe_it_neigh) {
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
        throw std::runtime_error("getCommonEdgeBetweenTriangles: commonEdge can't have the value INT_MAX!\n");
    }
    return commonEdge;
}

double
Crossfield::getKappa(const int refEdgeMain, const int refEdgeNeigh, const std::pair<int, int> commonEdge) {
    // there are 9 different scenarios on how reference edges can be placed in two adjacent triangles
    // refEdge shares common edge
    double alpha = 0.0, beta = alpha, kappa = alpha;
    if (refEdgeMain == commonEdge.first) {
        alpha = 0.0;
        if (refEdgeNeigh == commonEdge.second)
            beta = 0.0;
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
            beta = 0.0;
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
            beta = 0.0;
        if (refEdgeNeigh ==
            trimesh_.prev_halfedge_handle(trimesh_.halfedge_handle(commonEdge.second)).idx())
            beta = trimesh_.calc_sector_angle(trimesh_.halfedge_handle(refEdgeNeigh));
        if (refEdgeNeigh ==
            trimesh_.next_halfedge_handle(trimesh_.halfedge_handle(commonEdge.second)).idx())
            beta = M_PI - trimesh_.calc_sector_angle(trimesh_.halfedge_handle(commonEdge.second));
//    } else {
//        throw std::runtime_error("Something went wrong, there needs to be a ref edge in the main triangle!\n");
    }
    kappa = alpha + beta;
//    if (kappa <= 0 || kappa >= 2 * M_PI) {
//        throw std::runtime_error(
//                std::string("getKappa: kappa can't be bigger than ") + std::to_string(kappa * 180 / M_PI));
//    }
    return kappa;
}

void Crossfield::addKappaHeToMap(const std::pair<int, int> commonEdge, const double kappa,
                                 std::map<int, double> &edgeKappa) {
    if (edgeKappa.find(commonEdge.second) == edgeKappa.end()) {
        int RefHEdgeIndex = commonEdge.first;
        std::cout << "RefHEdgeIndex: " << RefHEdgeIndex << " kappa: " << kappa * 180 / M_PI << std::endl << std::endl;
        edgeKappa[RefHEdgeIndex] = kappa;        
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
    auto referenceHeIdx = OpenMesh::FProp<int>(trimesh_, "referenceHeIdx");
    OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i);
    if (trimesh_.is_boundary(heh)) {
        heh = trimesh_.opposite_halfedge_handle(heh);
    }
    OpenMesh::FaceHandle fh = trimesh_.face_handle(heh);
    if (!trimesh_.status(heh).tagged() && !trimesh_.status(fh).tagged()) {
        // add reference edge to face
        referenceHeIdx[fh] = heh.idx();
        temp = heh.idx();
        // set face and edge as used
        trimesh_.status(heh).set_tagged(true);
        trimesh_.status(fh).set_tagged(true);
        faces.push_back(fh.idx());
    }
}

void Crossfield::setFacesVec(const int i, int &temp, std::vector<int> &faces) {
    auto referenceHeIdx = OpenMesh::FProp<int>(trimesh_, "referenceHeIdx");
    auto heh_color = OpenMesh::HProp<int>(trimesh_, "heh_color");
    OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i);
    if (trimesh_.is_boundary(heh))
        heh = trimesh_.opposite_halfedge_handle(heh);
    heh_color(heh) = 2;
    OpenMesh::FaceHandle fh = trimesh_.face_handle(heh);
    // avoids including faces which only have one halfedge in heInRange
    OpenMesh::HalfedgeHandle nheh = trimesh_.next_halfedge_handle(heh);
    // halfedges or faces which are already tagged, won't be used anymore, no duplicates are allowed
    if ((std::find(heInRange_.begin(), heInRange_.end(), nheh.idx()) != heInRange_.end()) &&
        (!trimesh_.status(heh).tagged() && !trimesh_.status(fh).tagged())) {
        referenceHeIdx[fh] = heh.idx();
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
        for (TriMesh::FaceHalfedgeIter fh_it = trimesh_.fh_iter(fh); fh_it.is_valid(); ++fh_it) {
            // avoids duplicates with std::find
            if (std::find(constraints.begin(), constraints.end(), fh_it->idx()) == constraints.end())
                constraints.push_back(fh_it->idx());
        }
    }
}

void Crossfield::setlocalCoordFrame(const std::vector<int> &faces) {
    auto referenceHeIdx = OpenMesh::FProp<int>(trimesh_, "referenceHeIdx");
    auto constraint_angle = OpenMesh::FProp<double>(trimesh_, "constraint_angle");
    auto xa = OpenMesh::FProp<Point>(trimesh_, "xa");
    auto ya = OpenMesh::FProp<Point>(trimesh_, "ya");
    auto x_vec_field = OpenMesh::FProp<Point>(trimesh_, "x_vec_field");
    auto y_vec_field = OpenMesh::FProp<Point>(trimesh_, "y_vec_field");
    auto x_vec_field_m = OpenMesh::FProp<Point>(trimesh_, "x_vec_field_m");
    auto y_vec_field_m = OpenMesh::FProp<Point>(trimesh_, "y_vec_field_m");
    for (int i: faces) {
        double alpha = 0.0;
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(referenceHeIdx[fh]);
        Point a = trimesh_.calc_edge_vector(heh).normalize();
        double xComponent = OpenMesh::dot(a, {1, 0, 0});
        double yComponent = OpenMesh::dot(a, {0, 1, 0});
        alpha = std::atan2(yComponent, xComponent);
//        if (alpha < 0) {
//            alpha = 2 * M_PI + alpha;
//        }
        constraint_angle[fh] = alpha;
        std::cout << "theta_c (" << i << ") is: " << alpha * 180 / M_PI << " degrees." << std::endl;
        Point b = a % trimesh_.calc_face_normal(fh);
        xa[fh] = {1, 0, 0};
        ya[fh] = {0, 1, 0};
        x_vec_field[fh] = a;
        x_vec_field_m[fh] = -a;
        y_vec_field[fh] = b;
        y_vec_field_m[fh] = -b;
    }
    std::cout << std::endl;
}

void Crossfield::rotateLocalCoordFrame(const std::vector<int> &faces, const std::vector<double> _x) {
    auto pos_matrixA = OpenMesh::FProp<int>(trimesh_, "pos_matrixA");
    auto x_vec_field = OpenMesh::FProp<Point>(trimesh_, "x_vec_field");
    auto y_vec_field = OpenMesh::FProp<Point>(trimesh_, "y_vec_field");
    auto x_vec_field_rot = OpenMesh::FProp<Point>(trimesh_, "x_vec_field_rot");
    auto y_vec_field_rot = OpenMesh::FProp<Point>(trimesh_, "y_vec_field_rot");
    for (int i: faces) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        int position = pos_matrixA[fh];
        Point p_x = x_vec_field[fh];
        Point p_y = y_vec_field[fh];
        double radians = std::fmod(_x[position], 2 * M_PI);
        std::cout << "face (" << i << ") with (" << _x[position] << "): with angle: " << radians << "(rad) and "
                  << radians * 180 / M_PI << "(deg)"
                  << std::endl;
        x_vec_field_rot[fh] = p_x * cos(radians) - p_y * sin(radians);
        y_vec_field_rot[fh] = p_x * sin(radians) + p_y * cos(radians);
    }
}

double Crossfield::getTotalArea(const std::vector<int> &faces) {
    double totalarea = 0.0;
    for (int i: faces)
        totalarea += trimesh_.calc_face_area(trimesh_.face_handle(i));
    return totalarea;
}


void Crossfield::colorFaces(const std::vector<int> &faces) {
    auto face_color = OpenMesh::FProp<int>(trimesh_, "face_color");
    for (auto f_it = trimesh_.faces_begin(); f_it != trimesh_.faces_end(); ++f_it) {
        face_color[*f_it] = 0;
    }
    for (int i: faces) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        face_color[fh] = 1;
    }
}

void Crossfield::colorHEdges(const std::vector<int> &constrainedEdges) {
    auto heh_color = OpenMesh::HProp<int>(trimesh_, "heh_color");
    for (int i: heInRange_) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i);
        heh_color(heh) = 1;
    }
    for (int i: constrainedEdges) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i);
        heh_color(heh) = 2;
    }
}
