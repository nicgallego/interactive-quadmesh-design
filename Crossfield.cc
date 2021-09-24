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
    getBaryCenter(barycenters);
}

void Crossfield::getBaryCenter(std::vector<Point> &barycenters) {
    int numberVerticesInTriangle = 3;
    std::vector<OpenMesh::VertexHandle> vhs;
    std::vector<int> constraints;
    getConstraints(constraints);
    //use faces to assign edges and barcenters to it
    for (int i: constraints) {
        std::vector<Point> triangle;
        triangle.reserve(numberVerticesInTriangle);
        OpenMesh::EdgeHandle eh = trimesh_.edge_handle(i);
        // get random halfedge
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(eh, 1);
        if (trimesh_.is_boundary(heh))
            heh = trimesh_.opposite_halfedge_handle(heh);
        OpenMesh::FaceHandle fh = trimesh_.face_handle(heh);

        for (auto fv_it = trimesh_.fv_iter(fh); fv_it.is_valid(); ++fv_it) {
            triangle.push_back(trimesh_.point(*fv_it));
        }
        Point bCenter = (triangle[0] + triangle[1] + triangle[2]) / 3.0;
        barycenters.push_back(bCenter);
        trimesh_.property(inUse, eh) = true;
        trimesh_.property(associated_edge_to_face, fh) = eh;
        trimesh_.property(barycenter, fh) = bCenter;
    }
}

// convert array of halfedges to faces and fill constraints with them
void Crossfield::getConstraints(std::vector<int> &constraints) {

    std::cout << "hello wolrd\n";
}



