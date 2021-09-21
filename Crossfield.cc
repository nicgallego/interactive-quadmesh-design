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
    getCircumCenter();
}

void Crossfield::getCircumCenter() {
    std::vector<int> constraints = getConstraints();
    for (int i: constraints) {
        OpenMesh::EdgeHandle eh = trimesh_.edge_handle(i);
        // get random halfedge
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(eh, 1);
        if (trimesh_.is_boundary(heh))
            heh = trimesh_.opposite_halfedge_handle(heh);
        OpenMesh::FaceHandle fh = trimesh_.face_handle(heh);
        for (auto fv_it = trimesh_.fv_iter(fh); fv_it.is_valid(); ++fv_it) {
            std::cout << "This is vertex: " << *fv_it << " with " << std::endl;
        }
    }
}

std::vector<int> Crossfield::getConstraints() {
    std::vector<int> selectedVertices = MeshSelection::getVertexSelection(&trimesh_);
    std::vector<int> selectedEdges = MeshSelection::getEdgeSelection(&trimesh_);
    std::vector<int> selectedHEdges = MeshSelection::getHalfedgeSelection(&trimesh_);
    std::vector<int> selectedFaces = MeshSelection::getFaceSelection(&trimesh_);

    return selectedEdges;

}



