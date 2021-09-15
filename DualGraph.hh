//
// Created by wuschelbueb on 15.09.21.
//

#ifndef OPENFLIPPER_DUALGRAPH_HH
#define OPENFLIPPER_DUALGRAPH_HH

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <ObjectTypes/PolyMesh/PolyMesh.hh>
#include <ACG/Geometry/Types/PlaneT.hh>
#include <iostream>
#include <vector>
#include <float.h>

class DualGraph {
public:
    // temporary
    typedef OpenMesh::PolyMesh_ArrayKernelT<> MyMesh;
    /*
     * is the same as:
     * DualGraph(TriMesh &trimesh) {
     *      trimesh_ = trimesh;
     * }
     */
    DualGraph(TriMesh &trimesh)
            : trimesh_{trimesh} {
    }

    ~DualGraph() {}

public:

    void getCircumCenter();

private:
    TriMesh &trimesh_;
    PolyMesh dualGraph_;
};


#endif //OPENFLIPPER_DUALGRAPH_HH
