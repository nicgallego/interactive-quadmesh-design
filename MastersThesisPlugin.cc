#include "MastersThesisPlugin.hh"
#include <MeshTools/MeshSelectionT.hh>
#include "OpenFlipper/BasePlugin/PluginFunctions.hh"
#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <iostream>
#include <vector>

void MastersThesisPlugin::initializePlugin() {
    tool_ = new MastersThesisToolbar();
    QSize size(300, 300);
    tool_->resize(size);

    connect(tool_->get_mesh_boundaries, SIGNAL(clicked()), this, SLOT(slot_get_boundary()));

    emit addToolbox(tr("MastersThesis"), tool_);

}

void MastersThesisPlugin::pluginsInitialized() {}

void MastersThesisPlugin::slot_get_boundary() {

    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH);
         o_it != PluginFunctions::objectsEnd(); ++o_it) {
        // get mesh
        TriMesh *trimesh = PluginFunctions::triMesh(*o_it);

        if (trimesh) {
            MeshSelection::selectBoundaryVertices(trimesh);
            // get boundary vertices in vector
            std::vector<int> boundaryVertices = MeshSelection::getVertexSelection(trimesh);
            for (int i: boundaryVertices) {
//                std::cout << "i: " << i << std::endl;
                // add handle to vertices so you can work with them
                OpenMesh::VertexHandle vh = trimesh->vertex_handle(i);
                std::cout << "neighbouring vertices of " << vh << ": " << trimesh->point(vh) << " are:\n";
                // find neighbouring vertices
                for(auto vv_it = trimesh->vv_iter(vh); vv_it.is_valid(); vv_it++){
                    std::cout  << *vv_it << ": " << trimesh->point(*vv_it) << " ||\t";
                }
                std::cout << "\n=============\n";
            }


        }
    }

}

#if QT_VERSION < 0x050000
Q_EXPORT_PLUGIN2(tetctplugin, TetCTPlugin);
#endif