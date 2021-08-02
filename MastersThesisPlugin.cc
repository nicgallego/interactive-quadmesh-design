#include "MastersThesisPlugin.hh"
#include <MeshTools/MeshSelectionT.hh>
#include "OpenFlipper/BasePlugin/PluginFunctions.hh"
#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>
#include <iostream>
#include <vector>

typedef OpenMesh::TriMesh_ArrayKernelT<> MyMesh;

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
                OpenMesh::EPropHandleT<float> edge_length;
                TriMesh::Color green = {0.0,255.0,255.0,0};
                std::cout << "neighbouring vertices of " << vh << ": " << trimesh->point(vh) << " are:\n";
                // find adjacent edges
                for (TriMesh::VertexEdgeIter ve_it = trimesh->ve_iter(vh); ve_it.is_valid(); ++ve_it) {
                    trimesh->add_property(edge_length, "edge_length");
                    trimesh->property(edge_length,*ve_it) = trimesh->calc_edge_length(*ve_it);
                    trimesh->request_edge_colors();
                    trimesh->set_color(*ve_it, green);
                    float LL = trimesh->property(edge_length, *ve_it);
                    std::cout << *ve_it << " length is: " <<  LL << " ||\t";
                }
                trimesh->remove_property(edge_length);
                std::cout << "\n=============\n";
            }


        }
    }

}

#if QT_VERSION < 0x050000
Q_EXPORT_PLUGIN2(tetctplugin, TetCTPlugin);
#endif