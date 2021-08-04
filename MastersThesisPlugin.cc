#include "MastersThesisPlugin.hh"
#include <MeshTools/MeshSelectionT.hh>
#include "DijkstraDistance.hh"
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
        // create mesh
        TriMeshObject *tri_obj = PluginFunctions::triMeshObject(*o_it);
        TriMesh *trimesh = tri_obj->mesh();
        // different creation of mesh: TriMesh *trimesh = PluginFunctions::triMesh(*o_it);

        if (trimesh) {
            DijkstraDistance test{*trimesh};
            test.colorizeEdgeSelection();

            // change layer of display
            PluginFunctions::triMeshObject(*o_it)->meshNode()->drawMode(ACG::SceneGraph::DrawModes::EDGES_COLORED);
            emit updatedObject(o_it->id(), UPDATE_ALL);
        }
    }

}

#if QT_VERSION < 0x050000
Q_EXPORT_PLUGIN2(tetctplugin, TetCTPlugin);
#endif