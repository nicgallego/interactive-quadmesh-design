#ifndef OPENFLIPPER_DIJKSTRADISTANCE_HH
#define OPENFLIPPER_DIJKSTRADISTANCE_HH

#include <QObject>
#include <MeshTools/MeshSelectionT.hh>
#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <OpenFlipper/BasePlugin/BaseInterface.hh>


class DijkstraDistance {
public:

    DijkstraDistance();

    ~DijkstraDistance();

public:

    void getBoundaryVertices(TriMesh *_mesh);

    void calculateDijkstra(TriMesh &_mesh);


private:
    TriMesh *mesh_;


};


#endif //OPENFLIPPER_DIJKSTRADISTANCE_HH
