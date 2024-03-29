#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QFileDialog>
#include <QMainWindow>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

namespace Ui {
class MainWindow;
}

using namespace OpenMesh;
using namespace OpenMesh::Attributes;

struct MyTraits : public OpenMesh::DefaultTraits
{
    // use vertex normals and vertex colors
    VertexAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color );
    // store the previous halfedge
    HalfedgeAttributes( OpenMesh::Attributes::PrevHalfedge );
    // use face normals face colors
    FaceAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color );
    EdgeAttributes( OpenMesh::Attributes::Color );
    // vertex thickness
    VertexTraits{float thickness;};
    // edge thickness
    EdgeTraits{float thickness;};
};
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    // les 4 fonctions à compléter
    void showSelections(MyMesh* _mesh);
    void showSelectionsNeighborhood(MyMesh* _mesh);
    void showPath(MyMesh* _mesh, int v1, int v2, bool vertex);
    void showBorder(MyMesh* _mesh);
    void displayMesh(MyMesh *_mesh);
    void resetAllColorsAndThickness(MyMesh* _mesh);

private slots:

    void on_pushButton_chargement_clicked();
    void on_pushButton_vertexMoins_clicked();
    void on_pushButton_vertexPlus_clicked();
    void on_pushButton_edgeMoins_clicked();
    void on_pushButton_edgePlus_clicked();
    void on_pushButton_faceMoins_clicked();
    void on_pushButton_facePlus_clicked();
    void on_pushButton_afficherChemin_clicked();
    void on_pushButton_aretes_clicked();
    void on_pushButton_voisinage_clicked();
    void on_pushButton_bordure_clicked();

private:

    void calculEdgeWeight(MyMesh* _mesh, float dist[]);
    void dijkstraVertex(MyMesh* _mesh, int u, int v);
    void dijkstraEdge(MyMesh* _mesh, int u, int v);

    bool modevoisinage;

    MyMesh mesh;

    int vertexSelection;
    int edgeSelection;
    int faceSelection;

    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
