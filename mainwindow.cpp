#include "mainwindow.h"
#include "ui_mainwindow.h"

/* **** début de la partie à compléter **** */

/**
 * Affiche la sélection des point, arête et face
 * @brief MainWindow::showSelections
 * @param _mesh
 */
void MainWindow::showSelections(MyMesh* _mesh)
{
    // on réinitialise les couleurs de tout le maillage
    resetAllColorsAndThickness(_mesh);

    //Affiche le sommet sélectionné
    if(vertexSelection >= 0 && vertexSelection < _mesh->n_vertices()){
        VertexHandle vh =_mesh->vertex_handle(vertexSelection); // VertexHandle d'indice i
        _mesh->set_color(vh, MyMesh::Color(255, 0, 0));
        _mesh->data(vh).thickness = 8;
    }

    //Affiche l'arête sélectionné
    if(edgeSelection >= 0 && edgeSelection < _mesh->n_edges()){
        EdgeHandle eh = _mesh->edge_handle(edgeSelection); // EdgeHandle d'indice i
        _mesh->set_color(eh, MyMesh::Color(0, 255, 0));
        _mesh->data(eh).thickness = 4;

        HalfedgeHandle heh = _mesh->halfedge_handle(eh, 0);
        VertexHandle vh1 = _mesh->to_vertex_handle(heh);
        _mesh->set_color(vh1, MyMesh::Color(0, 255, 0));
        _mesh->data(vh1).thickness = 8;

        VertexHandle vh2 = _mesh->from_vertex_handle(heh);
        _mesh->set_color(vh2, MyMesh::Color(0, 255, 0));
        _mesh->data(vh2).thickness = 8;
    }

    //Affiche la face sélectionné
    if(faceSelection >= 0 && faceSelection < _mesh->n_faces()) {
        FaceHandle fh = _mesh->face_handle(faceSelection); // FaceHandle d'indice i
        _mesh->set_color(fh, MyMesh::Color(0, 0, 255));

        HalfedgeHandle heh = _mesh->halfedge_handle(fh);
        for(int i = 0; i < 3; ++i){
            EdgeHandle eh = _mesh->edge_handle(heh);
            _mesh->set_color(eh, MyMesh::Color(0, 0, 180));
            _mesh->data(eh).thickness = 4;

            VertexHandle vh = _mesh->to_vertex_handle(heh);
            _mesh->set_color(vh, MyMesh::Color(0, 0, 180));
            _mesh->data(vh).thickness = 8;

            heh = _mesh->next_halfedge_handle(heh);
        }
    }

    // on affiche le nouveau maillage
    displayMesh(_mesh);
}

/**
 * Affiche les points adjacents du point sélectionné,
 * les faces et les points adjacents de l'arête sélectionné et
 * les faces, les arêtes et les points adjacents de la face sélectionné
 * @brief MainWindow::showSelectionsNeighborhood
 * @param _mesh
 */
void MainWindow::showSelectionsNeighborhood(MyMesh* _mesh)
{
    // on réinitialise les couleurs de tout le maillage
    resetAllColorsAndThickness(_mesh);

    //Affiche les arêtes adjacentes au point sélectionné
    if(vertexSelection >= 0 && vertexSelection < _mesh->n_vertices()){
        VertexHandle vh =_mesh->vertex_handle(vertexSelection); // VertexHandle d'indice i
        _mesh->set_color(vh, MyMesh::Color(255, 0, 0));
        _mesh->data(vh).thickness = 8;

        HalfedgeHandle init_heh = _mesh->halfedge_handle(vh); // HalfedgeHandle d'indice i

        HalfedgeHandle hh = init_heh;
        while(true){
            EdgeHandle eh = _mesh->edge_handle(hh); // EdgeHandle d'indice i
            _mesh->set_color(eh, MyMesh::Color(200, 0, 0));
            _mesh->data(eh).thickness = 4;
            hh = _mesh->opposite_halfedge_handle(hh);
            hh = _mesh->next_halfedge_handle(hh);
            if(hh == init_heh) break;
        }
    }

    //Affiche les faces adjacentes à l'arrête sélectionné
    if(edgeSelection >= 0 && edgeSelection < _mesh->n_edges()){
        EdgeHandle eh = _mesh->edge_handle(edgeSelection); // EdgeHandle d'indice i
        _mesh->set_color(eh, MyMesh::Color(0, 255, 0));
        _mesh->data(eh).thickness = 4;

        HalfedgeHandle heh0 = _mesh->halfedge_handle(eh, 0); // la première demi-arête
        HalfedgeHandle heh1 = _mesh->halfedge_handle(eh, 1); // la seconde demi-arête

        FaceHandle fh0 = _mesh->face_handle(heh0);
        if(fh0.is_valid()) _mesh->set_color(fh0, MyMesh::Color(0, 200, 0));

        FaceHandle fh1 = _mesh->face_handle(heh1);
        if(fh1.is_valid()) _mesh->set_color(fh1, MyMesh::Color(0, 200, 0));

        VertexHandle vh0 = _mesh->to_vertex_handle(heh0);
        _mesh->set_color(vh0, MyMesh::Color(0, 255, 0));
        _mesh->data(vh0).thickness = 8;

        VertexHandle vh1 = _mesh->to_vertex_handle(heh1);
        _mesh->set_color(vh1, MyMesh::Color(0, 255, 0));
        _mesh->data(vh1).thickness = 8;
    }

    //Affiche les faces adjacentes à la face sélectionné
    if(faceSelection >= 0 && faceSelection < _mesh->n_faces()) {
        FaceHandle fh = _mesh->face_handle(faceSelection); // FaceHandle d'indice i
        _mesh->set_color(fh, MyMesh::Color(0, 0, 100));

        for (MyMesh::FaceVertexIter curVertex = _mesh->fv_iter(fh); curVertex.is_valid(); curVertex ++){
            VertexHandle vh = *curVertex;
            _mesh->set_color(vh, MyMesh::Color(255, 0, 0));
            _mesh->data(vh).thickness = 8;
        }

        for (MyMesh::FaceEdgeIter curEdge = _mesh->fe_iter(fh); curEdge.is_valid(); curEdge ++)
        {
            EdgeHandle eh = *curEdge;
            _mesh->set_color(eh, MyMesh::Color(0, 255, 0));
            _mesh->data(eh).thickness = 4;
        }

        for (MyMesh::FaceFaceIter curFace = _mesh->ff_iter(fh); curFace.is_valid(); curFace ++){
            FaceHandle fh = *curFace;
            _mesh->set_color(fh, MyMesh::Color(0, 0, 255));
        }
    }

    // on affiche le nouveau maillage
    displayMesh(_mesh);
}

/**
 * Affiche la bordure de l'objet s'il y en a une
 * @brief MainWindow::showBorder
 * @param _mesh
 */
void MainWindow::showBorder(MyMesh* _mesh)
{
    // on réinitialise l'affichage
    resetAllColorsAndThickness(_mesh);

    // parcours des arêtes
    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
    {
        EdgeHandle eh = *curEdge;

        if(_mesh->is_boundary(eh))
            _mesh->set_color(eh, MyMesh::Color(255, 0, 255));
    }

    // on affiche le nouveau maillage
    displayMesh(_mesh);
}

/**
 * Calcul le poids de chaque arêtes
 * @brief MainWindow::calculEdgeWeight
 * @param _mesh
 * @param dist
 */
void MainWindow::calculEdgeWeight(MyMesh* _mesh, float dist[]){
    EdgeHandle eh;
    HalfedgeHandle heh;
    VertexHandle vhA;
    VertexHandle vhB;
    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++){
        eh = *curEdge;
        heh = _mesh->halfedge_handle(eh, 0);
        vhA = _mesh->from_vertex_handle(heh);
        vhB = _mesh->to_vertex_handle(heh);
        dist[eh.idx()] = ( _mesh->point(vhB) - _mesh->point(vhA) ).norm();
    }
}

/**
 * Calcul un chemin le plus court en nombre de sommets de source à dest
 * @brief MainWindow::dijkstra
 * @param _mesh
 * @param source
 * @param dest
 */
void MainWindow::dijkstraVertex(MyMesh* _mesh, int source, int dest){
    int dist[_mesh->n_vertices()];
    bool visited[_mesh->n_vertices()];
    int prev[_mesh->n_vertices()];

    //Initialisation des tableaux dist et visited
    for(size_t i = 0; i < _mesh->n_vertices(); ++i){
        dist[i] = INT_MAX;
        visited[i] = false;
    }

    VertexHandle vh_src = _mesh->vertex_handle(source);
    VertexHandle vh_dest = _mesh->vertex_handle(dest);

    VertexHandle vh_cur = vh_src;
    dist[vh_cur.idx()] = 0;
    int minVoisin;
    VertexHandle voisin;
    VertexHandle next;

    //Calcul le poids de chaque sommet jusqu'à ce que la destination soit atteinte
    while(vh_cur != vh_dest){
        minVoisin = INT_MAX;

        //Parcours chaque voisin du sommet courant
        for (MyMesh::VertexVertexIter curVertex = _mesh->vv_iter(vh_cur); curVertex.is_valid(); curVertex ++){
            voisin = *curVertex;

            if(voisin == vh_dest){
                dist[voisin.idx()] = dist[vh_cur.idx()] + 1;
                next = voisin;
                break;
            }

            if(visited[voisin.idx()]) continue;

            if(dist[voisin.idx()] > dist[vh_cur.idx()] + 1)
                dist[voisin.idx()] = dist[vh_cur.idx()] + 1;

            if(dist[voisin.idx()] < minVoisin){
                minVoisin = dist[voisin.idx()];
                next = _mesh->vertex_handle(voisin.idx());
            }
        }
        visited[vh_cur.idx()] = true;

        if(vh_cur == next) next = _mesh->vertex_handle(prev[vh_cur.idx()]);
        else prev[next.idx()] = vh_cur.idx();

        vh_cur = next;
    }

    vh_cur = vh_dest;
    HalfedgeHandle heh_color;
    VertexHandle vh_color;
    EdgeHandle eh_color;
    //Parcours le chemin dans le sens inverse pour afficher le chemin
    while(vh_cur != vh_src){
        minVoisin = INT_MAX;

        //Parcours chaque voisin du sommet courant
        for (MyMesh::VertexVertexIter curVertex = _mesh->vv_iter(vh_cur); curVertex.is_valid(); curVertex ++){
            voisin = *curVertex;

            if(voisin == vh_src){
                next = voisin;
                break;
            }

            if(dist[voisin.idx()] < minVoisin){
                minVoisin = dist[voisin.idx()];
                next = _mesh->vertex_handle(voisin.idx());
            }
        }

        //Coloriage des arêtes
        for (MyMesh::VertexOHalfedgeIter curHalfedge = _mesh->voh_iter(vh_cur); curHalfedge.is_valid(); curHalfedge ++){
            heh_color = *curHalfedge;
            vh_color = _mesh->to_vertex_handle(heh_color);
            if(next == vh_color){
                eh_color = _mesh->edge_handle(heh_color);
                _mesh->set_color(eh_color, MyMesh::Color(255, 0, 0));
                _mesh->data(eh_color).thickness = 4;
                break;
            }
        }

        vh_cur = next;
    }
    displayMesh(&mesh);

}

/**
 * Calcul un chemin le plus court en distance d'arêtes de source à dest
 * @brief MainWindow::dijkstra
 * @param _mesh
 * @param source
 * @param dest
 */
void MainWindow::dijkstraEdge(MyMesh* _mesh, int source, int dest){
    float dist[_mesh->n_vertices()];
    float edgeWeight[_mesh->n_edges()];
    bool visited[_mesh->n_vertices()];
    int prev[_mesh->n_vertices()];

    calculEdgeWeight(_mesh, edgeWeight);

    //Initialisation des tableaux dist et visited
    for(size_t i = 0; i < _mesh->n_vertices(); ++i){
        dist[i] = INT_MAX;
        visited[i] = false;
    }

    VertexHandle vh_src = _mesh->vertex_handle(source);
    VertexHandle vh_dest = _mesh->vertex_handle(dest);

    VertexHandle vh_cur = vh_src;
    dist[vh_cur.idx()] = 0;
    float minVoisin;
    float weight = 0;
    VertexHandle voisin;
    VertexHandle next;
    HalfedgeHandle heh_weight;
    VertexHandle vh_weight;
    EdgeHandle eh_weight;

    //Calcul le poids de chaque sommet jusqu'à ce que la destination soit atteinte
    while(vh_cur != vh_dest){
        minVoisin = INT_MAX;

        for (MyMesh::VertexVertexIter curVertex = _mesh->vv_iter(vh_cur); curVertex.is_valid(); curVertex ++){
            voisin = *curVertex;

            if(vh_cur == vh_dest){
                next = vh_cur;
                break;
            }

            if(visited[voisin.idx()]) continue;

            //Cherche le poids de l'arête
            for (MyMesh::VertexOHalfedgeIter curHalfedge = _mesh->voh_iter(vh_cur); curHalfedge.is_valid(); curHalfedge ++){
                heh_weight = *curHalfedge;
                vh_weight = _mesh->to_vertex_handle(heh_weight);
                if(voisin == vh_weight){
                    eh_weight = _mesh->edge_handle(heh_weight);
                    weight = edgeWeight[eh_weight.idx()];
                    break;
                }
            }

            if(dist[voisin.idx()] > dist[vh_cur.idx()] + weight)
                dist[voisin.idx()] = dist[vh_cur.idx()] + weight;

            if(dist[voisin.idx()] < minVoisin){
                minVoisin = dist[voisin.idx()];
                next = _mesh->vertex_handle(voisin.idx());
            }
        }
        visited[vh_cur.idx()] = true;

        if(vh_cur == next) next = _mesh->vertex_handle(prev[vh_cur.idx()]);
        else prev[next.idx()] = vh_cur.idx();

        vh_cur = next;
    }

    vh_cur = vh_dest;
    HalfedgeHandle heh_color;
    VertexHandle vh_color;
    EdgeHandle eh_color;
    //Parcours le chemin dans le sens inverse pour afficher le chemin
    while(vh_cur != vh_src){
        minVoisin = INT_MAX;

        for (MyMesh::VertexVertexIter curVertex = _mesh->vv_iter(vh_cur); curVertex.is_valid(); curVertex ++){
            voisin = *curVertex;

            if(voisin == vh_src){
                next = voisin;
                break;
            }

            if(dist[voisin.idx()] < minVoisin){
                minVoisin = dist[voisin.idx()];
                next = _mesh->vertex_handle(voisin.idx());
            }
        }

        //Coloriage des arêtes
        for (MyMesh::VertexOHalfedgeIter curHalfedge = _mesh->voh_iter(vh_cur); curHalfedge.is_valid(); curHalfedge ++){
            heh_color = *curHalfedge;
            vh_color = _mesh->to_vertex_handle(heh_color);
            if(next == vh_color){
                eh_color = _mesh->edge_handle(heh_color);
                _mesh->set_color(eh_color, MyMesh::Color(255, 0, 0));
                _mesh->data(eh_color).thickness = 4;
                break;
            }
        }

        vh_cur = next;
    }
    displayMesh(&mesh);

}

void MainWindow::showPath(MyMesh* _mesh, int v1, int v2, bool vertex)
{
    // on réinitialise l'affichage
    resetAllColorsAndThickness(_mesh);

    // point de départ et point d'arrivée en vert et en gros
    _mesh->set_color(_mesh->vertex_handle(v1), MyMesh::Color(0, 255, 0));
    _mesh->set_color(_mesh->vertex_handle(v2), MyMesh::Color(0, 255, 0));
    _mesh->data(_mesh->vertex_handle(v1)).thickness = 12;
    _mesh->data(_mesh->vertex_handle(v2)).thickness = 12;

    /* **** à compléter ! **** */
    if(vertex){
        dijkstraVertex(_mesh, v1, v2);
    } else {
        dijkstraEdge(_mesh, v1, v2);
    }
    // on affiche le nouveau maillage
    displayMesh(_mesh);
}

/* **** fin de la partie à compléter **** */


/* **** début de la partie boutons et IHM **** */

void MainWindow::on_pushButton_bordure_clicked()
{
    showBorder(&mesh);
}

void MainWindow::on_pushButton_voisinage_clicked()
{
    // changement de mode entre avec et sans voisinage
    if(modevoisinage)
    {
        ui->pushButton_voisinage->setText("Repasser en mode normal");
        modevoisinage = false;
    }
    else
    {
        ui->pushButton_voisinage->setText("Passer en mode voisinage");
        modevoisinage = true;
    }

    // on montre la nouvelle selection
    if(!modevoisinage)
        showSelections(&mesh);
    else
        showSelectionsNeighborhood(&mesh);
}


void MainWindow::on_pushButton_vertexMoins_clicked()
{
    // mise à jour de l'interface
    vertexSelection = vertexSelection - 1;
    ui->labelVertex->setText(QString::number(vertexSelection));

    // on montre la nouvelle selection
    if(!modevoisinage)
        showSelections(&mesh);
    else
        showSelectionsNeighborhood(&mesh);
}

void MainWindow::on_pushButton_vertexPlus_clicked()
{
    // mise à jour de l'interface
    vertexSelection = vertexSelection + 1;
    ui->labelVertex->setText(QString::number(vertexSelection));

    // on montre la nouvelle selection
    if(!modevoisinage)
        showSelections(&mesh);
    else
        showSelectionsNeighborhood(&mesh);
}

void MainWindow::on_pushButton_edgeMoins_clicked()
{
    // mise à jour de l'interface
    edgeSelection = edgeSelection - 1;
    ui->labelEdge->setText(QString::number(edgeSelection));

    // on montre la nouvelle selection
    if(!modevoisinage)
        showSelections(&mesh);
    else
        showSelectionsNeighborhood(&mesh);
}

void MainWindow::on_pushButton_edgePlus_clicked()
{
    // mise à jour de l'interface
    edgeSelection = edgeSelection + 1;
    ui->labelEdge->setText(QString::number(edgeSelection));

    // on montre la nouvelle selection
    if(!modevoisinage)
        showSelections(&mesh);
    else
        showSelectionsNeighborhood(&mesh);
}

void MainWindow::on_pushButton_faceMoins_clicked()
{
    // mise à jour de l'interface
    faceSelection = faceSelection - 1;
    ui->labelFace->setText(QString::number(faceSelection));

    // on montre la nouvelle selection
    if(!modevoisinage)
        showSelections(&mesh);
    else
        showSelectionsNeighborhood(&mesh);
}

void MainWindow::on_pushButton_facePlus_clicked()
{
    // mise à jour de l'interface
    faceSelection = faceSelection + 1;
    ui->labelFace->setText(QString::number(faceSelection));

    // on montre la nouvelle selection
    if(!modevoisinage)
        showSelections(&mesh);
    else
        showSelectionsNeighborhood(&mesh);
}

void MainWindow::on_pushButton_afficherChemin_clicked()
{
    // on récupère les sommets de départ et d'arrivée
    int indexV1 = ui->spinBox_v1_chemin->value();
    int indexV2 = ui->spinBox_v2_chemin->value();

    showPath(&mesh, indexV1, indexV2, true);
}

void MainWindow::on_pushButton_aretes_clicked()
{
    // on récupère les sommets de départ et d'arrivée
    int indexV1 = ui->spinBox_v1_chemin->value();
    int indexV2 = ui->spinBox_v2_chemin->value();

    showPath(&mesh, indexV1, indexV2, false);
}


void MainWindow::on_pushButton_chargement_clicked()
{
    // fenêtre de sélection des fichiers
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.obj)"));

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&mesh);

    // on affiche le maillage
    displayMesh(&mesh);
}

/* **** fin de la partie boutons et IHM **** */



/* **** fonctions supplémentaires **** */

// permet d'initialiser les couleurs et les épaisseurs des élements du maillage
void MainWindow::resetAllColorsAndThickness(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        _mesh->data(*curVert).thickness = 1;
        _mesh->set_color(*curVert, MyMesh::Color(0, 0, 0));
    }

    for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
    {
        _mesh->set_color(*curFace, MyMesh::Color(150, 150, 150));
    }

    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
    {
        _mesh->data(*curEdge).thickness = 1;
        _mesh->set_color(*curEdge, MyMesh::Color(0, 0, 0));
    }
}

// charge un objet MyMesh dans l'environnement OpenGL
void MainWindow::displayMesh(MyMesh* _mesh)
{
    GLuint* triIndiceArray = new GLuint[_mesh->n_faces() * 3];
    GLfloat* triCols = new GLfloat[_mesh->n_faces() * 3 * 3];
    GLfloat* triVerts = new GLfloat[_mesh->n_faces() * 3 * 3];

    MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
    MyMesh::ConstFaceVertexIter fvIt;
    int i = 0;
    for (; fIt!=fEnd; ++fIt)
    {
        fvIt = _mesh->cfv_iter(*fIt);
        triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
        triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
        triIndiceArray[i] = i;

        i++; ++fvIt;
        triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
        triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
        triIndiceArray[i] = i;

        i++; ++fvIt;
        triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
        triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
        triIndiceArray[i] = i;

        i++;
    }

    ui->displayWidget->loadMesh(triVerts, triCols, _mesh->n_faces() * 3 * 3, triIndiceArray, _mesh->n_faces() * 3);

    delete[] triIndiceArray;
    delete[] triCols;
    delete[] triVerts;

    GLuint* linesIndiceArray = new GLuint[_mesh->n_edges() * 2];
    GLfloat* linesCols = new GLfloat[_mesh->n_edges() * 2 * 3];
    GLfloat* linesVerts = new GLfloat[_mesh->n_edges() * 2 * 3];

    i = 0;
    QHash<float, QList<int> > edgesIDbyThickness;
    for (MyMesh::EdgeIter eit = _mesh->edges_begin(); eit != _mesh->edges_end(); ++eit)
    {
        float t = _mesh->data(*eit).thickness;
        if(t > 0)
        {
            if(!edgesIDbyThickness.contains(t))
                edgesIDbyThickness[t] = QList<int>();
            edgesIDbyThickness[t].append((*eit).idx());
        }
    }
    QHashIterator<float, QList<int> > it(edgesIDbyThickness);
    QList<QPair<float, int> > edgeSizes;
    while (it.hasNext())
    {
        it.next();

        for(int e = 0; e < it.value().size(); e++)
        {
            int eidx = it.value().at(e);

            MyMesh::VertexHandle vh1 = _mesh->to_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh1)[0];
            linesVerts[3*i+1] = _mesh->point(vh1)[1];
            linesVerts[3*i+2] = _mesh->point(vh1)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;

            MyMesh::VertexHandle vh2 = _mesh->from_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh2)[0];
            linesVerts[3*i+1] = _mesh->point(vh2)[1];
            linesVerts[3*i+2] = _mesh->point(vh2)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;
        }
        edgeSizes.append(qMakePair(it.key(), it.value().size()));
    }

    ui->displayWidget->loadLines(linesVerts, linesCols, i * 3, linesIndiceArray, i, edgeSizes);

    delete[] linesIndiceArray;
    delete[] linesCols;
    delete[] linesVerts;

    GLuint* pointsIndiceArray = new GLuint[_mesh->n_vertices()];
    GLfloat* pointsCols = new GLfloat[_mesh->n_vertices() * 3];
    GLfloat* pointsVerts = new GLfloat[_mesh->n_vertices() * 3];

    i = 0;
    QHash<float, QList<int> > vertsIDbyThickness;
    for (MyMesh::VertexIter vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
    {
        float t = _mesh->data(*vit).thickness;
        if(t > 0)
        {
            if(!vertsIDbyThickness.contains(t))
                vertsIDbyThickness[t] = QList<int>();
            vertsIDbyThickness[t].append((*vit).idx());
        }
    }
    QHashIterator<float, QList<int> > vitt(vertsIDbyThickness);
    QList<QPair<float, int> > vertsSizes;

    while (vitt.hasNext())
    {
        vitt.next();

        for(int v = 0; v < vitt.value().size(); v++)
        {
            int vidx = vitt.value().at(v);

            pointsVerts[3*i+0] = _mesh->point(_mesh->vertex_handle(vidx))[0];
            pointsVerts[3*i+1] = _mesh->point(_mesh->vertex_handle(vidx))[1];
            pointsVerts[3*i+2] = _mesh->point(_mesh->vertex_handle(vidx))[2];
            pointsCols[3*i+0] = _mesh->color(_mesh->vertex_handle(vidx))[0];
            pointsCols[3*i+1] = _mesh->color(_mesh->vertex_handle(vidx))[1];
            pointsCols[3*i+2] = _mesh->color(_mesh->vertex_handle(vidx))[2];
            pointsIndiceArray[i] = i;
            i++;
        }
        vertsSizes.append(qMakePair(vitt.key(), vitt.value().size()));
    }

    ui->displayWidget->loadPoints(pointsVerts, pointsCols, i * 3, pointsIndiceArray, i, vertsSizes);

    delete[] pointsIndiceArray;
    delete[] pointsCols;
    delete[] pointsVerts;
}


MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    vertexSelection = -1;
    edgeSelection = -1;
    faceSelection = -1;

    modevoisinage = false;

    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

