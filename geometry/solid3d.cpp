#include <QFileInfo>
#include "solid3d.h"

namespace model = gmsh::model;
//namespace factory = gmsh::model::occ;

icy::Solid3D::Solid3D() { }


icy::Solid3D::Solid3D(QString fileName)
{
    ImportSTL(fileName);
}


void icy::Solid3D::ImportSTL(QString fileName)
{
    qDebug() << "importing 3D shell from STL file";
    if(fileName.isEmpty()) return;

    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::open(fileName.toStdString());

    // get tris
    std::vector<std::size_t> trisTags, nodeTagsInFaces;
    model::mesh::getElementsByType(2, trisTags, nodeTagsInFaces);
    std::cout << "triangles " << trisTags.size() << "nodetags " << nodeTagsInFaces.size() << std::endl;

    // nodes of tris
    std::vector<std::size_t> nodeTags4;
    std::vector<double> nodeCoords4, parametricCoords4;
    model::mesh::getNodesByElementType(2, nodeTags4, nodeCoords4, parametricCoords4, -1, false);
    std::cout << "nds of tris: " << nodeTags4.size() << "; coords " << nodeCoords4.size() << std::endl;

    gmsh::clear();

    //=============================== read vectors => ::nodes, ::tris

    std::map<std::size_t, int> mtags;

    std::size_t count=0;
    for(std::size_t i=0;i<nodeTags4.size();i++) {
        std::size_t gmshTag = nodeTags4[i];
        auto iter = mtags.find(gmshTag);
        if(iter != mtags.end()) continue;

        double x = nodeCoords4[i*3+0];
        double y = nodeCoords4[i*3+1];
        double z = nodeCoords4[i*3+2];

        nodes.push_back(Eigen::Vector3d(x,y,z));
        std::pair<std::map<std::size_t,int>::iterator,bool> result =
                mtags.insert({gmshTag,count++});
        if(result.second==false) throw std::runtime_error("duplicate node");
    }

//    trisTags, nodeTagsInFaces;
    tris.reserve(trisTags.size());
    for(std::size_t i=0;i<trisTags.size();i++) {
        int id1=mtags[nodeTagsInFaces[3*i+0]];
        int id2=mtags[nodeTagsInFaces[3*i+1]];
        int id3=mtags[nodeTagsInFaces[3*i+2]];
        tris.push_back(Eigen::Vector3i(id1,id2,id3));
    }

    CreateUGrid();

    QFileInfo fileInfo(fileName);
    Name = fileInfo.baseName();
    treeWidget.setText(0, Name);
    treeWidget.setData(0, Qt::UserRole, QVariant::fromValue(this));
    treeWidget.setData(1, Qt::UserRole, QVariant::fromValue(2));
}


void icy::Solid3D::CreateUGrid()
{
    points->Reset();
    points->Allocate(nodes.size());

    vtkIdType count=0;
    for(auto const &nd : nodes) points->InsertPoint(count++, nd.x(), nd.y(), nd.z());
    // triangular mesh
    ugrid->Reset();
    ugrid->Allocate(tris.size());
    ugrid->SetPoints(points);

    vtkIdType pts2[3];
    for(auto const &tr : tris)
    {
        pts2[0] = tr.coeff(0);
        pts2[1] = tr.coeff(1);
        pts2[2] = tr.coeff(2);
        ugrid->InsertNextCell(VTK_TRIANGLE,3,pts2);
    }

    dataSetMapper->SetInputData(ugrid);
    actor_mesh->SetMapper(dataSetMapper);
    actor_mesh->GetProperty()->SetColor(colors->GetColor3d("Pink").GetData());
    actor_mesh->GetProperty()->EdgeVisibilityOn();
}
