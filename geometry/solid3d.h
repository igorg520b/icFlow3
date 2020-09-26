#ifndef SOLID3D_H
#define SOLID3D_H

#include <QObject>
#include <QTreeWidgetItem>
#include <QDebug>
#include <vector>
#include <Eigen/Core>
#include <gmsh.h>

#include <vtkNew.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellType.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkNamedColors.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkPolyDataMapper.h>
#include <vtkDataSetMapper.h>
#include <vtkLookupTable.h>

namespace icy { class Solid3D; }

class icy::Solid3D : public QObject
{
    Q_OBJECT
    Q_PROPERTY(QString Name MEMBER Name NOTIFY propertyChanged)

public:
    Solid3D();
    Solid3D(QString fileName);  // calls ImportSTL
    void ImportSTL(QString fileName);
    void CreateUGrid();
    QTreeWidgetItem treeWidget;     // to represent in GUI

    std::vector<Eigen::Vector3d> nodes; // coordinates of 3D nodes
    std::vector<Eigen::Vector3i> tris;  // boundary of the object
    QString Name;

    vtkNew<vtkPoints> points;

    vtkNew<vtkUnstructuredGrid> ugrid;
    vtkNew<vtkDataSetMapper> dataSetMapper;
    vtkNew<vtkActor> actor_mesh;

    vtkNew<vtkDoubleArray> forces_nodes;
    vtkNew<vtkDataSetMapper> mapper;

    vtkNew<vtkNamedColors> colors;
    vtkNew<vtkLookupTable> hueLut;

signals:
    void propertyChanged();
};

#endif // SOLID3D_H
