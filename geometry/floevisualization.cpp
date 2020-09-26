#include "floevisualization.h"
#include "node.h"

icy::FloeVisualization::FloeVisualization()
{
    selectedPointId = -1;

    ugrid->SetPoints(points);
    ugrid_vertices->SetPoints(points);
    ugrid_boundary->SetPoints(points);
    ugrid_selection->SetPoints(points);

    actor_mesh->PickableOn();
    actor_boundary->PickableOff();
    actor_arrows->PickableOff();
    actor_labels->PickableOff();

    actor_mesh->GetProperty()->VertexVisibilityOff();
    actor_boundary->GetProperty()->VertexVisibilityOff();
    actor_arrows->GetProperty()->VertexVisibilityOff();

    dataSetMapper->SetInputData(ugrid);
    dataSetMapper->UseLookupTableScalarRangeOn();
    actor_mesh->SetMapper(dataSetMapper);
//    actor_mesh->GetProperty()->SetColor(colors->GetColor3d("Seashell").GetData());
    actor_mesh->GetProperty()->SetColor(218/255.0,228/255.0,242/255.0);
    actor_mesh->GetProperty()->SetEdgeColor(161.0/255.0, 176.0/255.0, 215.0/255.0);
    actor_mesh->GetProperty()->LightingOff();
    actor_mesh->GetProperty()->ShadingOff();
    actor_mesh->GetProperty()->SetInterpolationToFlat();

    dataSetMapper_boundary->SetInputData(ugrid_boundary);
    actor_boundary->SetMapper(dataSetMapper_boundary);
    actor_boundary->GetProperty()->SetColor(colors->GetColor3d("Black").GetData());
    actor_boundary->GetProperty()->EdgeVisibilityOn();
    actor_boundary->GetProperty()->SetLineWidth(2);

    visualized_values->SetName("visualized_values");

    edgeNumbers->SetNumberOfComponents(1);
    edgeNumbers->SetName("edgeNumbers");

    // arrows
    arrowCoords->SetNumberOfComponents(3);
    arrowCoords->SetName("arrowCoords");
    ugrid_vertices->GetPointData()->AddArray(arrowCoords);

    arrowSource->SetShaftRadius(0.02);
    arrowSource->SetTipRadius(0.03);
    arrowSource->SetTipLength(0.06);

    glyph3D->SetSourceConnection(arrowSource->GetOutputPort());
    glyph3D->SetVectorModeToUseVector();
    glyph3D->SetInputData(ugrid_vertices);
    glyph3D->OrientOn();
    glyph3D->ScalingOn();
    glyph3D->SetScaleModeToScaleByVector();
    glyph3D->SetScaleFactor(0.1);
    mapper_arrows->SetInputConnection(glyph3D->GetOutputPort());
    actor_arrows->SetMapper(mapper_arrows);
    actor_arrows->GetProperty()->SetColor(colors->GetColor3d("Black").GetData());

    actor_labels->VisibilityOff();

    InitializeLUT();

    // water level
    points_water->SetNumberOfPoints(gridSize*gridSize);
    grid_water->SetDimensions(gridSize, gridSize, 1);
    grid_water->SetPoints(points_water);
    mapper_water->SetInputData(grid_water);
    actor_water->SetMapper(mapper_water);
    actor_water->GetProperty()->SetRepresentationToWireframe();
    actor_water->GetProperty()->SetColor(colors->GetColor3d("Black").GetData());

    // mesh 3d
    ugrid_mesh3d->SetPoints(points_mesh3d);
    mapper_mesh3d->SetInputData(ugrid_mesh3d);
    actor_mesh_mesh3d->SetMapper(mapper_mesh3d);
    mapper_mesh3d->UseLookupTableScalarRangeOn();
}

void icy::FloeVisualization::UnsafeUpdateTopology(std::vector<icy::Node*> *nodes,
                                                  std::vector<icy::Element*> *elems,
                                                  std::vector<icy::Edge*> *edges)
{
    if(selectedPointId >= 0) UnsafeUpdateSelection(nodes, -1);

    cellArray->Reset();
    cellArray_boundary->Reset();
    cellArray_vertices->Reset();
    UnsafeUpdateDisplacements(nodes, elems);

    // create ugrid
    vtkIdType pts2[3];
    for(icy::Element *tr : *elems)
    {
        for(int j=0;j<3;j++) pts2[j] = tr->nds[j]->locId;
        cellArray->InsertNextCell(3, pts2);
    }
    ugrid->SetCells(VTK_TRIANGLE, cellArray);

    // ugrid_boundary
    for(icy::Edge *edge : *edges)
    {
        if(!edge->isBoundary) continue;
        pts2[0] = edge->nds[0]->locId;
        pts2[1] = edge->nds[1]->locId;
        cellArray_boundary->InsertNextCell(2, pts2);
    }
    ugrid_boundary->SetCells(VTK_LINE, cellArray_boundary);

    UnsafeUpdateArrows(nodes);
}

void icy::FloeVisualization::UnsafeUpdateDisplacements(std::vector<icy::Node*> *nodes,
                                                       std::vector<icy::Element*> *elems)
{
//    qDebug() << "UnsafeUpdateDisplacements()";

    points->SetNumberOfPoints(nodes->size());
    vtkIdType count=0;
    if(use_tentative_coordinates)
        for(icy::Node* nd : *nodes) points->SetPoint(count++, nd->xt.data());
    else
        for(icy::Node* nd : *nodes) points->SetPoint(count++, nd->xn.data());
    points->Modified();
    UnsafeUpdateValues(nodes, elems);


}

void icy::FloeVisualization::UnsafeUpdateValues(std::vector<icy::Node*> *nodes,
                                                std::vector<icy::Element*> *elems,
                                                int option)
{
    if(option >= 0) VisualizingVariable = (VisOpt)option;

    if(VisualizingVariable == VisOpt::none || nodes->size() == 0) {
        dataSetMapper->ScalarVisibilityOff();
        ugrid->GetPointData()->RemoveArray("visualized_values");
        ugrid->GetCellData()->RemoveArray("visualized_values");
        return;
    }

    visualized_values->SetNumberOfValues(nodes->size());

    switch(VisualizingVariable)
    {
    case VisOpt::vert_force:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->vertical_force);
        dataSetMapper->SetLookupTable(hueLut);
        break;

    case VisOpt::boundary:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->isBoundary ? 1 : 0);
        dataSetMapper->SetLookupTable(hueLut);
        break;

    case VisOpt::deflection:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->xn.z());
        dataSetMapper->SetLookupTable(hueLut);
        break;

    case VisOpt::max_normal_traction:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->max_normal_traction);
        dataSetMapper->SetLookupTable(hueLut);
        break;

    case VisOpt::fracture_support:
        for(icy::Node* nd : *nodes) {
            double value;
            if(nd->crack_tip) value = 1.3;
            else if(nd->core_node) value = 1;
            else if(nd->support_node) value = 0.5;
            else value = 0;
            visualized_values->SetValue(nd->locId, value);
        }
        dataSetMapper->SetLookupTable(defaultLut);
        defaultLut->SetTableRange(-1, 1.3);
        break;

    case VisOpt::AbsMx:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, std::abs(nd->str_b[0]));
        dataSetMapper->SetLookupTable(defaultLut);
        break;

    case VisOpt::Mx:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_b[0]);
        dataSetMapper->SetLookupTable(defaultLut);
        break;

    case VisOpt::My:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_b[1]);
        dataSetMapper->SetLookupTable(defaultLut);
        break;

    case VisOpt::Mxy:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_b[2]);
        dataSetMapper->SetLookupTable(defaultLut);
        break;

    case VisOpt::Mx_e:
        visualized_values->SetNumberOfValues(elems->size());
        for(std::size_t i=0;i<elems->size();i++) visualized_values->SetValue(i, (*elems)[i]->str_b[0]);
        dataSetMapper->SetLookupTable(defaultLut);
        break;

    case VisOpt::My_e:
        visualized_values->SetNumberOfValues(elems->size());
        for(std::size_t i=0;i<elems->size();i++) visualized_values->SetValue(i, (*elems)[i]->str_b[1]);
        dataSetMapper->SetLookupTable(defaultLut);
        break;

    case VisOpt::Mxy_e:
        visualized_values->SetNumberOfValues(elems->size());
        for(std::size_t i=0;i<elems->size();i++) visualized_values->SetValue(i, (*elems)[i]->str_b[2]);
        dataSetMapper->SetLookupTable(defaultLut);
        break;

    case VisOpt::Tx:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_m[0]);
        dataSetMapper->SetLookupTable(defaultLut);
        break;

    case VisOpt::Ty:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_m[1]);
        dataSetMapper->SetLookupTable(defaultLut);
        break;

    case VisOpt::Txy:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_m[2]);
        dataSetMapper->SetLookupTable(defaultLut);
        break;

    case VisOpt::Qx:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_s[0]);
        dataSetMapper->SetLookupTable(defaultLut);
        break;

    case VisOpt::Qy:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_s[1]);
        dataSetMapper->SetLookupTable(defaultLut);
        break;

    case VisOpt::stx:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_b_top[0]);
        dataSetMapper->SetLookupTable(defaultLut);
        break;

    case VisOpt::sty:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_b_top[1]);
        dataSetMapper->SetLookupTable(defaultLut);
        break;

    case VisOpt::stxy:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_b_top[2]);
        dataSetMapper->SetLookupTable(defaultLut);
        break;

    case VisOpt::st1:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_b_top_principal[0]);
        dataSetMapper->SetLookupTable(defaultLut);
        break;

    case VisOpt::st2:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_b_top_principal[1]);
        dataSetMapper->SetLookupTable(defaultLut);
        break;

    case VisOpt::sbx:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_b_bottom[0]);
        dataSetMapper->SetLookupTable(defaultLut);
        break;

    case VisOpt::sby:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_b_bottom[1]);
        dataSetMapper->SetLookupTable(defaultLut);
        break;

    case VisOpt::sbxy:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_b_bottom[2]);
        dataSetMapper->SetLookupTable(defaultLut);
        break;

    case VisOpt::sb1:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_b_bottom_principal[0]);
        dataSetMapper->SetLookupTable(defaultLut);
        break;

    case VisOpt::sb2:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_b_bottom_principal[1]);
        dataSetMapper->SetLookupTable(defaultLut);
        break;

    default:
        break;
    }

    visualized_values->Modified();

    if(VisualizingVariable == VisOpt::Mxy_e ||
            VisualizingVariable == VisOpt::Mx_e ||
            VisualizingVariable == VisOpt::My_e)
    {
        // elements
        ugrid->GetPointData()->RemoveArray("visualized_values");
        ugrid->GetCellData()->AddArray(visualized_values);
        ugrid->GetCellData()->SetActiveScalars("visualized_values");
        dataSetMapper->SetScalarModeToUseCellData();
    }
    else
    {
        // nodes
        ugrid->GetCellData()->RemoveArray("visualized_values");
        ugrid->GetPointData()->AddArray(visualized_values);
        ugrid->GetPointData()->SetActiveScalars("visualized_values");
        dataSetMapper->SetScalarModeToUsePointData();
    }

    dataSetMapper->ScalarVisibilityOn();
    dataSetMapper->SetColorModeToMapScalars();

    if(update_minmax && VisualizingVariable != VisOpt::fracture_support) {
        double minmax[2];
        visualized_values->GetValueRange(minmax);
        hueLut->SetTableRange(minmax[0], minmax[1]);
        defaultLut->SetTableRange(minmax[0], minmax[1]);
    }


    UnsafeUpdateArrows(nodes);
}

void icy::FloeVisualization::InitializeLUT()
{
    hueLut->SetNumberOfTableValues(257);
    for ( int i=0; i<257; i++) {
            hueLut->SetTableValue(i, (double)lutArray[i][0], (double)lutArray[i][1], (double)lutArray[i][2], 1.0);
    }
}

void icy::FloeVisualization::UnsafeUpdateSelection(std::vector<icy::Node*> *nodes,
                                                   vtkIdType selectedPoint)
{
    //    ugrid_selection
    edgeNumbers->Reset();
    if(selectedPointId >= 0 && selectedPoint < 0)
    {
        // remove visualization
        actor_labels->VisibilityOff();
        ugrid_selection->Reset();
    }
    else if(selectedPoint > 0)
    {
        ugrid_selection->Reset();
        ugrid_selection->SetPoints(points);
        icy::Node *nd = (*nodes)[selectedPoint];

        vtkIdType pts2[3];
        std::size_t nFan = nd->fan.size();
        for(std::size_t i=0;i<nFan;i++)
        {
            icy::Node::FanPrecomp &f = nd->fan[i];
            for(int j=0;j<3;j++) pts2[j] = f.face->nds[j]->locId;
            ugrid_selection->InsertNextCell(VTK_TRIANGLE,3,pts2);
            edgeNumbers->InsertNextValue(i);
        }
        ugrid_selection->GetCellData()->AddArray(edgeNumbers);
        ugrid_selection->GetCellData()->SetActiveScalars("edgeNumbers");
        geometryFilter->SetInputData(ugrid_selection);
        geometryFilter->Update();
        idfilter->SetInputConnection(geometryFilter->GetOutputPort());
        idfilter->PointIdsOff();
        idfilter->CellIdsOff();
        idfilter->FieldDataOn();
        idfilter->Update();
        cellCenters->SetInputConnection(idfilter->GetOutputPort());
        cellCenters->Update();
        labledDataMapper->SetInputConnection(cellCenters->GetOutputPort());
        labledDataMapper->SetLabelModeToLabelScalars();
        labledDataMapper->Update();
        actor_labels->SetMapper(labledDataMapper);
        actor_labels->GetProperty()->SetColor(colors->GetColor3d("Black").GetData());
        actor_labels->VisibilityOn();
    }
    selectedPointId = selectedPoint;
}



void icy::FloeVisualization::UnsafeUpdateArrows(std::vector<icy::Node*> *nodes)
{
    actor_arrows->SetVisibility(update_arrows);
    if(!update_arrows)
    {
        ugrid_vertices->GetPointData()->RemoveArray("arrowCoords");
        return;
    }
    std::size_t nNodes = nodes->size();
    arrowCoords->SetNumberOfTuples(nNodes);

    // vertices (for arrows)
    cellArray_vertices->Reset();
    vtkIdType pt;
    for(std::size_t i=0;i<nNodes;i++)
    {
        icy::Node *nd = (*nodes)[i];
        arrowCoords->SetTuple(i, nd->dir.data());
        pt=nd->locId;
        cellArray_vertices->InsertNextCell(1, &pt);
    }
    ugrid_vertices->SetCells(VTK_VERTEX, cellArray_vertices);

    ugrid_vertices->GetPointData()->AddArray(arrowCoords);
    ugrid_vertices->GetPointData()->SetActiveVectors("arrowCoords");
}


void icy::FloeVisualization::UnsafeUpdateWaterLine(int mode, double totalTime)
{
    if(mode == 1)
    {
        // update water level
        for(unsigned i=0;i<gridSize;i++)
            for(unsigned j=0;j<gridSize;j++) {
                double x = 10.0 * ((double)i/(double)gridSize-0.5);
                double y = 5.0 * ((double)j/(double)gridSize-0.5);
                double riverRapids = 0.5*icy::Node::RiverRapids(x, totalTime*0.3-4);
                points_water->SetPoint(j+i*gridSize, x,y,riverRapids);
            }
        points_water->Modified();
    }
    else if(mode == 2)
    {
        for(unsigned i=0;i<gridSize;i++)
            for(unsigned j=0;j<gridSize;j++) {
                double x = 10.0 * ((double)i/(double)gridSize-0.5);
                double y = 5.0 * ((double)j/(double)gridSize-0.5);
                points_water->SetPoint(j+i*gridSize, x,y,0);
            }
        points_water->Modified();
    }
    else if(mode == 3)
    {
        for(unsigned i=0;i<gridSize;i++)
            for(unsigned j=0;j<gridSize;j++) {
                double x = 10.0 * ((double)i/(double)gridSize-0.5);
                double y = 5.0 * ((double)j/(double)gridSize-0.5);
                double rest_position = sin(x*1.5);
                rest_position *= 0.07*sin(totalTime * 2 * M_PI / 10);
                if(totalTime < 10) rest_position *= totalTime/10;
                points_water->SetPoint(j+i*gridSize, x,y,rest_position);
            }
        points_water->Modified();

    }
    else if(mode == 4)
    {
        for(unsigned i=0;i<gridSize;i++)
            for(unsigned j=0;j<gridSize;j++) {
                double x = 10.0 * ((double)i/(double)gridSize-0.5);
                double y = 5.0 * ((double)j/(double)gridSize-0.5);
                double rest_position = sin(x*1)+sin(y*1+1);
                rest_position *= 0.05*sin(totalTime * 2 * M_PI / 10);
                if(totalTime < 10) rest_position *= totalTime/10;
                points_water->SetPoint(j+i*gridSize, x,y,rest_position);
            }
        points_water->Modified();
    }



}
