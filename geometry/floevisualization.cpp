#include "floevisualization.h"
#include "node.h"
#include "model.h"

icy::FloeVisualization::FloeVisualization()
{
    selectedPointId = -1;

    ugrid->SetPoints(points);
    ugrid_boundary->SetPoints(points);
    ugrid_selection->SetPoints(points);

    actor_mesh->PickableOn();
    actor_boundary->PickableOff();
//    actor_arrows->PickableOff();
    actor_labels->PickableOff();

    actor_mesh->GetProperty()->VertexVisibilityOff();
    actor_boundary->GetProperty()->VertexVisibilityOff();
//    actor_arrows->GetProperty()->VertexVisibilityOff();

    dataSetMapper->SetInputData(ugrid);
    dataSetMapper->UseLookupTableScalarRangeOn();
    dataSetMapper->SetLookupTable(hueLut);
    actor_mesh->SetMapper(dataSetMapper);
//    actor_mesh->GetProperty()->SetColor(218/255.0,228/255.0,242/255.0);
    actor_mesh->GetProperty()->SetColor(0.54938, 0.772213, 0.848103);


//    actor_mesh->GetProperty()->SetEdgeColor(161.0/255.0, 176.0/255.0, 215.0/255.0);
    actor_mesh->GetProperty()->SetEdgeColor(91.0/255.0, 116.0/255.0, 145.0/255.0);
    actor_mesh->GetProperty()->SetLineWidth(1.5);
    actor_mesh->GetProperty()->LightingOff();
    actor_mesh->GetProperty()->ShadingOff();
    actor_mesh->GetProperty()->SetInterpolationToFlat();

    dataSetMapper_boundary->SetInputData(ugrid_boundary);
    actor_boundary->SetMapper(dataSetMapper_boundary);
    actor_boundary->GetProperty()->SetColor(colors->GetColor3d("Black").GetData());
    actor_boundary->GetProperty()->EdgeVisibilityOn();
    actor_boundary->GetProperty()->SetLineWidth(4);
//    actor_boundary->GetProperty()->SetLineWidth(1.9);

    visualized_values->SetName("visualized_values");

    edgeNumbers->SetNumberOfComponents(1);
    edgeNumbers->SetName("edgeNumbers");

    // arrows
/*
    ugrid_vertices->SetPoints(points);

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
    glyph3D->SetScaleFactor(0.025);
    mapper_arrows->SetInputConnection(glyph3D->GetOutputPort());
    actor_arrows->SetMapper(mapper_arrows);
    actor_arrows->GetProperty()->SetColor(colors->GetColor3d("Black").GetData());
    actor_arrows->PickableOff();
*/
    actor_labels->VisibilityOff();

    // water level
    points_water->SetNumberOfPoints(gridSizeX*gridSizeY);
    grid_water->SetDimensions(gridSizeX, gridSizeY, 1);
    grid_water->SetPoints(points_water);
    mapper_water->SetInputData(grid_water);
    actor_water->SetMapper(mapper_water);
    actor_water->GetProperty()->SetRepresentationToWireframe();
    actor_water->GetProperty()->SetColor(colors->GetColor3d("Black").GetData());
    actor_water->PickableOff();

    // indenter
    sphereSource->SetCenter(0.0, 0.0, 0.0);
    sphereSource->SetRadius(1.0);
    // Make the surface smooth.
    sphereSource->SetPhiResolution(15);
    sphereSource->SetThetaResolution(15);
    sphereMapper->SetInputConnection(sphereSource->GetOutputPort());
    actor_sphere->SetMapper(sphereMapper);
    actor_sphere->GetProperty()->SetColor(colors->GetColor3d("Black").GetData());
    actor_sphere->GetProperty()->SetRepresentationToWireframe();
    actor_sphere->GetProperty()->SetLineWidth(1.5);
    actor_sphere->GetProperty()->LightingOff();
    actor_sphere->GetProperty()->ShadingOff();
    actor_sphere->GetProperty()->SetEdgeColor(colors->GetColor3d("Black").GetData());

    actor_sphere->VisibilityOff();
}

void icy::FloeVisualization::UnsafeUpdateTopology(std::vector<Node*> *nodes, std::vector<Element*> *elems,
                                                  std::vector<Edge> &boundaryEdges,
                                                  double temporalThreshold)
{
    if(selectedPointId >= 0) UnsafeUpdateSelection(nodes, -1);

    cellArray->Reset();
    cellArray_boundary->Reset();
//    cellArray_vertices->Reset();
    UnsafeUpdateDisplacements(nodes, elems, temporalThreshold);

    // create ugrid
    vtkIdType pts2[3];
    for(icy::Element *tr : *elems)
    {
        for(int j=0;j<3;j++) pts2[j] = tr->nds[j]->locId;
        cellArray->InsertNextCell(3, pts2);
    }
    ugrid->SetCells(VTK_TRIANGLE, cellArray);

    // ugrid_boundary
    for(icy::Edge &edge : boundaryEdges)
    {
        pts2[0] = edge.nds[0]->locId;
        pts2[1] = edge.nds[1]->locId;
        cellArray_boundary->InsertNextCell(2, pts2);
    }
    ugrid_boundary->SetCells(VTK_LINE, cellArray_boundary);

    UnsafeUpdateArrows(nodes);
}

void icy::FloeVisualization::UnsafeUpdateDisplacements(std::vector<Node*> *nodes,
                                                       std::vector<Element*> *elems, double temporalThreshold)
{
    points->SetNumberOfPoints(nodes->size());
    vtkIdType count=0;
    if(use_tentative_coordinates)
        for(icy::Node* nd : *nodes) points->SetPoint(count++, nd->xt.data());
    else
        for(icy::Node* nd : *nodes) points->SetPoint(count++, nd->xn.data());
    points->Modified();
    UnsafeUpdateValues(nodes, elems, temporalThreshold);
}

void icy::FloeVisualization::UnsafeUpdateValues(std::vector<Node*> *nodes,
                                                std::vector<Element*> *elems, double temporalThreshold,
                                                int option)
{
    if(option >= 0)
    {
        VisualizingVariable = (VisOpt)option;
        // set lookup table
        if(VisualizingVariable == VisOpt::vert_force ||
                VisualizingVariable == VisOpt::boundary ||
                VisualizingVariable == VisOpt::deflection ||
                VisualizingVariable == VisOpt::disp_x ||
                VisualizingVariable == VisOpt::max_normal_traction ||
                VisualizingVariable == VisOpt::time_loaded) InitializeLUT(1);
        else if(VisualizingVariable == VisOpt::fracture_support) InitializeLUT(4);
        else if(VisualizingVariable == VisOpt::region) InitializeLUT(5);
        else InitializeLUT(3);
    }

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
        break;

    case VisOpt::boundary:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->isBoundary ? 1 : 0);
        break;

    case VisOpt::deflection:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->xn.z());
        break;

    case VisOpt::disp_x:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->un.x());
        break;

    case VisOpt::max_normal_traction:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->max_normal_traction);
        break;

    case VisOpt::time_loaded:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->timeLoadedAboveThreshold);
        break;

    case VisOpt::region:
        visualized_values->SetNumberOfValues(elems->size());
        for(std::size_t i=0;i<elems->size();i++) visualized_values->SetValue(i, (*elems)[i]->region%40);
        break;

    case VisOpt::fracture_support:
        for(icy::Node* nd : *nodes) {
            double value;
            if(nd->crack_tip) value = 3;
            else if(nd->support_node) value = 2;
            else if(nd->reset_timing) value = 1;
            else value = 0;
            visualized_values->SetValue(nd->locId, value);
        }
        break;

    case VisOpt::AbsMx:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, std::abs(nd->str_b[0]));
        break;

    case VisOpt::Mx:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_b[0]);
        break;

    case VisOpt::My:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_b[1]);
        break;

    case VisOpt::Mxy:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_b[2]);
        break;

    case VisOpt::Mx_e:
        visualized_values->SetNumberOfValues(elems->size());
        for(std::size_t i=0;i<elems->size();i++) visualized_values->SetValue(i, (*elems)[i]->a_str_b[0]);
        break;

    case VisOpt::My_e:
        visualized_values->SetNumberOfValues(elems->size());
        for(std::size_t i=0;i<elems->size();i++) visualized_values->SetValue(i, (*elems)[i]->a_str_b[1]);
        break;

    case VisOpt::Mxy_e:
        visualized_values->SetNumberOfValues(elems->size());
        for(std::size_t i=0;i<elems->size();i++) visualized_values->SetValue(i, (*elems)[i]->a_str_b[2]);
        break;

    case VisOpt::Tx:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_m[0]);
        break;

    case VisOpt::Ty:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_m[1]);
        break;

    case VisOpt::Txy:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_m[2]);
        break;

    case VisOpt::Qx:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_s[0]);
        break;

    case VisOpt::Qy:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_s[1]);
        break;

    case VisOpt::stx:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_b_top[0]);
        break;

    case VisOpt::sty:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_b_top[1]);
        break;

    case VisOpt::stxy:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_b_top[2]);
        break;

    case VisOpt::sbx:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_b_bottom[0]);
        break;

    case VisOpt::sby:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_b_bottom[1]);
        break;

    case VisOpt::sbxy:
        for(icy::Node* nd : *nodes) visualized_values->SetValue(nd->locId, nd->str_b_bottom[2]);
        break;

    case VisOpt::stx_e:
        visualized_values->SetNumberOfValues(elems->size());
        for(std::size_t i=0;i<elems->size();i++) visualized_values->SetValue(i, (*elems)[i]->a_str_b_top[0]);
        break;

    default:
        break;
    }

    visualized_values->Modified();

    if(VisualizingVariable == VisOpt::Mxy_e ||
            VisualizingVariable == VisOpt::Mx_e ||
            VisualizingVariable == VisOpt::My_e ||
            VisualizingVariable == VisOpt::region ||
            VisualizingVariable == VisOpt::stx_e)
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

    if(VisualizingVariable == VisOpt::fracture_support) hueLut->SetTableRange(-0.5,5.5);
    else if(VisualizingVariable == VisOpt::region) hueLut->SetTableRange(-0.5,40.5);
    else if(VisualizingVariable == VisOpt::time_loaded) hueLut->SetTableRange(0,temporalThreshold);
    else if(update_minmax)
    {
        double minmax[2];
        visualized_values->GetValueRange(minmax);
        hueLut->SetTableRange(minmax[0], minmax[1]);
    }

    UnsafeUpdateArrows(nodes);
}

void icy::FloeVisualization::InitializeLUT(int table)
{
    const int n = 51;
    hueLut->SetNumberOfTableValues(n);

    if(table==0)
        for ( int i=0; i<n; i++)
                hueLut->SetTableValue(i, (double)lutArrayThermometer[i][0],
                        (double)lutArrayThermometer[i][1],
                        (double)lutArrayThermometer[i][2], 1.0);

    else if(table==1)
    for ( int i=0; i<n; i++)
            hueLut->SetTableValue(i, (double)lutArrayTerrain[i][0],
                    (double)lutArrayTerrain[i][1],
                    (double)lutArrayTerrain[i][2], 1.0);
    else if(table==2)
    for ( int i=0; i<n; i++)
            hueLut->SetTableValue(i, (double)lutArrayTemperature[i][0],
                    (double)lutArrayTemperature[i][1],
                    (double)lutArrayTemperature[i][2], 1.0);
    else if(table==3)
    for ( int i=0; i<n; i++)
            hueLut->SetTableValue(i, (double)lutArrayTemperatureAdj[i][0],
                    (double)lutArrayTemperatureAdj[i][1],
                    (double)lutArrayTemperatureAdj[i][2], 1.0);
    else if(table==4) {
        hueLut->SetNumberOfTableValues(6);
    for ( int i=0; i<6; i++)
            hueLut->SetTableValue(i, (double)lutArrayBands[i][0],
                    (double)lutArrayBands[i][1],
                    (double)lutArrayBands[i][2], 1.0);
    }
    else if(table==5) {
        hueLut->SetNumberOfTableValues(40);
    for ( int i=0; i<40; i++)
            hueLut->SetTableValue(i, (double)lutArrayPastel[i][0],
                    (double)lutArrayPastel[i][1],
                    (double)lutArrayPastel[i][2], 1.0);
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
        nd->PrintoutFan();

        vtkIdType pts2[3];
        std::size_t nFan = nd->fan.size();
        for(std::size_t i=0;i<nFan;i++)
        {
            icy::Node::Sector &f = nd->fan[i];
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
    /*
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
    double pts[3];
    for(std::size_t i=0;i<nNodes;i++)
    {
        icy::Node *nd = (*nodes)[i];
        pts[0]=nd->dir.x();
        pts[1]=nd->dir.y();
        pts[2]=0;
        arrowCoords->SetTuple(i, pts);
        pt=nd->locId;
        cellArray_vertices->InsertNextCell(1, &pt);
    }
    ugrid_vertices->SetCells(VTK_VERTEX, cellArray_vertices);

    ugrid_vertices->GetPointData()->AddArray(arrowCoords);
    ugrid_vertices->GetPointData()->SetActiveVectors("arrowCoords");


    ugrid_vertices->Modified();

//    glyph3D->OrientOn();
//    glyph3D->ScalingOn();
    glyph3D->SetScaleModeToScaleByVector();
//    glyph3D->SetScaleFactor(0.03);
    glyph3D->Modified();
    actor_arrows->GetProperty()->SetColor(colors->GetColor3d("Black").GetData());

*/
}


void icy::FloeVisualization::UnsafeUpdateWaterLine(double totalTime, SimParams &prms)
{
    int mode = prms.loadType;    

    for(unsigned i=0;i<gridSizeX;i++)
        for(unsigned j=0;j<gridSizeY;j++) {
            double x = 25.0 * ((double)i/(double)gridSizeX-0.5);
            double y = 6.0 * ((double)j/(double)gridSizeY-0.5);
            double rest_position = icy::Node::WaterLine(x, y, totalTime, prms);
            points_water->SetPoint(i+j*gridSizeX, x,y,rest_position);
        }

    points_water->Modified();

    actor_sphere->SetVisibility(mode==4);
    if(mode == icy::Model::LoadOpt::indentation)
    {
        sphereSource->SetCenter(0.0, 0.0, 1-totalTime*2.0/100-0.00);
        sphereSource->Modified();
        actor_sphere->VisibilityOn();
    } else actor_sphere->VisibilityOff();
}
