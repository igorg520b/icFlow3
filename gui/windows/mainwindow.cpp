#include <algorithm>
#include "mainwindow.h"
#include "./ui_mainwindow.h"

MainWindow::~MainWindow() {delete ui;}

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    worker = new BackgroundWorker(&controller);

    connect(&controller.model, SIGNAL(requestGeometryUpdate()), SLOT(render_results()));
    connect(&controller, SIGNAL(stepAborted()),SLOT(updateGUI()));
    connect(&controller, SIGNAL(progressUpdated()),SLOT(progress_updated()));
    connect(&controller, SIGNAL(stepCompleted()), SLOT(updateGUI()));
    connect(worker, SIGNAL(workerPaused()), SLOT(background_worker_paused()));

    // chart, Mohr's circle
    chart_line_mohr = new QChart;
    series_mohr = new QLineSeries;
    series_mohr->setName("ARCSim");
    chart_line_mohr->addSeries(series_mohr);

    mohr_sectors = new QScatterSeries;
    mohr_sectors->setName("sectors");
    chart_line_mohr->addSeries(mohr_sectors);
    chart_line_mohr->createDefaultAxes();

    chartView = new QChartView;
    chartView->setRenderHint(QPainter::Antialiasing);
    chartView->setChart(chart_line_mohr);
    chartView->hide();

    // pie chart for benchmarking
    chartView_benchmark = new QChartView;
    chartView_benchmark->setRenderHint(QPainter::Antialiasing);

    series_pie_benchmark = new QPieSeries;
    series_pie_benchmark->setHoleSize(0.35);
//    series_pie_benchmark->setPieSize(0.5);

    chart_pie_benchmark = new QChart;
    chart_pie_benchmark->addSeries(series_pie_benchmark);
    chart_pie_benchmark->legend()->setAlignment(Qt::AlignBottom);
    //chart_pie_benchmark->setTheme(QChart::ChartThemeBlueCerulean);
    chart_pie_benchmark->legend()->setFont(QFont("Arial", 14));
    chart_pie_benchmark->setFont(QFont("Arial", 16));
    chart_pie_benchmark->setBackgroundRoundness(0);

    chartView_benchmark->setChart(chart_pie_benchmark);
    chartView_benchmark->setContentsMargins(0,0,0,0);
    chartView_benchmark->hide();

    // tree
    tree = new QTreeWidget;
    tree->setHeaderHidden(true);
    tree->setContextMenuPolicy(Qt::ContextMenuPolicy::CustomContextMenu);
    connect(tree, SIGNAL(customContextMenuRequested(QPoint)), SLOT(customMenuRequested(QPoint)));
    connect(tree, SIGNAL(itemSelectionChanged()), SLOT(treeItemSelected()));

    // property browser
    pbrowser = new ObjectPropertyBrowser(this);
    splitter_left_panel = new QSplitter(Qt::Orientation::Vertical);
    splitter_left_panel->addWidget(tree);
    splitter_left_panel->addWidget(pbrowser);

    tiParams_sim = new QTreeWidgetItem(tree);
    tiParams_sim->setText(0, "sim");

    tiFloes = new QTreeWidgetItem(tree);
    tiFloes->setText(0, "floes");

    tiFloes->setExpanded(true);
//    tiSolids = new QTreeWidgetItem(tree);
//    tiSolids->setText(0, "solids");
//    tiSolids->setExpanded(true);

    // VTK
    qt_vtk_widget = new QVTKOpenGLNativeWidget();
    qt_vtk_widget->setRenderWindow(renderWindow);
    renderer->SetBackground(colors->GetColor3d("White").GetData());
    renderer->AddActor(controller.model.floes_vtk.actor_mesh);
    renderer->AddActor(controller.model.floes_vtk.actor_boundary);
    renderer->AddActor(controller.model.floes_vtk.actor_labels);
    renderer->AddActor(controller.model.floes_vtk.actor_water);
    renderer->AddActor(controller.model.floes_vtk.actor_sphere);
    renderer->AddActor(controller.model.floes_vtk.actor_arrows);

    renderer->AddActor(axes);
    renderer->AddActor(scalarBar);
    scalarBar->SetLookupTable(controller.model.floes_vtk.hueLut);

    scalarBar->SetMaximumWidthInPixels(220);
    scalarBar->SetBarRatio(0.07);
    scalarBar->SetMaximumHeightInPixels(800);
    scalarBar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedDisplay();
    scalarBar->GetPositionCoordinate()->SetValue(0.86,0.03, 0.0);
    scalarBar->SetLabelFormat("%.1e");
    scalarBar->GetLabelTextProperty()->BoldOff();
    scalarBar->GetLabelTextProperty()->ItalicOff();
    scalarBar->GetLabelTextProperty()->ShadowOff();
    scalarBar->GetLabelTextProperty()->SetColor(0.1,0.1,0.1);
    scalarBar->GetLabelTextProperty()->SetFontFamilyToTimes();

    renderWindow->AddRenderer(renderer);

    renderWindow->GetInteractor()->SetPicker(pointPicker);
    vtkSmartPointer<vtkCallbackCommand> pickCallback =
            vtkSmartPointer<vtkCallbackCommand>::New();
    pickCallback->SetCallback(MainWindow::PickCallbackFunction);
    pointPicker->AddObserver(vtkCommand::EndPickEvent, pickCallback);
    pickCallback->SetClientData((void*)this);

    // VTK - screenshot
    windowToImageFilter->SetInput(renderWindow);
    windowToImageFilter->SetScale(1); // image quality
    windowToImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
    windowToImageFilter->ReadFrontBufferOn(); // read from the back buffer
    writer->SetInputConnection(windowToImageFilter->GetOutputPort());

    // right frame
    right_side_container = new QWidget;
    right_side_layout = new QHBoxLayout;
    right_side_container->setLayout(right_side_layout);
    right_side_layout->setSpacing(0);
    right_side_layout->setMargin(0);
    right_side_layout->addWidget(qt_vtk_widget);
    right_side_layout->addWidget(chartView_benchmark, 1);
    right_side_layout->addWidget(chartView, 1);

    // splitter
    splitter_main = new QSplitter(Qt::Orientation::Horizontal);
    splitter_main->addWidget(splitter_left_panel);
    splitter_main->addWidget(right_side_container);
    setCentralWidget(splitter_main);

    // toolbar - comboboxes
    comboBox_visualizations = new QComboBox();
    comboBox_load_types = new QComboBox();
    ui->toolBar->addWidget(comboBox_visualizations);
    ui->toolBar->addWidget(comboBox_load_types);

    // populate combobox
    QMetaEnum qme = QMetaEnum::fromType<icy::FloeVisualization::VisOpt>();
    for(int i=0;i<qme.keyCount();i++) comboBox_visualizations->addItem(qme.key(i));

    qme = QMetaEnum::fromType<icy::Model::LoadOpt>();
    for(int i=0;i<qme.keyCount();i++) comboBox_load_types->addItem(qme.key(i));

    connect(comboBox_visualizations, QOverload<int>::of(&QComboBox::currentIndexChanged),
            [=](int index){ comboboxIndexChanged_visualizations(index); });
    connect(comboBox_load_types, QOverload<int>::of(&QComboBox::currentIndexChanged),
            [=](int index){ comboboxIndexChanged_load_types(index); });

    // slider
    slider = new QSlider(Qt::Horizontal);
    ui->toolBar->addWidget(slider);
    labelStepCount = new QLabel();
    ui->toolBar->addWidget(labelStepCount);

    connect(slider, SIGNAL(valueChanged(int)), this, SLOT(sliderValueChanged(int)));

    // statusbar
    statusLabel = new QLabel("-");
    progressBar = new QProgressBar(this);
    progressBar->hide();

    statusLabelAttempt = new QLabel("Att: ");
    statusLabelIteration = new QLabel("Iter: ");
    statusLabelTimestep = new QLabel("Ts: ");
    statusLabelSimTime = new QLabel("ST: ");

    QSizePolicy sp;
    sp.setHorizontalPolicy(QSizePolicy::Fixed);
    progressBar->setSizePolicy(sp);
    statusLabelAttempt->setSizePolicy(sp);
    statusLabelIteration->setSizePolicy(sp);
    statusLabelTimestep->setSizePolicy(sp);
    statusLabelAttempt->setFixedWidth(70);
    statusLabelIteration->setFixedWidth(70);
    statusLabelTimestep->setFixedWidth(100);
    statusLabelSimTime->setSizePolicy(sp);
    statusLabelSimTime->setFixedWidth(100);
    progressBar->setFixedWidth(200);
    progressBar->setFixedHeight(25);

    ui->statusbar->addWidget(statusLabel);
    ui->statusbar->addPermanentWidget(progressBar,0);
    ui->statusbar->addPermanentWidget(statusLabelAttempt);
    ui->statusbar->addPermanentWidget(statusLabelIteration);
    ui->statusbar->addPermanentWidget(statusLabelTimestep);
    ui->statusbar->addPermanentWidget(statusLabelSimTime);


    // read/restore saved settings
    QSettings settings(m_sSettingsFile);
    splitter_left_panel->restoreState(settings.value("splitter_left_panel").toByteArray());
    splitter_main->restoreState(settings.value("splitter_main").toByteArray());

    QVariant var = settings.value("camData");
    if(!var.isNull()) {
        QByteArray arr = var.toByteArray();
        double *vec = (double*)arr.constData();
        renderer->GetActiveCamera()->SetPosition(vec[0],vec[1],vec[2]);
        renderer->GetActiveCamera()->SetFocalPoint(vec[3],vec[4],vec[5]);
        renderer->GetActiveCamera()->SetViewUp(vec[6],vec[7],vec[8]);
        renderer->GetActiveCamera()->SetViewAngle(vec[9]);
    }

    prefsGUI.LoadState(settings);
    ui->action_Load_Recent_on_Startup->setChecked(prefsGUI.LoadLastScene);

    if(prefsGUI.LoadLastScene) {
        QFile sceneFile(prefsGUI.LastSceneFilename);
        if(sceneFile.exists()) OpenSceneFile(prefsGUI.LastSceneFilename);
    }

    // create context menus that are used with tree widget
    menuSolidGroup = new QMenu(this);
    menuSolidGroup->addAction(ui->action_import_3D_boundary);

    menuSolidIndividual = new QMenu(this);
    menuSolidIndividual->addAction(ui->action_mesh_remove);

    // create/initialize outputDirectory
    QString outputPathName = "output";
    outputDirectory.setPath(outputPathName);
    if(!outputDirectory.exists()) {
        bool result = QDir::current().mkdir(outputPathName);
        if(!result) throw std::runtime_error("could not create output directory");
        outputDirectory.setPath(QDir::current().path() + "/output");
        if(!outputDirectory.exists()) throw std::runtime_error("could not open output directory");
    }
    qDebug() << "output directory: " << outputDirectory.path();

    ui->action_show_axes->setChecked(prefsGUI.ShowAxes);
    on_action_show_axes_triggered(prefsGUI.ShowAxes);

    ui->action_Tentative->setChecked(prefsGUI.ShowTentative);
    on_action_Tentative_triggered(prefsGUI.ShowTentative);
    ui->action_Tentative->setVisible(false);

    ui->action_show_scalar_bar->setChecked(prefsGUI.ShowScalarBar);
    on_action_show_scalar_bar_triggered(prefsGUI.ShowScalarBar);

    ui->action_draw_Edges->setChecked(prefsGUI.DrawEdges);
    on_action_draw_Edges_triggered(prefsGUI.DrawEdges);

    ui->action_draw_Boundary->setChecked(prefsGUI.DrawBoundary);
    on_action_draw_Boundary_triggered(prefsGUI.DrawBoundary);

    ui->action_draw_Arrows->setChecked(prefsGUI.DrawArrows);
    on_action_draw_Arrows_triggered(prefsGUI.DrawArrows);

    ui->action_draw_water_Level->setChecked(prefsGUI.ShowWaterLevel);
    on_action_draw_water_Level_triggered(prefsGUI.ShowWaterLevel);

    comboBox_visualizations->setCurrentIndex(prefsGUI.VisualizationOption);
}

void MainWindow::showEvent( QShowEvent*)
{
    qDebug() << "size of element " << sizeof(icy::Element);
    qDebug() << "size of node " << sizeof(icy::Node);
    qDebug() << "size of edge " << sizeof(icy::Edge);
    pbrowser->setActiveObject(&controller.prms);
    updateGUI();
    renderWindow->Render();
}

// present context menu for tree widget
void MainWindow::customMenuRequested(QPoint pos)
{
    qDebug() << "context menu requested";
    QTreeWidgetItem *nd = tree->itemAt(pos);
    QVariant qv1 = nd->data(1, Qt::UserRole);

    //if(nd == tiSolids) menuSolidGroup->exec(tree->mapToGlobal(pos));
    if(qv1.toInt()==2) menuSolidIndividual->exec(tree->mapToGlobal(pos));
}

void MainWindow::on_action_scene_open_triggered()
{
    if(worker->running) return; // ensure that backgroundworker is stopped
    qDebug() << "select and open a scene file";
    QString fileName = QFileDialog::getOpenFileName(this, "Open Scene",
                                                    outputDirectory.path(),
                                                    "Scene Files (*.h5)");
    if (fileName.isEmpty()) return;

    Reset();
    prefsGUI.LastSceneFilename = fileName;
    OpenSceneFile(fileName);
}

void MainWindow::OpenSceneFile(QString fileName)
{
    if(worker->running) return;

    controller.Load(fileName);
    comboBox_load_types->setCurrentIndex(controller.prms.loadType);
    pbrowser->setActiveObject(&controller.prms);

    // update window title
    QFileInfo fi(fileName);
    setWindowTitle(wtitle + fi.baseName());
    renderWindow->Render();
    updateGUI();
}

void MainWindow::on_action_scene_save_as_triggered()
{
    if(worker->running) return;
    qDebug() << "menu action: save scene as...";

    QFileDialog qfd(this, "Save Scene", outputDirectory.path(), "Scene Files (*.h5)");
    qfd.setDefaultSuffix("h5");

    qfd.setAcceptMode(QFileDialog::AcceptSave);
    qfd.setFileMode(QFileDialog::AnyFile);

    bool result = qfd.exec();
    if(!result) return;

    QString fileName = qfd.selectedFiles()[0];

    if (fileName.isEmpty()) return;

    QFileInfo fileInfo(fileName);
    prefsGUI.LastSceneFilename = fileName;

    // save as new file
    controller.SaveAs(fileName);
    setWindowTitle(wtitle + fileInfo.baseName());
    updateGUI();    // current step becomes 0
}

void MainWindow::on_action_scene_reset_triggered() { Reset(); }

void MainWindow::on_action_import_3D_boundary_triggered()
{
    // this functionality is currently unused
    qDebug() << "STL import enclosed boundary for 3D object";
    QString fileName = QFileDialog::getOpenFileName(this, "Open 2D STL Mesh",
                                         prefsGUI.LastFolder3DGometry, "Mesh Files (*.stl)");
    if (fileName.isEmpty()) return;
    QFileInfo fileInfo(fileName);
    prefsGUI.LastFolder3DGometry = fileInfo.dir().path();

    icy::Solid3D *solid = new icy::Solid3D(fileName);
    controller.model.solids.push_back(solid);
//    tiSolids->addChild(&solid->treeWidget);
    renderer->AddActor(solid->actor_mesh);
    renderWindow->Render();
}

void MainWindow::on_action_import_floes_triggered()
{
    qDebug() << "on_actionImport_Floe_Patch_from_STL_triggered()";
    if(worker->running) return;

    QString fileName = QFileDialog::getOpenFileName(this, "Open 2D STL Patch",
                                                    prefsGUI.LastFolderFloes,
                                                    "Mesh Files (*.stl)");
    if (fileName.isEmpty()) return;
    QFileInfo fileInfo(fileName);
    prefsGUI.LastFolderFloes = fileInfo.dir().path();

    // IMPORT FLOE PATCH
    controller.ImportFloePatch(fileName);
    controller.model.floes_vtk.UnsafeUpdateTopology(controller.model.floes.nodes.get(),
                                                    controller.model.floes.elems.get(),
                                                    controller.model.floes.boundaryEdges,
                                                    controller.prms.temporal_attenuation);
    updateGUI();
    renderWindow->Render();
}

void MainWindow::on_action_floe_remesh_triggered()
{
    qDebug() << "on_action_floe_Remesh_triggered()";
    if(worker->running) return;
    controller.Remesh();
    renderWindow->Render();
}

void MainWindow::treeItemSelected()
{
    qDebug() << "tree selection changed";
    QList<QTreeWidgetItem*> sl = tree->selectedItems();
    if(sl.size() == 1) {
        pbrowser->setActiveObject(nullptr);
        QTreeWidgetItem* w = sl.first();
        // set property browser
        if(w == tiParams_sim) pbrowser->setActiveObject(&controller.prms);
        else if(w == tiFloes) pbrowser->setActiveObject(&controller.model);
        else {
            QVariant qv0 = w->data(0, Qt::UserRole);
            QVariant qv1 = w->data(1, Qt::UserRole);
            if(qv1.toInt()==2)
                pbrowser->setActiveObject(qvariant_cast<QObject *>(qv0));
        }
    }
}

void MainWindow::on_action_quit_triggered() { this->close(); }

void MainWindow::closeEvent( QCloseEvent* event )
{
    // save settings and stop simulation
    QSettings settings(m_sSettingsFile);
    qDebug() << "MainWindow: closing main window; " << settings.fileName();

    settings.setValue("splitter_left_panel", splitter_left_panel->saveState());
    settings.setValue("splitter_main", splitter_main->saveState());

    double data[10];
    renderer->GetActiveCamera()->GetPosition(&data[0]);
    renderer->GetActiveCamera()->GetFocalPoint(&data[3]);
    renderer->GetActiveCamera()->GetViewUp(&data[6]);
    data[9]=renderer->GetActiveCamera()->GetViewAngle();

    QByteArray arr((char*)&data[0], sizeof(double)*10);
    settings.setValue("camData", arr);

    prefsGUI.LoadLastScene = ui->action_Load_Recent_on_Startup->isChecked();
    prefsGUI.SaveState(settings);

    // kill backgroundworker
    worker->Finalize();

    for(auto &solid : controller.model.solids) {
        renderer->RemoveActor(solid->actor_mesh);
        tree->removeItemWidget(&solid->treeWidget,0);
    }

    controller.prms.Serialize();
    controller.serializer.SaveParams(controller.prms.serialization_buffer,
                                     icy::SimParams::buffer_size);
    controller.serializer.CloseFile();
    event->accept();
}

void MainWindow::background_worker_paused()
{
    ui->action_simulation_start->blockSignals(true);
    ui->action_simulation_start->setEnabled(true);
    ui->action_simulation_start->setChecked(false);
    ui->action_simulation_start->blockSignals(false);
}

void MainWindow::render_results()
{
    controller.model.UnsafeUpdateGeometry(controller.ts.SimulationTime, controller.prms);
    renderWindow->Render();
}

void MainWindow::Reset()
{
    qDebug() << "MainWindow::Reset()";
    if(worker->running) return; // ensure that backgroundworker is stopped
    pbrowser->setActiveObject(nullptr);

    for(auto &solid : controller.model.solids) {
        renderer->RemoveActor(solid->actor_mesh);
        tree->removeItemWidget(&solid->treeWidget,0);
    }

    controller.Reset();
    comboBox_load_types->setCurrentIndex(controller.prms.loadType);
    prefsGUI.LastSceneFilename = "";
    renderWindow->Render();
    updateGUI();
    setWindowTitle(wtitle);
    qDebug() << "finished MainWindow::Reset()";
}

void MainWindow::on_action_simulation_single_step_triggered()
{
    qDebug() << "take one step";
    controller.Prepare();
    controller.Step();
}

void MainWindow::on_action_camera_reset_triggered()
{
    vtkCamera* camera = renderer->GetActiveCamera();
    camera->SetClippingRange(1e1,1e3);
    camera->SetFocalPoint(0.0, 0.0, 0.0);
    camera->SetPosition(0.0, 0.0, 50.0);
    camera->SetViewUp(0.0, 1.0, 0.0);
    camera->Modified();
    renderWindow->Render();
}

void MainWindow::sliderValueChanged(int val)
{
    controller.GoToStep(val);
    controller.model.floes_vtk.UnsafeUpdateTopology(controller.model.floes.nodes.get(),
                                                    controller.model.floes.elems.get(),
                                                    controller.model.floes.boundaryEdges,
                                                    controller.prms.temporal_attenuation);
    renderWindow->Render();
    progress_updated();
}

void MainWindow::on_action_GotoStep0_triggered() { slider->setValue(0); }

void MainWindow::updateGUI()
{
    // scrollbar
    slider->blockSignals(true);

    unsigned max_range = controller.getTotalSteps();
    slider->setRange(0, max_range);
    slider->setValue(controller.getCurrentStep());
    slider->setEnabled(controller.getTotalSteps() != 0 && !worker->running);

    bool r = worker->running;
    ui->action_simulation_single_step->setEnabled(!r);
    ui->action_floe_remesh->setEnabled(!r);
    ui->action_scene_open->setEnabled(!r);
    ui->action_scene_reset->setEnabled(!r);
    ui->action_Trim->setEnabled(!r);
    ui->action_mesh_remove->setEnabled(!r);
    ui->action_import_floes->setEnabled(!r);
    ui->action_scene_save_as->setEnabled(!r);
    ui->action_import_3D_boundary->setEnabled(!r);
    ui->action_simulation_single_step->setEnabled(!r);
    ui->action_GotoStep0->setEnabled(!r);

    statusLabel->setText(r ? "running" : "paused");
    if(!r) ui->action_simulation_start->setEnabled(true);

    slider->blockSignals(false);
    progress_updated();
}

void MainWindow::progress_updated()
{
    progressBar->setValue(controller.ts.solverProgress);
    icy::FrameInfo &ts = controller.ts;
    statusLabelAttempt->setText(QString{"Att: %1"}.arg(ts.count_attempts));
    statusLabelIteration->setText(QString{"Iter: %1"}.arg(ts.count_iterations));
    statusLabelSimTime->setText(QString{ "ST: %1" }.arg(ts.SimulationTime,6,'f',4));
    statusLabelTimestep->setText(QString{"Ts: %1"}.arg(ts.TimeStep,6,'f',4));

    labelStepCount->setText(QString::number(controller.getCurrentStep()));
}

void MainWindow::on_action_show_axes_triggered(bool checked)
{
    // toggle visibility of axes
    prefsGUI.ShowAxes = checked;
    axes->SetVisibility(checked);
    renderWindow->Render();
}

void MainWindow::on_action_show_benchmark_triggered()
{
    series_pie_benchmark->clear();
    series_pie_benchmark->setPieStartAngle(100);
    series_pie_benchmark->setPieEndAngle(100+360);

    series_pie_benchmark->setName("Motion (per solve)");

    if(controller.stepStats.size() > 0)
    {
        long total, n_solves;
        std::vector<std::pair<std::string, long>> results_motion;
        icy::FrameInfo::BenchmarkingSummarize(controller.stepStats, results_motion, total, n_solves);
        total/=1000000;

        std::cout << "Benchmark printout:\n";
        for(const auto &p : results_motion) {
            long val = p.second;
            if(val > 0) {
                QString str= QString::fromStdString(p.first) + " " + QString::number(val)+ " ms";
                series_pie_benchmark->append(str, val);
                std::cout << "{\"" << p.first << "\"," << val << "},\n";
            }
        }
        std::cout << std::endl;

        QList<QPieSlice*> list1 = series_pie_benchmark->slices();
        foreach(QPieSlice *slice, list1)
        {
            slice->setLabelFont(QFont("Times", 20));
        }

        QString strTotal = QString::number(total);
        QString strSolves = QString::number(n_solves);
        QString strSteps = QString::number(controller.stepStats.size()-1);
        QString str = "Total " + strTotal + " s; solves " + strSolves + "; steps " + strSteps;
        chart_pie_benchmark->setTitle(str);
    } else {
        chart_pie_benchmark->setTitle("No data");
    }
    series_pie_benchmark->setLabelsVisible();

    qt_vtk_widget->hide();
    chartView_benchmark->show();
    chartView->hide();
}

void MainWindow::on_action_show_model_triggered()
{
    qt_vtk_widget->show();
    chartView_benchmark->hide();
    chartView->hide();
    renderWindow->Render();
}

void MainWindow::on_actionMohr_s_triggered()
{
    int idx = controller.model.floes_vtk.selectedPointId;
    if(idx<0) return;

    series_mohr->clear();
    mohr_sectors->clear();

    std::vector<std::tuple<double,double,double>> res_curve;
    double ymin = DBL_MAX, ymax = -DBL_MAX;
    double xmin = DBL_MAX, xmax = -DBL_MAX;
    icy::Node *nd = (*controller.model.floes.nodes)[idx];
    nd->InitializeFan();

    std::cout << std::endl << "{";
    const unsigned nPoints = 500;
    for(unsigned i=0;i<nPoints;i++)
    {
        double angle = nd->fan_angle_span*(double)i/(double)(nPoints-1);

        icy::Node::SepStressResult tmpSsr;
        nd->evaluate_tractions(angle, tmpSsr, controller.prms.weakening_coeff);

        double val_x = tmpSsr.trac_normal_top;
        double val_y = tmpSsr.trac_tangential_top;
        if(std::isnan(val_x)) qDebug() << "trac_normal_top is NaN";
        if(std::isnan(val_y)) qDebug() << "trac_tangential_top is NaN";
        xmin = std::min(xmin, val_x);
        xmax = std::max(xmax, val_x);
        ymin = std::min(ymin, val_y);
        ymax = std::max(ymax, val_y);
        series_mohr->append(val_x,val_y);
        res_curve.push_back(std::make_tuple(angle, val_x, val_y));
        std::cout << "{" << angle << "," << val_x/1000 << "," << val_y/1000 << "}";
        if(i!=nPoints-1) std::cout << ",";
    }
    std::cout << "}" << std::endl;
    double span_y = ymax-ymin;
    double span_x = xmax-xmin;
    double span = std::max(span_y, span_x)*0.65;
    double avgx = (xmax+xmin)/2;

    QList<QAbstractAxis*> axes = chart_line_mohr->axes();
    axes[0]->setRange(avgx-span, avgx+span);
    axes[1]->setRange(-span, span);

    // sectors

    std::vector<std::tuple<double,double,double>> res_sectors;
    std::cout << std::endl << "{";
    for(icy::Node::Sector &s : nd->fan)
    {
        double angle = s.angle1;
        icy::Node::SepStressResult tmpSsr;
        nd->evaluate_tractions(angle, tmpSsr, controller.prms.weakening_coeff);

        double val_x = tmpSsr.trac_normal_top;
        double val_y = tmpSsr.trac_tangential_top;
        mohr_sectors->append(val_x, val_y);
        res_sectors.push_back(std::make_tuple(angle, val_x, val_y));
        std::cout << "{" << angle << "," << val_x/1000 << "," << val_y/1000 << "},";
    }
    icy::Node::SepStressResult tmpSsr;
    nd->evaluate_tractions(0, tmpSsr, controller.prms.weakening_coeff);
    double val_x = tmpSsr.trac_normal_top;
    double val_y = tmpSsr.trac_tangential_top;
    mohr_sectors->append(val_x, val_y);
    res_sectors.push_back(std::make_tuple(0, val_x, val_y));
    std::cout << "{" << 0 << "," << val_x/1000 << "," << val_y/1000 << "}";
    std::cout << "}" << std::endl;

    /*
    // export fan data for testing
    std::cout << "fan:" << std::endl;
    std::size_t nFan = nd->fan.size();
    for (std::size_t f=0; f < nFan; f++) {
        icy::Node::Sector &fp = nd->fan[f];
        Eigen::Vector3d u = fp.e[0]->getVec(nd);
        Eigen::Vector3d v = fp.e[1]->getVec(nd);
        std::cout << "model->AddSector(" << u.x() << "," << u.y() << ",";
        std::cout << v.x() << "," << v.y() << ",";
        double sx = fp.face->str_b_bottom.coeff(0);
        double sy = fp.face->str_b_bottom.coeff(1);
        double txy = fp.face->str_b_bottom.coeff(2);
        std::cout << sx << "," << sy << "," << txy << ");" << '\n';
    }
    std::cout << std::endl;
    */

    chartView->show();
    qt_vtk_widget->hide();
    chartView_benchmark->hide();
}


void MainWindow::comboboxIndexChanged_visualizations(int index)
{
    controller.model.floes_vtk.UnsafeUpdateValues(controller.model.floes.nodes.get(),
                                                  controller.model.floes.elems.get(),
                                                  controller.prms.temporal_attenuation,
                                                  (icy::FloeVisualization::VisOpt)index);
    prefsGUI.VisualizationOption = index;
    renderWindow->Render();
    scalarBar->SetVisibility(prefsGUI.ShowScalarBar && prefsGUI.VisualizationOption!=0);
}

void MainWindow::comboboxIndexChanged_load_types(int index)
{
    controller.prms.loadType = index;
    controller.model.floes_vtk.UnsafeUpdateValues(controller.model.floes.nodes.get(),
                                                  controller.model.floes.elems.get(),
                                                  controller.prms.temporal_attenuation,
                                                  (icy::FloeVisualization::VisOpt)prefsGUI.VisualizationOption);
    renderWindow->Render();
}


void MainWindow::on_action_show_scalar_bar_triggered(bool checked)
{
    prefsGUI.ShowScalarBar = checked;
    scalarBar->SetVisibility(prefsGUI.ShowScalarBar && prefsGUI.VisualizationOption!=0);
    renderWindow->Render();
}

void MainWindow::on_action_Trim_triggered()
{
    // for testing
    controller.Trim();
    renderWindow->Render();
    updateGUI();
}

void MainWindow::PickCallbackFunction(vtkObject* caller,
                          long unsigned int vtkNotUsed(eventId),
                          void* clientData,
                          void* vtkNotUsed(callData))
{
    vtkPointPicker* pp = static_cast<vtkPointPicker*>(caller);
    MainWindow *mw = (MainWindow*)clientData;
    mw->controller.model.floes_vtk.UnsafeUpdateSelection(mw->controller.model.floes.nodes.get(),
                                                     pp->GetPointId());
    mw->renderWindow->Render();
}





void MainWindow::on_action_Tentative_triggered(bool checked)
{
    // this action is hidden from GUI
    prefsGUI.ShowTentative = checked;
    controller.model.floes_vtk.use_tentative_coordinates = checked;
    controller.model.floes_vtk.UnsafeUpdateDisplacements(controller.model.floes.nodes.get(),
                                                         controller.model.floes.elems.get(),
                                                         controller.prms.temporal_attenuation);
    renderWindow->Render();
}

void MainWindow::on_action_draw_Edges_triggered(bool checked)
{
    prefsGUI.DrawEdges = checked;
    controller.model.floes_vtk.actor_mesh->GetProperty()->SetEdgeVisibility(checked);
    renderWindow->Render();
}

void MainWindow::on_action_draw_Boundary_triggered(bool checked)
{
    prefsGUI.DrawBoundary = checked;
    controller.model.floes_vtk.actor_boundary->SetVisibility(checked);
    renderWindow->Render();
}

void MainWindow::on_action_draw_Arrows_triggered(bool checked)
{
    prefsGUI.DrawArrows = checked;
    controller.model.floes_vtk.update_arrows = checked;
    controller.model.floes_vtk.UnsafeUpdateArrows(controller.model.floes.nodes.get());
    renderWindow->Render();
}

void MainWindow::on_action_simulation_start_triggered(bool checked)
{
    if(!worker->running && checked){
        qDebug() << "start button - starting";
        statusLabel->setText("starting simulation");
        controller.Prepare();
        worker->Resume();
    }
    else if(worker->running && !checked)
    {
        statusLabel->setText("pausing simulation");
        qDebug() << "start button - pausing";
        worker->Pause();
        ui->action_simulation_start->setEnabled(false);
    }
}

void MainWindow::on_action_draw_water_Level_triggered(bool checked)
{
    prefsGUI.ShowWaterLevel = checked;
    controller.model.floes_vtk.actor_water->SetVisibility(checked);
    renderWindow->Render();
}

void MainWindow::on_action_Screenshot_triggered()
{
    windowToImageFilter->SetInput(renderWindow);
    windowToImageFilter->SetScale(1); // image quality
    windowToImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
    windowToImageFilter->ReadFrontBufferOn(); // read from the back buffer
    windowToImageFilter->Update();
    windowToImageFilter->Modified();

    writer->Modified();
    QString outputPath = QString::number(controller.getCurrentStep()) + ".png";
    qDebug() << "taking screenshot: " << outputPath;
    writer->SetFileName(outputPath.toUtf8().constData());
    writer->Write();
}
