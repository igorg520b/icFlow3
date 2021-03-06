#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

#include <QFileDialog>
#include <QSizePolicy>
#include <QPushButton>
#include <QSplitter>
#include <QLabel>
#include <QVBoxLayout>
#include <QTreeWidget>
#include <QProgressBar>
#include <QMenu>
#include <QList>
#include <QDebug>
#include <QComboBox>
#include <QMetaEnum>
#include <QDir>

#include <QtCharts>
#include <QVTKOpenGLNativeWidget.h>

#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkNew.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkDataSetMapper.h>
#include <vtkCamera.h>

#include <vtkOBJReader.h>
#include <vtkNamedColors.h>
#include <vtkProperty.h>
#include <vtkVersion.h>
#include <vtkWindowToImageFilter.h>
//#include <vtkPointSource.h>
//#include <vtkLineSource.h>
//#include <vtkOBBTree.h>
#include <vtkPolyDataMapper.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkCompositePolyDataMapper2.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkAxesActor.h>

#include <vtkDataSetSurfaceFilter.h>
#include <vtkScalarBarActor.h>
#include <vtkTextProperty.h>

#include <vtkAreaPicker.h>
#include <vtkPointPicker.h>
#include <vtkProp3DCollection.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCallbackCommand.h>
#include <vtkInteractorStyleRubberBandPick.h>

#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkJPEGWriter.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <stdexcept>
#include <cstdint>

#include "objectpropertybrowser.h"
#include "preferences_gui.h"

#include "modelcontroller.h"
#include "backgroundworker.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT
private:
    Ui::MainWindow *ui;

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

    void showEvent( QShowEvent* event ) override;
    void closeEvent( QCloseEvent* event ) override;

    static void PickCallbackFunction(vtkObject* caller,
                              long unsigned int vtkNotUsed(eventId),
                              void* vtkNotUsed(clientData),
                              void* vtkNotUsed(callData));


private slots:

    void sliderValueChanged(int val);
    void comboboxIndexChanged_visualizations(int index);
    void comboboxIndexChanged_load_types(int index);

    void customMenuRequested(QPoint pos);
    void treeItemSelected();
    void progress_updated();
    void render_results();
    void background_worker_paused();

    void on_action_quit_triggered();

    void on_action_scene_open_triggered();
    void on_action_scene_save_as_triggered();
    void on_action_scene_reset_triggered();

    void on_action_import_3D_boundary_triggered();
    void on_action_import_floes_triggered();

    void on_action_floe_remesh_triggered();

    void on_action_simulation_single_step_triggered();

    void on_action_camera_reset_triggered();
    void on_action_show_axes_triggered(bool checked);
    void on_action_show_benchmark_triggered();
    void on_action_show_model_triggered();

    void on_action_show_scalar_bar_triggered(bool checked);
    void on_action_Trim_triggered();

    void on_actionMohr_s_triggered();

    void on_action_Tentative_triggered(bool checked);

    void on_action_draw_Edges_triggered(bool checked);

    void on_action_draw_Boundary_triggered(bool checked);

    void on_action_draw_Arrows_triggered(bool checked);

    void on_action_simulation_start_triggered(bool checked);

    void updateGUI();   // when simulation is started/stopped or when a step is advanced


    void on_action_draw_water_Level_triggered(bool checked);

    void on_action_Screenshot_triggered();

    void on_action_GotoStep0_triggered();

    void on_actionScreenshot_every_step_triggered();

private:
    PreferencesGUI prefsGUI;
    QDir outputDirectory;
    icy::ModelController controller;   // simulation algorithms
    BackgroundWorker *worker;

    QString m_sSettingsFile = "ic_config";
    QLabel *statusLabel;                    // statusbar
    QLabel *statusLabelAttempt;             // attempt # at advancing a time step
    QLabel *statusLabelIteration;           // Newton-Raphson iteration #
    QLabel *statusLabelTimestep;            // time step for current attempt
    QLabel *statusLabelSimTime;             // current simulated time
    QProgressBar *progressBar;      // statusbar
    QSlider *slider;                // in toolbar
    QComboBox *comboBox_visualizations;
    QComboBox *comboBox_load_types;
    QLabel *labelStepCount;

    // splitter and the right side window
    QSplitter *splitter_left_panel;
    QSplitter *splitter_main;
    QHBoxLayout *right_side_layout;
    QWidget *right_side_container;

    // properties
    ObjectPropertyBrowser *pbrowser;

    // charts_pie
    QChartView *chartView_benchmark;
    QChart *chart_pie_benchmark;
    QPieSeries *series_pie_benchmark;

    // series
    QChartView *chartView;
    QChart *chart_line_mohr;
    QLineSeries *series_mohr;
    QScatterSeries *mohr_sectors;

    // VTK
    QVTKOpenGLNativeWidget *qt_vtk_widget;

    vtkNew<vtkGenericOpenGLRenderWindow> renderWindow;
    vtkNew<vtkRenderer> renderer;
    vtkNew<vtkNamedColors> colors;
    vtkNew<vtkAxesActor> axes;
    vtkNew<vtkScalarBarActor> scalarBar;
    vtkNew<vtkPointPicker> pointPicker;
    vtkNew<vtkWindowToImageFilter> windowToImageFilter;
    vtkNew<vtkPNGWriter> writer;
    vtkNew<vtkJPEGWriter> writer2;

    void Reset();
    void OpenSceneFile(QString fileName);

    // tree widget and context menus
    QTreeWidget *tree;
    QMenu *menuSolidGroup;      // for solids group
    QMenu *menuSolidIndividual;
    QTreeWidgetItem *tiParams_sim;
    QTreeWidgetItem *tiFloes;
//    QTreeWidgetItem *tiSolids;

    const QString wtitle = "icFlow3: Finite Element Simulation of Ice - ";

    void save_fan_data_for_testing();
    unsigned screenshot_counter = 0; // for recording video
    void TakeScreenshotForVideo();
};
#endif // MAINWINDOW_H
