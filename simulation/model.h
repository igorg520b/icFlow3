#ifndef MESHCOLLECTION_H
#define MESHCOLLECTION_H

#include <QFileInfo>
#include <QObject>
#include <QMutex>

#include <vector>
#include <algorithm>
#include <chrono>
#include <unordered_set>
#include <string>

#include "parameters_sim.h"
#include "parameters_beam.h"

#include "mesh.h"
#include "modelcontrollerinterface.h"
#include "linearsystem.h"
#include "floevisualization.h"
#include "serializer.h"

namespace icy { class Model; class Node; class Element;}

class icy::Model : public QObject, public ModelControllerInterface
{
    Q_OBJECT

    // CONTROLLER
public:
    void Reset(unsigned setup);
    void Prepare() override;        // invoked once, at simulation start
    bool Step() override;           // either invoked by Worker or via GUI
    void RequestAbort() override;   // invoked from GUI

    void GoToStep(int step) override;                // only works if serializer has a file open
    void SaveAs(std::string fileName) override;
    void Load(std::string fileName) override;

private:
    Serializer serializer;
    bool abortRequested = false;
    void Aborting();       // called before exiting Step() if aborted
    constexpr static int colWidth = 12;    // table column width when logging

signals:
    void stepCompleted();
    void stepAborted();
    void fractureProgress();    // signals that the VTK view of mesh is not in sync with internal representation


    // MODEL
public:
    SimParams prms;
    BeamParams prms_beam;
    icy::Mesh mesh;
    icy::FloeVisualization vtk_representation;
    icy::LinearSystem ls;

    int currentStep;
    double timeStepFactor, simulationTime;

    // loading options
    enum LoadOpt { stretch_x, stretch_xy, indentation, waterfall, waves_x, waves_xy, waves_diag, waves_wind, L_beam};
    Q_ENUM(LoadOpt)

    void SetIndenterPosition(double position);

private:
    void Fracture(double timeStep);
    void Fracture_LocalSubstep();    // part of Fracture()
    void InitialGuess(double timeStep, double timeStepFactor);
    bool AssembleAndSolve(double timeStep, bool enable_collisions, bool enable_spring,
                          std::vector<icy::Node*> &nodes, std::vector<icy::Element*> &elems);  // return true if solved
    bool AcceptTentativeValues(double timeStep);    // return true if plastic deformation occurred
    void IdentifyDisconnectedRegions();

    // Visualization
public:
    void UnsafeSynchronizeVTK();    // synchronize what VTK shows with internal mesh representation; invoke from the main thread
    void ChangeVisualizationOption(int option);  // invoke from the main thread
    bool topologyInvalid = true;
    bool displacementsInvalid = true;

private:
    QMutex vtk_update_mutex; // to prevent modifying mesh data while updating VTK representation
    bool vtk_update_requested = false;  // true when signal has been already emitted to update vtk geometry

private:
    void UnsafeUpdateGeometry(double simulationTime, SimParams &prms); // to be called from the main thread
    void RestoreFromSerializationBuffers(SimParams &prms); // called from controller after loading data from a file


};

#endif // MESHCOLLECTION_H


/*
    std::vector<icy::FrameInfo> stepStats;  // information about each time step
    icy::FrameInfo ts;  // tentative step info (the step itself may fail)
    bool requestToStop = false; // ModelController request backgroundworker to stop calling Step()

    void ImportFloePatch(QString fileName);
    void Remesh();
    void Reset();                       // reset simulation to pristine state; erase model
    void Trim();                        // remove subsequent steps
    void Prepare();                     // compute constant matrices - call once before first Step()
    void Step();                        // perform one computation step
    void Fracture();
    void RequestAbort();                // cancel current step; invoked by GUI thread

    // progress summary
    int getTotalSteps() {
        std::size_t n = stepStats.size();
        if(n>0) n--;
        return n;
    }
    int getCurrentStep() { return current_step; }

private:
    unsigned current_step = 0;
    bool abort_requested = false;
    void Aborting();       // called before exiting Step() if aborted

    // parts of Step();
    void _BeginStep();  // make initial guess, initialize FrameInfo
    void _Write();

signals:
    void stepCompleted();
    void fractureUpdated();
    void stepAborted();
    void progressUpdated();
*/
