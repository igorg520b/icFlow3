#ifndef MESHCOLLECTION_H
#define MESHCOLLECTION_H

#include <QFileInfo>
#include <QObject>
#include <QMutex>

#include <vector>
#include <algorithm>
#include <chrono>
#include <unordered_set>

#include <gmsh.h>

#include "parameters_sim.h"
#include "geometry.h"
#include "linearsystem.h"
#include "solid3d.h"
#include "floevisualization.h"

namespace icy { class Model; class Node; class Element;}

class icy::Model : public QObject
{
    Q_OBJECT

    Q_PROPERTY(int in_Elems READ getElemCount)
    Q_PROPERTY(int in_Nodes READ getNodeCount)
    Q_PROPERTY(int in_Regions READ getRegionCount)
    Q_PROPERTY(double in_length READ getLength)
    Q_PROPERTY(double in_width READ getWidth)
    Q_PROPERTY(double in_area READ getArea)

public:    
    // visualization options
    enum LoadOpt { stretch_x, stretch_xy, indentation, waterfall, waves_x, waves_xy, waves_diag, waves_wind};
    Q_ENUM(LoadOpt)

    void Reset();   // erases all geometry

    icy::Geometry floes;
    icy::FloeVisualization floes_vtk;
    std::vector<Solid3D*> solids;

    void InitialGuessTentativeVals(double timeStep, double beta, double gamma);

    // moved from controller
    void AssembleAndSolve(long &time_clear, long &time_forces, long &time_structure,
                          long &time_assemble, long &time_solve, long &time_pull,
                          SimParams &prms, double timeStep, double totalTime, double &resultSqNorm);

    long PullFromLinearSystem(double timeStep, double beta, double gamma);
    void AcceptTentativeValues(SimParams &prms);
    void FractureStep(SimParams &prms, double timeStep, double totalTime, long &b_substep, long &b_directions, long &b_split, long &b_support);
    long IdentifyDisconnectedRegions();
    void UnsafeUpdateGeometry(double simulationTime, SimParams &prms); // to be called from the main thread
    void RestoreFromSerializationBuffers(SimParams &prms); // called from controller after loading data from a file
    QMutex mutex; // to prevent modifying mesh data while updating VTK representation

private:
    icy::LinearSystem ls;
    long ComputeElasticForcesAndAssemble(SimParams &prms, double timeStep, double totalTime);
    long LocalSubstep(SimParams &prms, double timeStep, double totalTime);

    // synchronize VTK visualization and internal representation
    bool topologyInvalid = true;        // topology changed (means that displacements and values also changed)
    bool displacementsInvalid = true;   // displacements changed since VTK last updated
    bool valuesInvalid = true;          // visualization changed (e.g. via GUI comobobox)

    // signal has been sent to the main thread asking to invoke UnsafeUpdateGeometry()
    // this is to prevent emitting multiple requests before the existing request is processed
    bool updateRequested = false;

    double getElemCount() { return floes.elems->size(); }
    double getNodeCount() { return floes.nodes->size(); }
    double getRegionCount() {return floes.regions.size(); }
    double getLength() {return floes.length; }
    double getWidth() {return floes.width; }
    double getArea() {return floes.area; }

signals:
    void requestGeometryUpdate(); // this goes to the main thread, which calls UnsafeUpdateGeometry()
    void propertyChanged();
};

#endif // MESHCOLLECTION_H
