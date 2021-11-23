#include <vtkPointData.h>
#include <QtGlobal>
#include "model.h"
#include <spdlog/spdlog.h>

namespace model = gmsh::model;
namespace factory = gmsh::model::occ;

void icy::Model::Reset()
{
    qDebug() << "icy::Model::Reset()";
    for(auto &solid : solids) delete solid;
    solids.clear();
    floes.Reset();
    topologyInvalid = displacementsInvalid = valuesInvalid = true;
}

void icy::Model::InitialGuessTentativeVals(double timeStep, double beta, double gamma)
{
    double &h = timeStep;
    double hsq = timeStep*timeStep;
    double c1 = 1.0 - 1.0/(2.0*beta);
    double c2 = 1.0/(h*beta);
    double c3 = 1.0/(hsq*beta);
    double c4 = h-(h*gamma)/(2.0*beta);
    double c5 = 1.0-(gamma/beta);
    double c6 = gamma*c2;

    std::size_t nNodes = floes.getNodeCount();
#pragma omp parallel for
    for(std::size_t i=0;i<nNodes;i++) {
        icy::Node *nd = (*floes.nodes)[i];
        nd->ut = nd->un + nd->vn*h + nd->an*(hsq/2);
        nd->xt = nd->ut;
        nd->xt.x() += nd->x_initial.x();
        nd->xt.y() += nd->x_initial.y();

        // per free node, calculate velocity and acceleration at step n+1 from given xt
        Eigen::Matrix<double,LinearSystem::DOFS,1> xt_xn = nd->ut - nd->un;
        nd->at = nd->an*c1 - nd->vn*c2 + xt_xn*c3;
        nd->vt = nd->an*c4 + nd->vn*c5 + xt_xn*c6;
    }
}

long icy::Model::PullFromLinearSystem(double timeStep, double beta, double gamma)
{
    auto t1 = std::chrono::high_resolution_clock::now();

    double &h = timeStep;
    double hsq = timeStep*timeStep;

    // these constants arise from Newmark-Beta equations
    double c1 = 1.0 - 1.0/(2.0*beta);
    double c2 = 1.0/(h*beta);
    double c3 = 1.0/(hsq*beta);
    double c4 = h-(h*gamma)/(2.0*beta);
    double c5 = 1.0-(gamma/beta);
    double c6 = gamma*c2;

    std::size_t nNodes = floes.getNodeCount();
#pragma omp parallel for
    for(std::size_t i=0;i<nNodes;i++) {
        icy::Node *nd = (*floes.nodes)[i];
        if(nd->lsId < 0) continue;
        ls.AdjustCurrentGuess(nd->lsId, nd->ut);
        nd->xt = nd->ut;
        nd->xt.x() += nd->x_initial.x();
        nd->xt.y() += nd->x_initial.y();
        // per free node, calculate velocity and acceleration at step n+1 from given xt
        Eigen::Matrix<double,LinearSystem::DOFS,1> xt_xn = nd->ut - nd->un;
        nd->at = nd->an*c1 - nd->vn*c2 + xt_xn*c3;
        nd->vt = nd->an*c4 + nd->vn*c5 + xt_xn*c6;
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
}

long icy::Model::ComputeElasticForcesAndAssemble(SimParams &prms, double timeStep, double totalTime)
{
    std::size_t nNodes = floes.nodes->size();
    std::size_t nElems = floes.elems->size();

    auto t1 = std::chrono::high_resolution_clock::now();
#pragma omp parallel for
    for(std::size_t i=0;i<nElems;i++)
        (*floes.elems)[i]->ComputeElasticForce(ls, prms, timeStep, floes.elasticityMatrix, floes.D_mats);

#pragma omp parallel for
    for(std::size_t i=0;i<nNodes;i++)
        (*floes.nodes)[i]->ComputeElasticForce(ls, prms, timeStep, totalTime);

    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
}

void icy::Model::AcceptTentativeValues(SimParams &prms)
{
    floes.EvaluateStresses(prms, *floes.elems);
    mutex.lock();
    std::size_t nNodes = floes.nodes->size();
#pragma omp parallel for
    for(std::size_t i=0;i<nNodes;i++)
    {
        icy::Node *nd = (*floes.nodes)[i];
        nd->AcceptTentativeValues();
        nd->normal_n.setZero(); // allow to accumulate
    }
    floes.DistributeStresses(); // also distribute normals

#pragma omp parallel for
    for(std::size_t i=0;i<nNodes;i++)
        (*floes.nodes)[i]->normal_n.normalize();

    mutex.unlock();
    displacementsInvalid = true;
    if(!updateRequested) { updateRequested = true; emit requestGeometryUpdate();}
}

void icy::Model::AssembleAndSolve(long &time_clear, long &time_forces, long &time_structure,
                                  long &time_assemble, long &time_solve, long &time_pull,
                                  SimParams &prms, double timeStep, double totalTime,
                                  double &resultSqNorm)
{
    time_clear += ls.ClearAndResize(floes.getFreeNodeCount());

    auto t1 = std::chrono::high_resolution_clock::now();
    std::size_t nElems = floes.elems->size();

#pragma omp parallel for
    for(std::size_t i=0;i<nElems;i++) (*floes.elems)[i]->UpdateSparseSystemEntries(ls);
    auto t2 = std::chrono::high_resolution_clock::now();
    time_structure+= std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();

    time_structure += ls.CreateStructure();
    time_forces += ComputeElasticForcesAndAssemble(prms, timeStep, totalTime);
    time_solve += ls.Solve();
    time_pull += PullFromLinearSystem(timeStep, prms.NewmarkBeta, prms.NewmarkGamma);
    resultSqNorm = ls.SqNormOfDx();
}


long icy::Model::LocalSubstep(SimParams &prms, double timeStep, double totalTime)
{
    auto t1 = std::chrono::high_resolution_clock::now();

    double localTimeStep = timeStep*prms.substepping_timestep_factor;
    if(floes.local_support.size() == 0)
    {
        // this may occur in manual "testing" situation
        qDebug() << "LocalSubstep: support size is zero, aborting";
        return 0;
    }

    for(icy::Node *nd : *floes.nodes) nd->lsId=-1;
    int count = 0;
    for(icy::Node *nd : floes.local_support) nd->lsId = count++;

    //support_range1
    // similar to AssembleAndSolve, but only run on the local domain
    std::size_t nNodesLocal = floes.local_support.size();
    std::size_t nElemsLocal = floes.local_elems.size();
    for(int i=0;i<prms.substep_iterations;i++)
    {
        ls.ClearAndResize(count);

#pragma omp parallel for
    for(std::size_t i=0;i<nElemsLocal;i++) floes.local_elems[i]->UpdateSparseSystemEntries(ls);

        ls.CreateStructure();

#pragma omp parallel for
        for(std::size_t i=0;i<nElemsLocal;i++)
            floes.local_elems[i]->ComputeElasticForce(ls, prms, timeStep, floes.elasticityMatrix, floes.D_mats);

#pragma omp parallel for
        for(std::size_t i=0;i<nNodesLocal;i++)
            floes.local_support[i]->ComputeElasticForce(ls, prms, localTimeStep, totalTime);

        ls.Solve();
        PullFromLinearSystem(localTimeStep, prms.NewmarkBeta, prms.NewmarkGamma);
    }

    floes.EvaluateStresses(prms, floes.local_elems);
    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
}

void icy::Model::FractureStep(SimParams &prms, double timeStep, double totalTime,
                              long &b_substep, long &b_compute_fracture_directions, long &b_split, long &b_support)
{
    // ComputeFractureDirections must be invoked prior to this
    if(floes.maxNode == nullptr) throw std::runtime_error("FractureStep");

    mutex.lock();
    b_split += floes.SplitNodeAlt(prms);
    mutex.unlock();
    topologyInvalid = true;
    b_support += floes.InferLocalSupport(prms);

    b_substep += LocalSubstep(prms, timeStep, totalTime);

    mutex.lock();
    b_compute_fracture_directions += floes.ComputeFractureDirections(prms);
    mutex.unlock();

    displacementsInvalid = true;

    if(!updateRequested) {updateRequested = true; emit requestGeometryUpdate(); }
}

void icy::Model::UnsafeUpdateGeometry(double simulationTime, SimParams &prms)
{
    mutex.lock();   // prevent modifying nodes, elems and edges while updating VTK arrays
    updateRequested = false;    // reset flag, so that new requests will be issued

    if(topologyInvalid)
    {
        // re-create mesh topology
        topologyInvalid = valuesInvalid = displacementsInvalid = false; // update displacements and values
        floes_vtk.UnsafeUpdateTopology(floes.nodes.get(), floes.elems.get(), floes.boundaryEdges, prms.temporal_attenuation);
    }
    else if(displacementsInvalid)
    {
        valuesInvalid = displacementsInvalid = false;
        floes_vtk.UnsafeUpdateDisplacements(floes.nodes.get(), floes.elems.get(), prms.temporal_attenuation);
    }
    else if(valuesInvalid)
    {
        valuesInvalid = false;
        floes_vtk.UnsafeUpdateValues(floes.nodes.get(), floes.elems.get(), prms.temporal_attenuation);
    }
    mutex.unlock();
    floes_vtk.UnsafeUpdateWaterLine(simulationTime, prms);
}

void icy::Model::RestoreFromSerializationBuffers(SimParams &prms)
{
    mutex.lock();
    floes.EvaluateStresses(prms, *floes.elems);
    floes.DistributeStresses();
    floes.EvaluateAllNormalTractions(prms);
    mutex.unlock();
    topologyInvalid = displacementsInvalid = valuesInvalid = true;
    if(!updateRequested) { updateRequested = true; emit requestGeometryUpdate(); }
}

long icy::Model::IdentifyDisconnectedRegions()
{
    auto t1 = std::chrono::high_resolution_clock::now();
    floes.IdentifyDisconnectedRegions();
    topologyInvalid = displacementsInvalid = valuesInvalid = true;
    if(!updateRequested) { updateRequested = true; emit requestGeometryUpdate(); }
    //mutex.lock();
    //floes.CreateEdges2();
    //mutex.unlock();
    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
}
/*
#include "modelcontroller.h"


void icy::ModelController::Reset()
{
    qDebug() << "icy::ModelController::Reset()";
    prms.Serialize();
    serializer.SaveParams(prms.serialization_buffer, SimParams::buffer_size);
    serializer.CloseFile();
    model.Reset();
    stepStats.clear();
    current_step = 0;
    ts.Reset();
    prms.Reset();
}

void icy::ModelController::SaveAs(QString fileName)
{
    // Current step becomes the initial setup(!)
    prms.Serialize();
    serializer.SaveParams(prms.serialization_buffer, SimParams::buffer_size);   // save to previous file
    serializer.CloseFile();
    serializer.CreateFile(fileName.toLocal8Bit().data(), SimParams::buffer_size);
    serializer.SaveParams(prms.serialization_buffer, SimParams::buffer_size);   // save to the new file
//    model.floes.WriteToSerializationBuffers();
//    serializer.Write(model.floes.node_buffer, model.floes.elems_buffer, 0, 0);
    model.floes.WriteToHD5(0,0, serializer.ds_nodes_handle, serializer.ds_elems_handle);

    stepStats.clear();
    ts.Reset();
    ts.nElems = model.floes.getElemCount();
    ts.nNodes = model.floes.getNodeCount();
    stepStats.push_back(ts);

    serializer.WriteSteps(1, &stepStats.front());
    current_step = 0;
}

void icy::ModelController::Load(QString fileName)
{
    qDebug() << "icy::ModelController::Load " << fileName;
    serializer.SaveParams(prms.serialization_buffer, SimParams::buffer_size);
    serializer.CloseFile();
    serializer.OpenFile(fileName.toLocal8Bit().data());
    serializer.LoadParams(prms.serialization_buffer, SimParams::buffer_size);
    prms.Deserialize();
    serializer.ReadSteps(stepStats);
    stepStats.front().BenchmarkingClear();
    GoToStep(0);
}

void icy::ModelController::ImportFloePatch(QString fileName)
{
    model.floes.ImportFloePatch(fileName, prms.CharacteristicLengthMax);
    current_step = 0;
    stepStats.clear();
    ts.Reset();
    ts.nElems = model.floes.getElemCount();
    ts.nNodes = model.floes.getNodeCount();
    stepStats.push_back(ts);
    _Write();
}

void icy::ModelController::Remesh()
{
    GoToStep(0);
    model.floes.Remesh(prms.CharacteristicLengthMax);
    current_step = 0;
    stepStats.clear();
    ts.Reset();
    ts.nElems = model.floes.getElemCount();
    ts.nNodes = model.floes.getNodeCount();
    stepStats.push_back(ts);
    _Write();
}


void icy::ModelController::GoToStep(int step)
{
    // this should not be invoked while the simulation is running in a worker thread
    qDebug() << "GoToStep " << step;
    if(serializer.fileIsOpen)
    {
        if(step < (int)stepStats.size()) {
            // load from file
            // qDebug() << "loading from file step " << step;
            current_step = step;
            ts = stepStats[step];
            model.floes.RestoreFromHD5(ts.nodeOffset, ts.elemOffset, ts.nNodes, ts.nElems,
                                       serializer.ds_nodes_handle, serializer.ds_elems_handle);
//            serializer.Read(model.floes.node_buffer, model.floes.elems_buffer,
//                            ts.nodeOffset, ts.elemOffset, ts.nNodes, ts.nElems);
            model.RestoreFromSerializationBuffers(prms);
        }
    }
}

void icy::ModelController::Trim()
{
    qDebug() << "Trim to " << current_step+1;
    stepStats.resize(current_step+1);
    icy::FrameInfo &fi = stepStats[current_step];
    unsigned long nodes_extent = fi.nodeOffset+fi.nNodes;
    unsigned long elems_extent = fi.elemOffset+fi.nElems;
    serializer.Trim(stepStats.size(), nodes_extent, elems_extent);
}

void icy::ModelController::Prepare()
{
    // add frame zero to the list
    if(stepStats.size() == 0) {
        ts.Reset();
        stepStats.push_back(ts);
        current_step = 0;
    } else {
        // stepStats.size() > 0
        // if current_step is less than stepStats.size(), then trim stepStats and data file
        if(current_step+1 < stepStats.size()) Trim();
        else if(current_step >= stepStats.size()) throw std::runtime_error("current step is not in the stats");
    }
}

void icy::ModelController::_BeginStep()
{
    model.floes.AssignLsIds();
    ts.nActiveNodes = model.floes.getFreeNodeCount();

    abort_requested = false;
    ts.BenchmarkingClear();
    ts.count_iterations = 0;
    ts.count_attempts = 0;
    ts.count_solves = 0;
    ts.solverProgress = 0;
    ts.solution_reached = false;
    ts.StepNumber = current_step+1;     // number of the tentative step (current plus one)

    icy::FrameInfo &cs = stepStats[current_step];
    ts.elemOffset = cs.elemOffset + cs.nElems;
    ts.nodeOffset = cs.nodeOffset + cs.nNodes;
    ts.TimeScaleFactor = cs.TimeScaleFactor;
}


void icy::ModelController::Step()
{
    requestToStop = false;
    auto t1 = std::chrono::high_resolution_clock::now();

    _BeginStep();
    do {
        ts.count_iterations = 0;
        icy::FrameInfo &cs = stepStats[current_step];
        // calculate time step
        if(ts.TimeScaleFactor > 0) ts.TimeScaleFactor--;
        ts.TimeStep = prms.InitialTimeStep / pow(2.0, (double)ts.TimeScaleFactor/4.0);
        ts.SimulationTime = cs.SimulationTime+ts.TimeStep;
        //qDebug() << "att: " << ts.count_attempts << "; ts: " << ts.TimeStep<<"; ts.TimeScaleFactor" << ts.TimeScaleFactor;

        model.InitialGuessTentativeVals(ts.TimeStep, prms.NewmarkBeta, prms.NewmarkGamma);

        bool converged = false;
        bool diverges = false;
        do {
            emit progressUpdated();

            double resultSqNorm;
            model.AssembleAndSolve(ts.b_clear_ls, ts.b_force_elem, ts.b_create_structure,
                                   ts.b_assemble, ts.b_solve, ts.b_pull_from_ls, prms, ts.TimeStep,
                                   ts.SimulationTime, resultSqNorm);
            if(abort_requested) {Aborting(); return;}
            ts.count_solves++;
            ts.count_iterations++;
            //qDebug() << "iter: "<< ts.count_iterations<<"; resultSqNorm: " << resultSqNorm;

            if(ts.count_iterations == 1) ts.Error0 = resultSqNorm;
            else if(ts.count_iterations >= prms.IterationsMin)
            {
                if(resultSqNorm < prms.ConvergenceCutoff) { converged = true; break; }
                // evaluate "converged" and "diverges"
                double ratio = resultSqNorm/ts.Error0;
                converged = ratio < prms.ConvergenceEpsilon;
                diverges = ratio > 1.0;
                //qDebug() << "ratio: " << ratio << "; diverges: " << diverges;
            }

        }while(!diverges && !converged && ts.count_iterations<prms.IterationsMax);
        ts.count_attempts++;

        if(converged || ts.TimeScaleFactor >=16) ts.solution_reached=true;
        else ts.TimeScaleFactor += 4;

    } while(!ts.solution_reached);

    model.AcceptTentativeValues(prms);
    Fracture();
    if(abort_requested) {Aborting(); return;}

    ts.nElems = model.floes.getElemCount();
    ts.nNodes = model.floes.getNodeCount();
    stepStats.push_back(ts);
    current_step = ts.StepNumber;
    auto t2 = std::chrono::high_resolution_clock::now();
    stepStats.back().b_total += std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
    _Write();
    emit stepCompleted();
}

void icy::ModelController::_Write()
{
    if(!serializer.fileIsOpen) return;
    serializer.WriteSteps(stepStats.size(), &stepStats.back());
    model.floes.WriteToHD5(ts.nodeOffset, ts.elemOffset, serializer.ds_nodes_handle, serializer.ds_elems_handle);


//    model.floes.WriteToSerializationBuffers();
//    serializer.WriteAll(model.floes.node_buffer, model.floes.elems_buffer,
//                        ts.nodeOffset, ts.elemOffset, stepStats.size(), &stepStats.back());
}

void icy::ModelController::RequestAbort()
{
    abort_requested = true;
    // if needed, abort the solver
}

void icy::ModelController::Aborting()
{
    //perform any cleanup if step was aborted
    qDebug() << "icy::ModelController::Aborting()";
    abort_requested = false;
    ts.solverProgress = 0;
    emit stepAborted();
}

void icy::ModelController::Fracture()
{
    if(!prms.fracture_enable)
    {
        model.floes.EvaluateAllNormalTractions(prms);
        return;
    }

    model.mutex.lock();
    ts.b_compute_fracture_directions += model.floes.ComputeFractureDirections(prms, ts.TimeStep, true);
    model.mutex.unlock();
    int count=0;
    model.floes_vtk.update_minmax = false;

    while(model.floes.maxNode != nullptr && count < prms.fracture_max_substeps && !abort_requested)
    {
        model.FractureStep(prms, ts.TimeStep, ts.SimulationTime, ts.b_local_substep,
                           ts.b_compute_fracture_directions, ts.b_split, ts.b_infer_support);
        count++;
        emit fractureUpdated();
    }
    model.floes_vtk.update_minmax = true;

    if(count>0) ts.b_identify_regions += model.IdentifyDisconnectedRegions();
}

*/