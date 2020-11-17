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
    model.floes.RecomputeElasticityMatrix(prms);

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
