#include <vtkPointData.h>
#include <QtGlobal>
#include "model.h"

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
        Eigen::Matrix<double,DOFS,1> xt_xn = nd->ut - nd->un;
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
        Eigen::Matrix<double,DOFS,1> xt_xn = nd->ut - nd->un;
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
        nd->normal_n = Eigen::Vector3d::Zero(); // allow to accumulate
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
    floes.RecomputeElasticityMatrix(prms);
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
