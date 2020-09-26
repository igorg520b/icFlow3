#include "geometry.h"
#include "model.h"
#include <numeric>
#include <algorithm>
#include <iterator>

void icy::Geometry::EvaluateStresses(SimParams &prms)
{
#pragma omp parallel for
    for(std::size_t i=0;i<elems->size();i++) (*elems)[i]->EvaluateStresses(prms, elasticityMatrix, D_mats);
}

void icy::Geometry::DistributeStresses()
{
    std::size_t nElems = elems->size();
    std::size_t nNodes = nodes->size();
#pragma omp parallel for
    for(std::size_t i=0;i<nNodes;i++) {
        icy::Node *nd = (*nodes)[i];
        nd->str_b = nd->str_m = nd->str_b_top = nd->str_b_bottom = Eigen::Vector3d::Zero();
        nd->str_s = nd->str_b_top_principal = nd->str_b_bottom_principal = Eigen::Vector2d::Zero();
    }
    // this step is performed sequentially
    for(std::size_t i=0;i<nElems;i++) (*elems)[i]->DistributeStresses();
}


void icy::Geometry::ComputeFractureDirections(SimParams &prms, double timeStep, bool startingFracture)
{
    double temporal_attenuation = prms.temporal_attenuation;
    EvaluateStresses(prms);
    DistributeStresses();

    double threashold = prms.normal_traction_threshold;

    std::size_t nNodes = nodes->size();
    // compute max traction and potential fracture direction
#pragma omp parallel for
    for(std::size_t i=0;i<nNodes;i++)
    {
        icy::Node *nd = (*nodes)[i];
        nd->ComputeFanVariables(prms);

        if(startingFracture)
        {
            nd->core_node = (nd->max_normal_traction > threashold && nd->timeLoadedAboveThreshold >= temporal_attenuation);
            if(nd->max_normal_traction > threashold) nd->timeLoadedAboveThreshold+=timeStep;
            else nd->timeLoadedAboveThreshold=0;
        }
        else
        {
            if(nd->max_normal_traction < threashold || nd->timeLoadedAboveThreshold < temporal_attenuation) nd->core_node = false;
        }

        if(nd->crack_tip && nd->max_normal_traction > threashold)
            nd->core_node = true;
        if(!nd->core_node) nd->dir=Eigen::Vector3d::Zero();
    }

    // put "core_nodes" into breakable_range
    breakable_range.clear();
    std::copy_if(nodes->begin(), nodes->end(), std::back_inserter(breakable_range),
                 [](icy::Node *nd){return nd->core_node;});

    if(breakable_range.size()==0)
    {
        maxNode=nullptr;
        return;
    }

    auto it_nd = std::max_element(breakable_range.begin(), breakable_range.end(),
                                  [](Node *nd1, Node *nd2) {
            return nd1->max_normal_traction < nd2->max_normal_traction; });

    maxNode = *it_nd;
    maxNode->dir*=3;

    local_support.clear();
    local_support.push_back(maxNode);

    // prevent the nodes surrounging maxNode from fracturing
    CreateSupportRange(prms.substep_radius, local_support);
    for(Node *nd : local_support) nd->timeLoadedAboveThreshold = 0;

    CreateSupportRange(prms.substep_radius, breakable_range);

    // for visualization - mark support range (stored in breakable_range)
    for(icy::Node *nd : *nodes) nd->support_node = false;
    for(icy::Node *nd : breakable_range) nd->support_node = true; // for visualization
}

void icy::Geometry::CreateSupportRange(int neighborLevel, std::vector<Node*> &initial_set)
{
    tmp_range0->clear();
    tmp_range1->clear();

    std::copy(initial_set.begin(),
              initial_set.end(),
              std::inserter(*tmp_range0, tmp_range0->end()));

    for(int i=0;i<neighborLevel;i++)
    {
        std::swap(tmp_range0, tmp_range1);
        for(icy::Node *nd : *tmp_range1) {
            tmp_range0->insert(nd);
            for(icy::Node *nd2 : nd->adjacent_nodes)
                tmp_range0->insert(nd2);
        }
    }
    initial_set.clear();
    std::copy(tmp_range0->begin(), tmp_range0->end(), std::back_inserter(initial_set));
}



void icy::Geometry::SplitNode(SimParams &prms)
{
    if(maxNode == nullptr) {
        // qDebug() << "SplitNode: nothing to split";
        return;
    }

    icy::Node* nd = maxNode;
    nd->core_node = false;

    // qDebug() << "max nd: " << nd->locId << "; traction: " << nd->max_normal_traction;

    // the central node nd splits into itself an mainSplit
    icy::Node *mainSplit = nullptr;

    // subsequent calculations are based on the fracture direction where the traction is maximal
    icy::Node::SepStressResult &ssr = nd->sep_stress_results[nd->idxSepStressResult];

    // make sure that the interior node has two split faces
    bool isBoundary = (ssr.faces[1] == nullptr);
    if(isBoundary != nd->isBoundary) std::runtime_error("isBoundary != nd->isBoundary");

    // iterate over the "fan" and replace nd with mainSplit on the elements between angle_fwd and angle_bwd
    if(isBoundary)
    {
        for(icy::Node::FanPrecomp &f : nd->fan)
            if(f.angle0 > ssr.angle_fwd && f.angle1 > ssr.angle_fwd) {
                if(mainSplit == nullptr) mainSplit = AddNode(nd);
                f.face->ReplaceNode(nd, mainSplit);
            }
    }
    else if(ssr.angle_fwd > ssr.angle_bwd)
    {
        for(icy::Node::FanPrecomp &f : nd->fan)
            if((f.angle0 > ssr.angle_fwd || f.angle0 < ssr.angle_bwd)
                    && (f.angle1 > ssr.angle_fwd || f.angle1 < ssr.angle_bwd)) {
                if(mainSplit == nullptr) mainSplit = AddNode(nd);
                f.face->ReplaceNode(nd, mainSplit);
            }
    }
    else
    {
        // (ssr.angle_fwd < ssr.angle_bwd)
        for(icy::Node::FanPrecomp &f : nd->fan)
            if(f.angle0 > ssr.angle_fwd && f.angle1 > ssr.angle_fwd &&
                    f.angle0 < ssr.angle_bwd && f.angle1 < ssr.angle_bwd) {
                if(mainSplit == nullptr) mainSplit = AddNode(nd);
                f.face->ReplaceNode(nd, mainSplit);
            }
    }

    // split face #0
    icy::Node *nd0 = ssr.e[0]->getOtherNode(nd);
    icy::Node *nd1 = ssr.e[1]->getOtherNode(nd);
    int nd0idx = nd0->locId;
    int nd1idx = nd1->locId;
    bool forwardEdge = (nd0idx < nd1idx);

    double factor0 = sin(ssr.phi[0])*ssr.e[0]->getVec(nd).norm();
    double factor1 = sin(ssr.theta[0])*ssr.e[1]->getVec(nd).norm();
    double whereToSplit = forwardEdge ? factor1/(factor0+factor1) : factor0/(factor0+factor1);

    icy::Edge *splitEdge0 = getEdgeByNodalIdx(nd0idx, nd1idx);
    Eigen::Vector3d dir = ssr.tn;
    if(nd->crack_tip)
    {
        dir+=2*nd->weakening_direction;
        dir.normalize();
    }
    SplitEdge(splitEdge0, whereToSplit, nd, mainSplit, forwardEdge, prms, ssr.tn);

    if(ssr.faces[1] != nullptr)
    {
        // center node is interior, split face #0
        icy::Node *nd0 = ssr.e[2]->getOtherNode(nd);
        icy::Node *nd1 = ssr.e[3]->getOtherNode(nd);
        int nd0idx = nd0->locId;
        int nd1idx = nd1->locId;
        bool forwardEdge = (nd0idx < nd1idx);
        double factor0 = sin(ssr.phi[1])*ssr.e[2]->getVec(nd).norm();
        double factor1 = sin(ssr.theta[1])*ssr.e[3]->getVec(nd).norm();
        double whereToSplit = forwardEdge ? factor1/(factor0+factor1) : factor0/(factor0+factor1);

        icy::Edge *splitEdge0 = getEdgeByNodalIdx(nd0idx, nd1idx);
        SplitEdge(splitEdge0, whereToSplit, nd, mainSplit, !forwardEdge, prms, -ssr.tn);
    }

    nd->weakening_direction = Eigen::Vector3d::Zero();
    nd->crack_tip = false;

    CreateEdges();
}




void icy::Geometry::SplitEdge(icy::Edge *edge, double where,
                              icy::Node *centerNode, icy::Node* &splitNode,
                              bool forwardDirection,
                              SimParams &prms, Eigen::Vector3d dir)
{
    if(where < prms.fracture_epsilon) {
        SplitAlongExistingEdge(edge, centerNode, splitNode, 1, forwardDirection, dir);
        return;
    }
    else if(where > 1.0-prms.fracture_epsilon)
    {
        SplitAlongExistingEdge(edge, centerNode, splitNode, 0, forwardDirection, dir);
        return;
    }

    // qDebug() << "SplitEdge";
    if(splitNode == nullptr) splitNode = AddNode(centerNode);

    if(edge->isBoundary)
    {
        // break edge with 2 additional nodes; add one element on the "nd" side

        icy::Node *split0 = AddNode();
        icy::Node *split1 = AddNode();

        split0->InitializeFromAdjacent(edge->nds[0], edge->nds[1], where);
        split1->InitializeFromAnother(split0);

        icy::Element *existingFace = edge->getTheOnlyElement();
        bool orientation = existingFace->Orientation(edge->nds[0], edge->nds[1]);

        icy::Element *additionalFace = AddElement();

        if(forwardDirection)
        {
            existingFace->Initialize(edge->nds[0], split0, centerNode, orientation);
            additionalFace->Initialize(split1, edge->nds[1], splitNode, orientation);
        } else
        {
            existingFace->Initialize(edge->nds[0], split0, splitNode, orientation);
            additionalFace->Initialize(split1, edge->nds[1], centerNode, orientation);
        }
        existingFace->PrecomputeStiffnessMatrix(prms, elasticityMatrix, D_mats);
        additionalFace->PrecomputeStiffnessMatrix(prms, elasticityMatrix, D_mats);
    }
    else
    {
        icy::Element *existingFace0 = edge->getElementWithNode(centerNode);
        bool orientation0 = existingFace0->Orientation(edge->nds[0], edge->nds[1]);

        // insert one point to the edge; add element to each side of the edge (2 total)
        // CRACK TIP NODE
        icy::Node *split0 = AddNode();
        split0->InitializeFromAdjacent(edge->nds[0], edge->nds[1], where);
//        split0->weakening_direction = (split0->xt - centerNode->xt).block(0,0,3,1).normalized();
        split0->weakening_direction = dir;
        split0->crack_tip = true;

        icy::Node *farNode = edge->getFarNode(centerNode);
        icy::Element *existingFace1 = edge->getElementWithNode(farNode);
        bool orientation1 = existingFace1->Orientation(edge->nds[0], edge->nds[1]);

        icy::Element *additionalFace0 = AddElement();
        icy::Element *additionalFace1 = AddElement();

        existingFace1->Initialize(edge->nds[0], split0, farNode, orientation1);
        additionalFace1->Initialize(split0, edge->nds[1], farNode, orientation1);
        if(forwardDirection)
        {
            existingFace0->Initialize(edge->nds[0], split0, centerNode, orientation0);
            additionalFace0->Initialize(split0, edge->nds[1], splitNode, orientation0);
        }
        else
        {
            existingFace0->Initialize(edge->nds[0], split0, splitNode, orientation0);
            additionalFace0->Initialize(split0, edge->nds[1], centerNode, orientation0);
        }
        existingFace0->PrecomputeStiffnessMatrix(prms, elasticityMatrix, D_mats);
        additionalFace0->PrecomputeStiffnessMatrix(prms, elasticityMatrix, D_mats);
        existingFace1->PrecomputeStiffnessMatrix(prms, elasticityMatrix, D_mats);
        additionalFace1->PrecomputeStiffnessMatrix(prms, elasticityMatrix, D_mats);
    }
}

void icy::Geometry::SplitAlongExistingEdge(Edge *edge, Node *centerNode, Node* &splitNode,
                                           int oppositeNodeIdx, bool forwardDirection, Eigen::Vector3d dir)
{
    // qDebug() << "SplitAlongExistingEdge; oppositeNodeIdx: " << oppositeNodeIdx;

    icy::Element *existingFace = edge->getElementWithNode(centerNode);
    bool orientation = existingFace->Orientation(edge->nds[0], edge->nds[1]);

    icy::Node *oppositeNode = edge->nds[oppositeNodeIdx];

    icy::Edge *edgeBeingSplit = getEdgeByNodalIdx(oppositeNode->locId, centerNode->locId);
    if(edgeBeingSplit->isBoundary) return; // trying to split a boundary edge

    if(splitNode == nullptr) splitNode = AddNode(centerNode);

    bool boundary_detected = oppositeNode->isBoundary;
    //    qDebug() << "edge: " << edge->nds[0]->locId << " -- " << edge->nds[1]->locId;
    //    qDebug() << "elem0: " << edge->elems[0] << "; elem1: " << edge->elems[1];
    //    qDebug() << "boundary: " << boundary_detected << "; fwd: " << forwardDirection;

    auto fan_iter = std::find_if(oppositeNode->fan.begin(),
                                 oppositeNode->fan.end(),
                                 [existingFace](icy::Node::FanPrecomp &f)
    {return f.face == existingFace;});
    if(fan_iter == oppositeNode->fan.end()) throw std::runtime_error("existing face not in the fan");
    icy::Node::FanPrecomp &fanItem = *fan_iter;
    int whichEdge;
    if(fanItem.nd[0] == centerNode) whichEdge = 0;
    else if(fanItem.nd[1] == centerNode) whichEdge = 1;
    else throw std::runtime_error("split edge not detected");

    // reconnect the element
    if((oppositeNodeIdx == 0 && forwardDirection) || (oppositeNodeIdx == 1 && !forwardDirection)) {
        existingFace->Initialize(edge->nds[0], edge->nds[1], splitNode, orientation);
    }
    else
    {
        // nothing changes
        //existingFace->Initialize(edge->nds[0], edge->nds[1], centerNode, orientation);
    }

    if(!boundary_detected)
    {
        // CRACK TIP NODE
        oppositeNode->weakening_direction = dir;
//        oppositeNode->weakening_direction = (oppositeNode->xt - centerNode->xt).block(0,0,3,1).normalized();
        oppositeNode->crack_tip = true;
    }
    else
    {
        // break boundary
        icy::Node *split = AddNode(oppositeNode);

        //        qDebug() << "whichEdge " << whichEdge;
        if(whichEdge == 1) fan_iter++;
        if(fan_iter == oppositeNode->fan.end()) throw std::runtime_error("not supposed to happen");
        while(fan_iter != oppositeNode->fan.end())
        {
            fan_iter->face->ReplaceNode(oppositeNode, split);
            fan_iter++;
        }
    }

    //existingFace->PrecomputeStiffnessMatrix(prms, elasticityMatrix, D_mats);
}
