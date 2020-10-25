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
        for(int k=0;k<3;k++) nd->str_b[k] = nd->str_m[k] = nd->str_b_top[k] = nd->str_b_bottom[k] = 0;
        nd->str_s[0] = nd->str_s[1] = 0;
    }

#pragma omp parallel for
    for(std::size_t i=0;i<nElems;i++) (*elems)[i]->DistributeStresses();
}

long icy::Geometry::ComputeFractureDirections(SimParams &prms, double timeStep, bool startingFracture)
{
    auto t1 = std::chrono::high_resolution_clock::now();

    maxNode=nullptr;
    breakable_range.clear();
    double temporal_attenuation = prms.temporal_attenuation;
    EvaluateStresses(prms);
    DistributeStresses();

    float threashold = prms.normal_traction_threshold;

    std::size_t nNodes = nodes->size();

    // first, try to find a "breakable" crack tip
    // no need for multi-core, because there are few existing crack tips
    if(!startingFracture)
    {
        for(std::size_t i=0;i<nNodes;i++)
        {
            icy::Node *nd = (*nodes)[i];
            if(!nd->crack_tip) continue;
            nd->InitializeFan();
            nd->ComputeFanVariablesAlt(prms);
            if(nd->crack_tip && nd->max_normal_traction > threashold)
            {
                maxNode = nd;
                breakable_range.push_back(maxNode);
                break;
            }
        }
    }

    // crack tip not found - do a more elaborate search
    if(maxNode == nullptr)
    {
        // compute max traction and potential fracture direction
#pragma omp parallel for
        for(std::size_t i=0;i<nNodes;i++)
        {
            icy::Node *nd = (*nodes)[i];
            nd->InitializeFan();
            nd->ComputeFanVariablesAlt(prms);

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
            if(!nd->core_node) nd->dir=Eigen::Vector2f::Zero();
        }

        // put "core_nodes" into breakable_range

        std::copy_if(nodes->begin(), nodes->end(), std::back_inserter(breakable_range),
                     [](icy::Node *nd){return nd->core_node;});

        if(breakable_range.size()==0)
        {
            maxNode=nullptr;
            auto t2 = std::chrono::high_resolution_clock::now();
            return std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
        }

        // give priority to an existing crack tip
        std::vector<Node*>::iterator it_nd = std::find_if(breakable_range.begin(),
                                                          breakable_range.end(),
                                                          [](Node *nd) {return nd->crack_tip;});

        // if no tips, break the node under strongest load
        if(it_nd == breakable_range.end())
        {
            it_nd = std::max_element(breakable_range.begin(), breakable_range.end(),
                                     [](Node *nd1, Node *nd2) {
                    return nd1->max_normal_traction < nd2->max_normal_traction; });
        }

        maxNode = *it_nd;
    }
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

    // create a list of elements corresponding to breakable_range
    local_elems_set.clear();
    local_elems.clear();
    for(icy::Node *nd : breakable_range) for(icy::Element *elem : nd->adjacent_elems) local_elems_set.insert(elem);
    if(local_elems_set.size()>0) std::copy(local_elems_set.begin(), local_elems_set.end(), std::back_inserter(local_elems));

    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
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


long icy::Geometry::SplitNode(SimParams &prms)
{
    auto t1 = std::chrono::high_resolution_clock::now();

    if(maxNode == nullptr) {
        // qDebug() << "SplitNode: nothing to split";
        return 0;
    }

    icy::Node* nd = maxNode;
    nd->core_node = false;

    // qDebug() << "max nd: " << nd->locId << "; traction: " << nd->max_normal_traction;

    // the central node nd splits into itself an mainSplit
    icy::Node *mainSplit = nullptr;

    // subsequent calculations are based on the fracture direction where the traction is maximal
    icy::Node::SepStressResult &ssr = nd->result_with_max_traction;

    // make sure that the interior node has two split faces
    bool isBoundary = (ssr.faces[1] == nullptr);
    if(isBoundary != nd->isBoundary) std::runtime_error("isBoundary != nd->isBoundary");

    // iterate over the "fan" and replace nd with mainSplit on the elements between angle_fwd and angle_bwd
    if(isBoundary)
    {
        for(icy::Node::Sector &f : nd->fan)
            if(f.angle0 > ssr.angle_fwd && f.angle1 > ssr.angle_fwd) {
                if(mainSplit == nullptr) mainSplit = AddNode(nd);
                f.face->ReplaceNode(nd, mainSplit);
            }
    }
    else if(ssr.angle_fwd > ssr.angle_bwd)
    {
        for(icy::Node::Sector &f : nd->fan)
            if((f.angle0 > ssr.angle_fwd || f.angle0 < ssr.angle_bwd)
                    && (f.angle1 > ssr.angle_fwd || f.angle1 < ssr.angle_bwd)) {
                if(mainSplit == nullptr) mainSplit = AddNode(nd);
                f.face->ReplaceNode(nd, mainSplit);
            }
    }
    else
    {
        // (ssr.angle_fwd < ssr.angle_bwd)
        for(icy::Node::Sector &f : nd->fan)
            if(f.angle0 > ssr.angle_fwd && f.angle1 > ssr.angle_fwd &&
                    f.angle0 < ssr.angle_bwd && f.angle1 < ssr.angle_bwd) {
                if(mainSplit == nullptr) mainSplit = AddNode(nd);
                f.face->ReplaceNode(nd, mainSplit);
            }
    }

    // split face #0
    icy::Node *nd0 = ssr.e[0].getOtherNode(nd);
    icy::Node *nd1 = ssr.e[1].getOtherNode(nd);
    int nd0idx = nd0->locId;
    int nd1idx = nd1->locId;
    bool forwardEdge = (nd0idx < nd1idx);

    double factor0 = sin(ssr.phi[0])*ssr.e[0].getVec(nd).norm();
    double factor1 = sin(ssr.theta[0])*ssr.e[1].getVec(nd).norm();
    double whereToSplit = forwardEdge ? factor1/(factor0+factor1) : factor0/(factor0+factor1);

    icy::Edge splitEdge0 = getEdgeByNodalIdx(nd0idx, nd1idx);

    /*
    Eigen::Vector3d dir = ssr.tn;
    if(nd->crack_tip)
    {
        dir+=2*nd->weakening_direction;
        dir.normalize();
    }
    */

    SplitEdge(splitEdge0, whereToSplit, nd, mainSplit, forwardEdge, prms, ssr.tn);

    if(ssr.faces[1] != nullptr)
    {
        // center node is interior, split face #0
        icy::Node *nd0 = ssr.e[2].getOtherNode(nd);
        icy::Node *nd1 = ssr.e[3].getOtherNode(nd);
        int nd0idx = nd0->locId;
        int nd1idx = nd1->locId;
        bool forwardEdge = (nd0idx < nd1idx);
        double factor0 = sin(ssr.phi[1])*ssr.e[2].getVec(nd).norm();
        double factor1 = sin(ssr.theta[1])*ssr.e[3].getVec(nd).norm();
        double whereToSplit = forwardEdge ? factor1/(factor0+factor1) : factor0/(factor0+factor1);

        icy::Edge splitEdge0 = getEdgeByNodalIdx(nd0idx, nd1idx);
        SplitEdge(splitEdge0, whereToSplit, nd, mainSplit, !forwardEdge, prms, -ssr.tn);
    }

    nd->weakening_direction = Eigen::Vector2f::Zero();
    nd->crack_tip = false;

    CreateEdges2();
    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
}


void icy::Geometry::SplitEdge(icy::Edge edge, double where,
                              icy::Node *centerNode, icy::Node* &splitNode,
                              bool forwardDirection,
                              SimParams &prms, Eigen::Vector2f dir)
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

    if(splitNode == nullptr) splitNode = AddNode(centerNode);

    if(edge.isBoundary)
    {
        // break edge with 2 additional nodes; add one element on the "nd" side

        icy::Node *split0 = AddNode();
        icy::Node *split1 = AddNode();

        split0->InitializeFromAdjacent(edge.nds[0], edge.nds[1], where);
        split1->InitializeFromAnother(split0);

        icy::Element *existingFace = edge.getTheOnlyElement();
        bool orientation = existingFace->Orientation(edge.nds[0], edge.nds[1]);

        icy::Element *additionalFace = AddElement();

        if(forwardDirection)
        {
            existingFace->Initialize(edge.nds[0], split0, centerNode, orientation);
            additionalFace->Initialize(split1, edge.nds[1], splitNode, orientation);
        } else
        {
            existingFace->Initialize(edge.nds[0], split0, splitNode, orientation);
            additionalFace->Initialize(split1, edge.nds[1], centerNode, orientation);
        }
        existingFace->PrecomputeStiffnessMatrix(prms, elasticityMatrix, D_mats);
        additionalFace->PrecomputeStiffnessMatrix(prms, elasticityMatrix, D_mats);
    }
    else
    {
        icy::Element *existingFace0 = edge.getElementWithNode(centerNode);
        bool orientation0 = existingFace0->Orientation(edge.nds[0], edge.nds[1]);

        // insert one point to the edge; add element to each side of the edge (2 total)
        // CRACK TIP NODE
        icy::Node *split0 = AddNode();
        split0->InitializeFromAdjacent(edge.nds[0], edge.nds[1], where);
        split0->weakening_direction = dir;
        split0->crack_tip = true;

        icy::Node *farNode = edge.getFarNode(centerNode);
        icy::Element *existingFace1 = edge.getElementWithNode(farNode);
        bool orientation1 = existingFace1->Orientation(edge.nds[0], edge.nds[1]);

        icy::Element *additionalFace0 = AddElement();
        icy::Element *additionalFace1 = AddElement();

        existingFace1->Initialize(edge.nds[0], split0, farNode, orientation1);
        additionalFace1->Initialize(split0, edge.nds[1], farNode, orientation1);
        if(forwardDirection)
        {
            existingFace0->Initialize(edge.nds[0], split0, centerNode, orientation0);
            additionalFace0->Initialize(split0, edge.nds[1], splitNode, orientation0);
        }
        else
        {
            existingFace0->Initialize(edge.nds[0], split0, splitNode, orientation0);
            additionalFace0->Initialize(split0, edge.nds[1], centerNode, orientation0);
        }
        existingFace0->PrecomputeStiffnessMatrix(prms, elasticityMatrix, D_mats);
        additionalFace0->PrecomputeStiffnessMatrix(prms, elasticityMatrix, D_mats);
        existingFace1->PrecomputeStiffnessMatrix(prms, elasticityMatrix, D_mats);
        additionalFace1->PrecomputeStiffnessMatrix(prms, elasticityMatrix, D_mats);
    }
}

void icy::Geometry::SplitAlongExistingEdge(Edge edge, Node *centerNode, Node* &splitNode,
                                           int oppositeNodeIdx, bool forwardDirection, Eigen::Vector2f dir)
{
    // qDebug() << "SplitAlongExistingEdge; oppositeNodeIdx: " << oppositeNodeIdx;

    icy::Element *existingFace = edge.getElementWithNode(centerNode);
    bool orientation = existingFace->Orientation(edge.nds[0], edge.nds[1]);

    icy::Node *oppositeNode = edge.nds[oppositeNodeIdx];

    icy::Edge edgeBeingSplit = getEdgeByNodalIdx(oppositeNode->locId, centerNode->locId);
    if(edgeBeingSplit.isBoundary) return; // trying to split a boundary edge

    if(splitNode == nullptr) splitNode = AddNode(centerNode);

    bool boundary_detected = oppositeNode->isBoundary;
    //    qDebug() << "edge: " << edge->nds[0]->locId << " -- " << edge->nds[1]->locId;
    //    qDebug() << "elem0: " << edge->elems[0] << "; elem1: " << edge->elems[1];
    //    qDebug() << "boundary: " << boundary_detected << "; fwd: " << forwardDirection;

    auto fan_iter = std::find_if(oppositeNode->fan.begin(),
                                 oppositeNode->fan.end(),
                                 [existingFace](icy::Node::Sector &f)
    {return f.face == existingFace;});
    if(fan_iter == oppositeNode->fan.end()) throw std::runtime_error("existing face not in the fan");
    icy::Node::Sector &fanItem = *fan_iter;
    int whichEdge;
    if(fanItem.nd[0] == centerNode) whichEdge = 0;
    else if(fanItem.nd[1] == centerNode) whichEdge = 1;
    else throw std::runtime_error("split edge not detected");

    // reconnect the element
    if((oppositeNodeIdx == 0 && forwardDirection) || (oppositeNodeIdx == 1 && !forwardDirection)) {
        existingFace->Initialize(edge.nds[0], edge.nds[1], splitNode, orientation);
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


long icy::Geometry::IdentifyDisconnectedRegions()
{
    auto t1 = std::chrono::high_resolution_clock::now();
    regions.clear();
    for(icy::Element *e : *elems) e->traversed = false;  // set to not-traversed

    unsigned current_region = 0;
    wave.clear();
    wave.reserve(elems->size());
    area = 0;
    for(icy::Element *e : *elems)
    {
        if(e->traversed) continue;

        wave.push_back(e);
        unsigned count_elems = 0;
        double region_area = 0;
        while(wave.size() > 0)
        {
            icy::Element *elem = wave.back();
            wave.pop_back();
            count_elems++;
            region_area += elem->area_initial;
            elem->traversed = true;
            elem->region = current_region;
            for(int i=0;i<3;i++)
            {
                icy::Element *adj_e = elem->adj_elems[i];
                if(adj_e!= nullptr && !adj_e->traversed) wave.push_back(adj_e);
            }
        }
        regions.push_back(std::make_tuple(current_region, region_area, count_elems));
        current_region++;
        area+=region_area;
    }

    // for testing
//    std::cout << "printing regions:\n";
//    for(std::tuple<unsigned, double, unsigned> &r : regions)
//        std::cout << std::get<0>(r) << ": " << std::get<1>(r) << "; " << std::get<2>(r) <<  std::endl;
//    std::cout << "============= \n";

    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
}
/*
long icy::Geometry::RemoveDegenerateFragments()
{
    auto t1 = std::chrono::high_resolution_clock::now();
    double avg_elem_area = area/elems->size();
    bool region_was_removed;
    do
    {
        region_was_removed = false;
        auto iter = std::min_element(regions.begin(), regions.end(),
                                      [](std::tuple<unsigned, double, unsigned> r1,
                                      std::tuple<unsigned, double, unsigned> r2)
        {return std::get<2>(r1) < std::get<2>(r2);});
        if(std::get<2>(*iter) <= 2)
        {
            unsigned idx = std::get<0>(*iter);
            RemoveRegion(idx);
            regions.erase(iter);
            region_was_removed = true;
        }

        iter = std::min_element(regions.begin(), regions.end(),
                                [](std::tuple<unsigned, double, unsigned> r1,
                                std::tuple<unsigned, double, unsigned> r2)
        {return std::get<1>(r1) < std::get<1>(r2);});

        if(std::get<1>(*iter) < avg_elem_area*1.5)
        {
            unsigned idx = std::get<0>(*iter);
            RemoveRegion(idx);
            regions.erase(iter);
            region_was_removed = true;
        }

    } while(region_was_removed && regions.size() > 0);

    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
}

void icy::Geometry::RemoveRegion(unsigned idx)
{
    std::unordered_set<icy::Node*>nds;
    for(icy::Element *elem : *elems)
    {
        if(elem->region == idx) {
            for(int k=0;k<3;k++) nds.insert(elem->nds[k]);
            pool_elems.free(elem);
        }
    }

    elems->erase(std::remove_if(elems->begin(), elems->end(),
                                [idx](icy::Element *elem){return elem->region==idx;}), elems->end());

    for(icy::Node *nd : nds) pool_nodes.destroy(nd);

    nodes->erase(std::remove_if(nodes->begin(), nodes->end(),
                                [nds](icy::Node *nd){return nds.find(nd)!=nds.end();}), nodes->end());

    for(std::size_t i=0;i<nodes->size();i++) (*nodes)[i]->locId=i;

    CreateEdges2();
}
*/
