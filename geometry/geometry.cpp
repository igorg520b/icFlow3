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
        nd->potentially_can_fracture = false;
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
            if(nd->potentially_can_fracture)
            {
                nd->InitializeFan();
                nd->ComputeFanVariablesAlt(prms);
            }
            else
            {
                nd->max_normal_traction = 0;
            }

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


//========================== ALTERNATIVE ALGORITHM FOR SPLITTING
long icy::Geometry::SplitNodeAlt(SimParams &prms)
{
    auto t1 = std::chrono::high_resolution_clock::now();

    if(maxNode == nullptr) throw std::runtime_error("trying to split nullptr");

    affected_elements_during_split.clear();
    affected_nodes_during_split.clear();

    icy::Node* nd = maxNode;
    nd->core_node = false;
    affected_nodes_during_split.insert(nd);

    // subsequent calculations are based on the fracture direction where the traction is maximal
    icy::Node::SepStressResult &ssr = nd->result_with_max_traction;

    // make sure that the interior node has two split faces
    bool isBoundary = (ssr.faces[1] == nullptr);
    if(isBoundary != nd->isBoundary) std::runtime_error("isBoundary != nd->isBoundary");

    icy::Edge splitEdge_fw;
    EstablishSplittingEdge(splitEdge_fw, nd,
                                ssr.phi[0], ssr.theta[0], prms.fracture_epsilon,
                                ssr.e[0], ssr.e[1], ssr.e_opposite[0], ssr.faces[0],prms);

    if(isBoundary)
    {
        Fix_X_Topology(nd); // split as fan.front().e[0] --- splitEdge_fw --- fan.back().e[1]
    }
    else
    {
        // determine splitEdge_bw (create if necessary)
        icy::Edge splitEdge_bw;
        EstablishSplittingEdge(splitEdge_bw, nd,
                                    ssr.phi[1], ssr.theta[1], prms.fracture_epsilon,
                                    ssr.e[2], ssr.e[3], ssr.e_opposite[1], ssr.faces[1],prms);

        Fix_X_Topology(nd);
        // split between splitEdge_fw and splitEdge_bw
        Node *split1=splitEdge_bw.getOtherNode(nd);
        if(split1->isBoundary)
        {
            Fix_X_Topology(split1);
        }
        else
        {
            split1->crack_tip = true;
            split1->weakening_direction = Eigen::Vector2f(split1->xt.x()-nd->xt.x(), split1->xt.y()-nd->xt.y());
            split1->weakening_direction.normalize();
        }
    }
    Node *split0=splitEdge_fw.getOtherNode(nd);
    if(split0->isBoundary)
    {
        Fix_X_Topology(split0);
    }
    else
    {
        split0->crack_tip = true;
        split0->weakening_direction = Eigen::Vector2f(split0->xt.x()-nd->xt.x(), split0->xt.y()-nd->xt.y());
        split0->weakening_direction.normalize();
    }

    nd->weakening_direction = Eigen::Vector2f::Zero();
    nd->crack_tip = false;

    // TODO: replace with UpdateEdges
    CreateEdges2();

    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
}

void icy::Geometry::Fix_X_Topology(Node *nd)
{
    nd->PrepareFan2(false);
    Node *split=AddNode();
    split->InitializeFromAnother(nd);
    bool replacing = true;
    for(Node::Sector &s : nd->fan)
    {
        if(replacing)s.face->ReplaceNode(nd, split);
        if(s.e[1].toSplit) replacing=!replacing;
    }
}


void icy::Geometry::EstablishSplittingEdge(Edge &splitEdge, Node* nd,
                            const float phi, const float theta, const float fracture_epsilon,
                            const Edge e0, const Edge e1, const Edge e_opposite, Element *elem, SimParams &prms)
{
    icy::Node *nd0 = e0.getOtherNode(nd);
    icy::Node *nd1 = e1.getOtherNode(nd);

    Eigen::Vector2f nd_vec((float)nd->xt.x(),(float)nd->xt.y());
    Eigen::Vector2f nd0_vec((float)nd0->xt.x(),(float)nd0->xt.y());
    Eigen::Vector2f nd1_vec((float)nd1->xt.x(),(float)nd1->xt.y());

    float factor0 = sin(phi)*(nd0_vec-nd_vec).norm();
    float factor1 = sin(theta)*(nd1_vec-nd_vec).norm();
    float whereToSplit = factor1/(factor0+factor1);  // ~1 means the split is near nd0, ~0 means it is near nd1

    if(whereToSplit < fracture_epsilon && !e1.isBoundary)
    {
        splitEdge = e1;
    }
    else if(whereToSplit > 1-fracture_epsilon && !e0.isBoundary)
    {
        splitEdge = e0;
    }
    else
    {
        if(e_opposite.isBoundary)
        {
            CarefulSplitBoundaryElem(elem, nd, nd0, nd1, whereToSplit, splitEdge,prms);
        }
        else
        {
            icy::Element *elem_adj = e_opposite.getOtherElement(elem);
            CarefulSplitNonBoundaryElem(elem, elem_adj, nd, nd0, nd1, whereToSplit, splitEdge,prms);
        }

    }
    splitEdge.toSplit = true;
    splitEdge.elems[0]->edges[splitEdge.edge_in_elem_idx[0]]=splitEdge;
    splitEdge.elems[1]->edges[splitEdge.edge_in_elem_idx[1]]=splitEdge;
}

void icy::Geometry::CarefulSplitNonBoundaryElem(Element *originalElem, Element *adjElem,
                                 Node *nd, Node *nd0, Node *nd1, float where, Edge &insertedEdge, SimParams &prms)
{
    short ndIdx_orig = originalElem->getNodeIdx(nd);
    short nd0Idx_orig = originalElem->getNodeIdx(nd0);
    short nd1Idx_orig = originalElem->getNodeIdx(nd1);

    Node *oppositeNode = adjElem->getOppositeNode(originalElem->edges[ndIdx_orig]);
    short nd0Idx_adj = adjElem->getNodeIdx(nd0);
    short nd1Idx_adj = adjElem->getNodeIdx(nd1);
    short oppIdx_adj = adjElem->getNodeIdx(oppositeNode);

    Element *insertedFace = AddElement();
    nd->adjacent_elems.push_back(insertedFace);

    Element *insertedFace_adj = AddElement();
    oppositeNode->adjacent_elems.push_back(insertedFace_adj);

    Node *split=AddNode();
    split->InitializeFromAdjacent(nd0, nd1, where);
    split->isBoundary=false;
    split->adjacent_elems.push_back(originalElem);
    split->adjacent_elems.push_back(insertedFace);
    split->adjacent_elems.push_back(adjElem);
    split->adjacent_elems.push_back(insertedFace_adj);

    originalElem->nds[nd1Idx_orig] = split;
    insertedFace->nds[ndIdx_orig] = nd;
    insertedFace->nds[nd1Idx_orig] = nd1;
    insertedFace->nds[nd0Idx_orig] = split;
    insertedFace->edges[nd0Idx_orig] = originalElem->edges[nd0Idx_orig];

    adjElem->nds[nd1Idx_adj] = split;
    insertedFace_adj->nds[oppIdx_adj] = oppositeNode;
    insertedFace_adj->nds[nd1Idx_adj] = nd1;
    insertedFace_adj->nds[nd0Idx_adj] = split;
    insertedFace_adj->edges[nd0Idx_adj] = adjElem->edges[nd0Idx_adj];

    insertedEdge = Edge(nd, split);
    insertedEdge.AddElement(insertedFace, nd1Idx_orig);
    insertedEdge.AddElement(originalElem, nd0Idx_orig);
    insertedEdge.isBoundary = false;
    insertedEdge.toSplit = true;

    Edge insertedEdge_adj = Edge(oppositeNode, split);
    insertedEdge_adj.AddElement(insertedFace_adj, nd1Idx_adj);
    insertedEdge_adj.AddElement(adjElem, nd0Idx_adj);
    insertedEdge_adj.isBoundary = false;
    insertedEdge_adj.toSplit = false;
    insertedFace_adj->edges[nd1Idx_adj] = insertedEdge_adj;
    adjElem->edges[nd0Idx_adj] = insertedEdge_adj;

    Edge exteriorEdge1 = Edge(split, nd1);
    exteriorEdge1.isBoundary = false;
    exteriorEdge1.AddElement(insertedFace, ndIdx_orig);
    exteriorEdge1.AddElement(insertedFace_adj, oppIdx_adj);
    insertedFace->edges[ndIdx_orig] = exteriorEdge1;
    insertedFace_adj->edges[oppIdx_adj] = exteriorEdge1;

    Edge exteriorEdge2 = Edge(split, nd0);
    exteriorEdge2.isBoundary = false;
    exteriorEdge2.AddElement(originalElem, ndIdx_orig);
    exteriorEdge2.AddElement(adjElem, oppIdx_adj);
    originalElem->edges[ndIdx_orig] = exteriorEdge2;
    adjElem->edges[oppIdx_adj] = exteriorEdge2;

    originalElem->InitializePersistentVariables();
    insertedFace->InitializePersistentVariables();
    adjElem->InitializePersistentVariables();
    insertedFace_adj->InitializePersistentVariables();

    originalElem->PrecomputeStiffnessMatrix(prms, elasticityMatrix, D_mats);
    insertedFace->PrecomputeStiffnessMatrix(prms, elasticityMatrix, D_mats);
    adjElem->PrecomputeStiffnessMatrix(prms, elasticityMatrix, D_mats);
    insertedFace_adj->PrecomputeStiffnessMatrix(prms, elasticityMatrix, D_mats);
}



void icy::Geometry::CarefulSplitBoundaryElem(Element *originalElem, Node *nd,
                                             Node *nd0, Node *nd1, float where, Edge &insertedEdge,SimParams &prms)
{
    short ndIdx = originalElem->getNodeIdx(nd);
    short nd0Idx = originalElem->getNodeIdx(nd0);
    short nd1Idx = originalElem->getNodeIdx(nd1);

    Element *insertedFace = AddElement();
    nd->adjacent_elems.push_back(insertedFace);

    Node *split=AddNode();
    split->InitializeFromAdjacent(nd0, nd1, where);
    split->isBoundary=true;
    split->adjacent_elems.push_back(originalElem);
    split->adjacent_elems.push_back(insertedFace);

    originalElem->nds[nd1Idx] = split;
    insertedFace->nds[ndIdx] = nd;
    insertedFace->nds[nd1Idx] = nd1;
    insertedFace->nds[nd0Idx] = split;
    insertedFace->edges[nd0Idx] = originalElem->edges[nd0Idx];

    insertedEdge = Edge(nd, split);
    insertedEdge.AddElement(insertedFace, nd1Idx);
    insertedEdge.AddElement(originalElem, nd0Idx);
    insertedEdge.isBoundary = false;
    insertedEdge.toSplit = true;

//    insertedFace->edges[nd1Idx] = insertedEdge;
//    originalElem->edges[nd0Idx] = insertedEdge;


    Edge exteriorEdge1 = Edge(split, nd1);
    exteriorEdge1.isBoundary = true;
    exteriorEdge1.AddElement(insertedFace, ndIdx);
    insertedFace->edges[ndIdx] = exteriorEdge1;

    Edge exteriorEdge2 = Edge(split, nd0);
    exteriorEdge2.isBoundary = true;
    exteriorEdge2.AddElement(originalElem, ndIdx);
    originalElem->edges[ndIdx] = exteriorEdge2;

    originalElem->InitializePersistentVariables();
    insertedFace->InitializePersistentVariables();

    originalElem->PrecomputeStiffnessMatrix(prms, elasticityMatrix, D_mats);
    insertedFace->PrecomputeStiffnessMatrix(prms, elasticityMatrix, D_mats);
}

