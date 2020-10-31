#include "geometry.h"
#include "model.h"
#include <numeric>
#include <algorithm>
#include <iterator>


long icy::Geometry::ComputeFractureDirections(SimParams &prms, double timeStep, bool startingFracture)
{

    auto t1 = std::chrono::high_resolution_clock::now();

    maxNode=nullptr;
    double temporal_attenuation = prms.temporal_attenuation;

    float threashold = prms.normal_traction_threshold;

    std::size_t nNodes = nodes->size();


    if(startingFracture)
    {
        // evaluate all nodes to compute breakable range
        breakable_range_concurrent.clear();

        EvaluateStresses(prms, (*elems));
        DistributeStresses();

//#pragma omp parallel for
        for(std::size_t i=0;i<nNodes;i++)
        {
            icy::Node *nd = (*nodes)[i];
            if(nd->potentially_can_fracture)
            {
                if(nd->timeLoadedAboveThreshold >= temporal_attenuation)
                {
                    nd->PrepareFan2();
                    nd->ComputeFanVariablesAlt(prms);
                    nd->VerifySSR();
                    if(nd->max_normal_traction > threashold) breakable_range_concurrent.push_back(nd);
                    else nd->timeLoadedAboveThreshold = 0;
                }
                else
                {
                    nd->timeLoadedAboveThreshold+=timeStep;
                    nd->max_normal_traction = 0;
                    nd->dir=Eigen::Vector2f::Zero();
                }
            }
            else
            {
                nd->timeLoadedAboveThreshold = 0;
                nd->max_normal_traction = 0;
                nd->dir=Eigen::Vector2f::Zero();
            }
        }
        breakable_range.clear();
        std::copy(breakable_range_concurrent.begin(), breakable_range_concurrent.end(), std::back_inserter(breakable_range));
        std::sort(breakable_range.begin(), breakable_range.end(), [](Node *nd1, Node *nd2)
        {return nd1->max_normal_traction > nd2->max_normal_traction;});

        const unsigned max_breakable_range = 100;
        if(breakable_range.size() > max_breakable_range) breakable_range.resize(max_breakable_range);
    }
    else
    {

        // insert the recently created crack tips into the breakable range
        for(Node *nct : new_crack_tips)
        {
            nct->ComputeFanVariablesAlt(prms);
            auto find_result = std::find(breakable_range.begin(), breakable_range.end(),nct);
            bool already_contains = find_result!=breakable_range.end();

            if(nct->max_normal_traction > threashold && !already_contains)
                breakable_range.push_back(nct);
            nct->timeLoadedAboveThreshold = temporal_attenuation;
        }
        new_crack_tips.clear();

        // remove the nodes that were affected by the crack on the previous step
        breakable_range.erase(std::remove_if(breakable_range.begin(), breakable_range.end(),
                                          [temporal_attenuation](Node *nd)
                {return nd->max_normal_traction==0 || (nd->timeLoadedAboveThreshold < temporal_attenuation && !nd->crack_tip);}),
                breakable_range.end());

        // update Sector in case if topology changed around this node
        for(Node *nd : breakable_range)
        {
            nd->PrepareFan2();
            nd->ComputeFanVariablesAlt(prms);
            if(nd->max_normal_traction>0) {
                std::cout << "VerifyFan1" <<std::endl;
                nd->VerifyFan();
                nd->VerifySSR();
            }
        }

    }


    if(breakable_range.size() > 0)
    {
        // take out maximal node from breakable_range
        auto it_nd = std::max_element(breakable_range.begin(), breakable_range.end(),
                                      [](Node *nd1, Node *nd2) {
                if(nd2->crack_tip && nd2->max_normal_traction>0 && !nd1->crack_tip) return true;
                return nd1->max_normal_traction < nd2->max_normal_traction; });

        if((*it_nd)->max_normal_traction > 0)
        {
            maxNode = *it_nd;

            // make sure that the Sector information is updated

            maxNode->PrepareFan2();
            maxNode->ComputeFanVariablesAlt(prms);
            maxNode->timeLoadedAboveThreshold = 0;

#ifdef QT_DEBUG
            std::cout << "\n\nselected node " << maxNode->locId << std::endl;
            std::cout << "breakable range " << breakable_range.size() << "\n";
            for(Node *nd : breakable_range)
                std::cout << nd->locId << "; " << nd->max_normal_traction << (nd->crack_tip ? " *" : "") << std::endl;
#endif
            breakable_range.erase(it_nd);

            std::cout << "maxNode->VerifyFan()"<<std::endl;
            maxNode->VerifyFan();
            std::cout << "maxNode->VerifySSR()"<<std::endl;
            maxNode->VerifySSR();
            std::cout << "VerifySSR passed"<<std::endl;
        }
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
}





//========================== ALTERNATIVE ALGORITHM FOR SPLITTING
long icy::Geometry::SplitNodeAlt(SimParams &prms)
{
    auto t1 = std::chrono::high_resolution_clock::now();

    if(maxNode == nullptr) throw std::runtime_error("trying to split nullptr");


    new_crack_tips.clear();
    affected_elements_during_split.clear();

    icy::Node* nd = maxNode;
    std::cout << "verify1 " << std::endl;
    nd->VerifyFan();
    nd->timeLoadedAboveThreshold = 0;
    for(Element *e : nd->adjacent_elems) affected_elements_during_split.insert(e);

    // subsequent calculations are based on the fracture direction where the traction is maximal
    icy::Node::SepStressResult &ssr = nd->result_with_max_traction;

    // make sure that the interior node has two split faces
    bool isBoundary = (ssr.faces[1] == nullptr);
    if(isBoundary != nd->isBoundary) std::runtime_error("isBoundary != nd->isBoundary");
//    std::cout << "\n SplitNode " << nd->locId << (isBoundary ? " boundary" : " nb ") << "\n";

    icy::Edge splitEdge_fw;
    std::cout << "forward\n";
    if(!ssr.faces[0]->ContainsNode(nd)) throw std::runtime_error("SplitNode: mesh toplogy error 0");
    nd->VerifyFan();
    nd->VerifySSR();
    EstablishSplittingEdge(splitEdge_fw, nd,
                                ssr.phi[0], ssr.theta[0], prms.fracture_epsilon,
                                ssr.e[0], ssr.e[1], ssr.faces[0],prms);

    Node *split0=nullptr,*split1=nullptr;
//    std::cout << "split0" << std::endl;
    split0=splitEdge_fw.getOtherNode(nd);
//    std::cout << "split0 done: " << nd->locId << "-" << split0->locId << std::endl;

    icy::Edge splitEdge_bw;
    if(!isBoundary)
    {
        // determine splitEdge_bw (create if necessary)
        std::cout << "backward\n";
        if(!ssr.faces[1]->ContainsNode(nd)) throw std::runtime_error("SplitNode: mesh toplogy error 1");
        EstablishSplittingEdge(splitEdge_bw, nd,
                                    ssr.phi[1], ssr.theta[1], prms.fracture_epsilon,
                                    ssr.e[2], ssr.e[3], ssr.faces[1],prms);

//        std::cout << "split1" << std::endl;
        split1=splitEdge_bw.getOtherNode(nd);
//        std::cout << "split1 done: " << nd->locId << "-" << split1->locId  << (split1->isBoundary ? " b" : " nb" ) << std::endl;
    }


//    std::cout << "Fix_X_Topology(nd)" << std::endl;
    Fix_X_Topology(nd); // split as fan.front().e[0] --- splitEdge_fw --- fan.back().e[1]
//    std::cout << "Fix_X_Topology(nd) done" << std::endl;

    if(!isBoundary)
    {
        if(split1==nullptr) throw std::runtime_error("split1==nullptr");
        if(split1->isBoundary)
        {
//            std::cout << "Fix_X_Topology(split1) " << split1->locId << std::endl;
            Fix_X_Topology(split1);
//            std::cout << "Fix_X_Topology(split1) done" << std::endl;
        }
        else
        {
//            std::cout << "split1" << std::endl;
            split1->crack_tip = true;
            new_crack_tips.push_back(split1);
            split1->weakening_direction = Eigen::Vector2f(split1->xt.x()-nd->xt.x(), split1->xt.y()-nd->xt.y());
            split1->weakening_direction.normalize();
            std::cout << "split1 done" << std::endl;
            for(Element *e : split1->adjacent_elems) affected_elements_during_split.insert(e);
        }
    }

    if(split0==nullptr)  throw std::runtime_error("split0==nullptr");
    if(split0->isBoundary)
    {
//        std::cout << "Fix_X_Topology(split0)" << std::endl;
        Fix_X_Topology(split0);
//        std::cout << "Fix_X_Topology(split0) done" << std::endl;
    }
    else
    {
        std::cout << "split0" << std::endl;
        split0->crack_tip = true;
        new_crack_tips.push_back(split0);
        split0->weakening_direction = Eigen::Vector2f(split0->xt.x()-nd->xt.x(), split0->xt.y()-nd->xt.y());
        split0->weakening_direction.normalize();
        std::cout << "split0 done" << std::endl;
        for(Element *e : split0->adjacent_elems) affected_elements_during_split.insert(e);
    }

//    std::cout << "UpdateEdges" << std::endl;
    UpdateEdges();
//    CreateEdges2();


    nd->weakening_direction = Eigen::Vector2f::Zero();
    nd->crack_tip = false;

//    std::cout << "SplitNode done" << std::endl;
    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
}

void icy::Geometry::Fix_X_Topology(Node *nd)
{
    Eigen::Vector3d nd_vec = nd->x_initial.block(0,0,3,1);

    nd->fan.clear();
    for(unsigned k=0;k<nd->adjacent_elems.size();k++)
    {
        icy::Element *elem = nd->adjacent_elems[k];

        Node::Sector s;
        s.face = elem;
        Eigen::Vector3d tcv = elem->getCenter() - nd_vec;
        s.centerAngle = atan2(tcv.y(), tcv.x());

        short thisIdx, CWIdx, CCWIdx;
        elem->getIdxs(nd, thisIdx, CWIdx, CCWIdx);

        s.nd[0] = elem->nds[CWIdx];
        s.nd[1] = elem->nds[CCWIdx];

        // note that the indices are swapped
        s.e[0] = elem->edges[CCWIdx];
        s.e[1] = elem->edges[CWIdx];
        s.e[2] = elem->edges[thisIdx];
        nd->fan.push_back(s);
    }

    std::sort(nd->fan.begin(), nd->fan.end(),
              [](const Node::Sector &f0, const Node::Sector &f1)
    {return f0.centerAngle < f1.centerAngle; });


    auto cw_boundary = std::find_if(nd->fan.begin(), nd->fan.end(), [](const Node::Sector &f){return f.e[0].isBoundary;});
    if(cw_boundary != nd->fan.end()) {
        std::rotate(nd->fan.begin(), cw_boundary, nd->fan.end());
//        std::cout << "Fix_X_Topology: rotated" << std::endl;
    }


    Node *split=nullptr;
    bool replacing = false;
    for(unsigned i=1;i<nd->fan.size();i++)
    {
        Node::Sector &s = nd->fan[i];
        if(s.e[0].toSplit || s.e[0].isBoundary) replacing=!replacing;
        if(replacing) {
            if(split==nullptr) split=AddNode();
            s.face->ReplaceNode(nd, split);
        }
    }

    if(split!=nullptr) split->InitializeFromAnother(nd);
    else {
        qDebug()<<"nothing to split! " << nd->locId;
    }

    for(Node::Sector &s : nd->fan) affected_elements_during_split.insert(s.face);
//    std::cout << "Fix X\n";
//    nd->PrintoutFan();
    std::cout << endl;
}


void icy::Geometry::EstablishSplittingEdge(Edge &splitEdge, Node* nd,
                            const float phi, const float theta, const float fracture_epsilon,
                            const Edge e0, const Edge e1, Element *elem, SimParams &prms)
{
//    std::cout << "Establsh SplittingEdge" << std::endl;
    icy::Node *nd0, *nd1;
    nd0= e0.getOtherNode(nd);
    nd1= e1.getOtherNode(nd);

    Eigen::Vector2f nd_vec((float)nd->xt.x(),(float)nd->xt.y());
    Eigen::Vector2f nd0_vec((float)nd0->xt.x(),(float)nd0->xt.y());
    Eigen::Vector2f nd1_vec((float)nd1->xt.x(),(float)nd1->xt.y());

    float factor0 = sin(phi)*(nd0_vec-nd_vec).norm();
    float factor1 = sin(theta)*(nd1_vec-nd_vec).norm();
    float whereToSplit = factor1/(factor0+factor1);  // ~1 means the split is near nd0, ~0 means it is near nd1

    if((whereToSplit < 1e-5 && e1.isBoundary) || (whereToSplit > 1-1e-5 && e0.isBoundary))
    {
        std::cout << "about to create a generate element; nd " << nd->locId << std::endl;
        std::cout << "wheretosplit " << whereToSplit << std::endl;
        throw std::runtime_error("degenerate element 1");
        //TODO: just don't split anything if this happens
    }

    if(whereToSplit < fracture_epsilon && !e1.isBoundary)
    {
//        std::cout << nd->locId << " aligning with edge1 " << e1.nds[0]->locId << "-" << e1.nds[1]->locId << '\n';
        splitEdge = e1;
    }
    else if(whereToSplit > 1-fracture_epsilon && !e0.isBoundary)
    {
//        std::cout << nd->locId << "aligning with edge0 " << e0.nds[0]->locId << "-" << e0.nds[1]->locId << '\n';
        splitEdge = e0;
    }
    else
    {
        icy::Element *elem_adj = elem->getAdjacentElementOppositeToNode(nd);
        if(elem_adj==nullptr)
        {
            CarefulSplitBoundaryElem(elem, nd, nd0, nd1, whereToSplit, splitEdge,prms);
        }
        else
        {
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
//    std::cout << "NonBoundaryElem; nd " << nd->locId << "; nd0 " << nd0->locId << "; nd1 " << nd1->locId << '\n';
//    std::cout << "originalElem: " << originalElem->nds[0]->locId << ", "<< originalElem->nds[1]->locId;
//    std::cout << ", " << originalElem->nds[2]->locId << std::endl;
    originalElem->AssertEdges();

    if(!originalElem->ContainsNode(nd)) throw std::runtime_error("originalElem does not contain nd");
    if(!originalElem->ContainsNode(nd0)) throw std::runtime_error("originalElem does not contain nd0");
    if(!originalElem->ContainsNode(nd1)) throw std::runtime_error("originalElem does not contain nd1");
    short ndIdx_orig = originalElem->getNodeIdx(nd);
    short nd0Idx_orig = originalElem->getNodeIdx(nd0);
    short nd1Idx_orig = originalElem->getNodeIdx(nd1);

    Node *oppositeNode = adjElem->getOppositeNode(nd0, nd1);
    if(!adjElem->ContainsNode(oppositeNode)) throw std::runtime_error("adjElem does not contain opposite node");
    if(!adjElem->ContainsNode(nd0)) throw std::runtime_error("adjElem does not contain nd0");
    if(!adjElem->ContainsNode(nd1)) throw std::runtime_error("adjElem does not contain nd1");
    short nd0Idx_adj = adjElem->getNodeIdx(nd0);
    short nd1Idx_adj = adjElem->getNodeIdx(nd1);
    short oppIdx_adj = adjElem->getNodeIdx(oppositeNode);

    Element *insertedFace = AddElement();

    Element *insertedFace_adj = AddElement();

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
    nd->adjacent_elems.push_back(insertedFace);
    nd1->adjacent_elems.push_back(insertedFace);
    split->adjacent_elems.push_back(insertedFace);

    adjElem->nds[nd1Idx_adj] = split;
    insertedFace_adj->nds[oppIdx_adj] = oppositeNode;
    insertedFace_adj->nds[nd1Idx_adj] = nd1;
    insertedFace_adj->nds[nd0Idx_adj] = split;
    insertedFace_adj->edges[nd0Idx_adj] = adjElem->edges[nd0Idx_adj];
    oppositeNode->adjacent_elems.push_back(insertedFace_adj);
    nd1->adjacent_elems.push_back(insertedFace_adj);
    split->adjacent_elems.push_back(insertedFace_adj);

    insertedEdge = Edge(nd, split);
    insertedEdge.AddElement(insertedFace, nd1Idx_orig);
    insertedEdge.AddElement(originalElem, nd0Idx_orig);
    insertedEdge.isBoundary = false;
    insertedEdge.toSplit = true;
//    insertedFace->edges[nd1Idx_orig]=insertedEdge;
//    originalElem->edges[nd0Idx_orig]=insertedEdge;

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

    if(originalElem->normal_initial.z() < 0) throw std::runtime_error("NonBoundaryElem: normal inconsistent in originalElem");
    if(insertedFace->normal_initial.z() < 0) throw std::runtime_error("NonBoundaryElem: normal inconsistent in insertedFace");
    if(adjElem->normal_initial.z() < 0) throw std::runtime_error("NonBoundaryElem: normal inconsistent in adjElem");
    if(insertedFace_adj->normal_initial.z() < 0) throw std::runtime_error("NonBoundaryElem: normal inconsistent in insertedFace_adj");

    originalElem->PrecomputeStiffnessMatrix(prms, elasticityMatrix, D_mats);
    insertedFace->PrecomputeStiffnessMatrix(prms, elasticityMatrix, D_mats);
    adjElem->PrecomputeStiffnessMatrix(prms, elasticityMatrix, D_mats);
    insertedFace_adj->PrecomputeStiffnessMatrix(prms, elasticityMatrix, D_mats);

    affected_elements_during_split.insert(originalElem);
    affected_elements_during_split.insert(insertedFace);
    affected_elements_during_split.insert(adjElem);
    affected_elements_during_split.insert(insertedFace_adj);
}



void icy::Geometry::CarefulSplitBoundaryElem(Element *originalElem, Node *nd,
                                             Node *nd0, Node *nd1, float where, Edge &insertedEdge,SimParams &prms)
{
//    std::cout << "BoundaryElem; nd " << nd->locId << "; nd0 " << nd0->locId << "; nd1 " << nd1->locId << '\n';
//    std::cout << originalElem->nds[0]->locId << " - " << originalElem->nds[1]->locId << " - "<< originalElem->nds[2]->locId << '\n';

    if(!originalElem->ContainsNode(nd)) throw std::runtime_error("originalElem does not contain nd");
    if(!originalElem->ContainsNode(nd0)) throw std::runtime_error("originalElem does not contain nd0");
    if(!originalElem->ContainsNode(nd1)) throw std::runtime_error("originalElem does not contain nd1");
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
//    insertedFace->edges[nd1Idx]=insertedEdge;
//    originalElem->edges[nd0Idx]=insertedEdge;

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

    if(originalElem->normal_initial.z() < 0) throw std::runtime_error("CarefulSplitBoundaryElem: normal inconsistent in originalElem");
    if(insertedFace->normal_initial.z() < 0) throw std::runtime_error("CarefulSplitBoundaryElem: normal inconsistent in insertedFace");

    originalElem->PrecomputeStiffnessMatrix(prms, elasticityMatrix, D_mats);
    insertedFace->PrecomputeStiffnessMatrix(prms, elasticityMatrix, D_mats);

    affected_elements_during_split.insert(originalElem);
    affected_elements_during_split.insert(insertedFace);

//    std::cout << "BoundaryElem replaced with:\n";
//    std::cout << originalElem->nds[0]->locId << " - " << originalElem->nds[1]->locId << " - "<< originalElem->nds[2]->locId << '\n';
//    std::cout << insertedFace->nds[0]->locId << " - " << insertedFace->nds[1]->locId << " - "<< insertedFace->nds[2]->locId << '\n';

}

void icy::Geometry::UpdateEdges()
{
    std::unordered_set<Node*> affected_nodes_during_split; // their neighbors are also affected
    std::unordered_set<Element *> expanded_set_elems1;
    std::unordered_set<Element *> expanded_set_elems2;

    for(Element *elem : affected_elements_during_split)
    {
        for(int k=0;k<3;k++)
        {
            affected_nodes_during_split.insert(elem->nds[k]);
            if(elem->adj_elems[k]!=nullptr) expanded_set_elems1.insert(elem->adj_elems[k]);
            for(Element *elem2 : elem->nds[k]->adjacent_elems) {
                expanded_set_elems1.insert(elem2);
                for(int m=0;m<3;m++) affected_nodes_during_split.insert(elem2->nds[m]);
            }
        }
        expanded_set_elems1.insert(elem);
    }


    for(Node *nd : affected_nodes_during_split)
    {
        for(Element *elem : nd->adjacent_elems) expanded_set_elems2.insert(elem);
        nd->adjacent_elems.clear();
        nd->area=0;
        nd->isBoundary=false;
    }

    for(Element *elem : expanded_set_elems1) {
        for(int k=0;k<3;k++) {
            if(elem->adj_elems[k]!=nullptr) expanded_set_elems2.insert(elem->adj_elems[k]);
            elem->adj_elems[k]=nullptr;
        }
        expanded_set_elems2.insert(elem);
    }



    std::unordered_map<uint64_t, Edge> edges_map;

    for(Element *elem : expanded_set_elems2)
    {
        for(int i=0;i<3;i++)
        {
            Node *nd = elem->nds[i];
            // only work with the nodes in the affected_nodes set
            if(affected_nodes_during_split.find(nd)!=affected_nodes_during_split.end())
            {
                nd->adjacent_elems.push_back(elem);
            }
            // process edges
            int nd0idx = elem->nds[i]->locId;
            int nd1idx = elem->nds[(i+1)%3]->locId;
            if(nd0idx > nd1idx) std::swap(nd0idx, nd1idx);
            uint64_t key = ((uint64_t)nd0idx << 32) | nd1idx;

            icy::Node *nd0 = (*nodes)[nd0idx];
            icy::Node *nd1 = (*nodes)[nd1idx];

            Edge edge(nd0, nd1);

            edges_map.insert({key,edge});
        }
    }

    // note that edges_map may contain edges outside of affected_elements
    for(Element *elem : expanded_set_elems2)
    {
        for(int i=0;i<3;i++)
        {
            int nd0idx = elem->nds[i]->locId;
            int nd1idx = elem->nds[(i+1)%3]->locId;
            if(nd0idx > nd1idx) std::swap(nd0idx, nd1idx);
            uint64_t key = ((uint64_t)nd0idx << 32) | nd1idx;

            icy::Edge &existing_edge = edges_map.at(key);
            existing_edge.AddElement(elem, (i+2)%3);
        }
    }

    std::unordered_map<uint64_t, Edge> correctly_inferred_edges;
    // only take edges in the affected_elements set
    for(Element *elem : expanded_set_elems1)
    {
        for(int i=0;i<3;i++)
        {
            int nd0idx = elem->nds[i]->locId;
            int nd1idx = elem->nds[(i+1)%3]->locId;
            if(nd0idx > nd1idx) std::swap(nd0idx, nd1idx);
            uint64_t key = ((uint64_t)nd0idx << 32) | nd1idx;

            icy::Edge &existing_edge = edges_map.at(key);
            correctly_inferred_edges.insert({key,existing_edge});
        }
    }

//    qDebug() << "affected elems " << affected_elements_during_split.size();
//    qDebug() << "affected nodes " << affected_nodes_during_split.size();
 //   qDebug() << "expanded set of elems " << expanded_set_elems.size();
//    qDebug() << "edges_map " << edges_map.size();
//    qDebug() << "correctly inferred edges " << correctly_inferred_edges.size();


    boundaryEdges.erase(std::remove_if(boundaryEdges.begin(),boundaryEdges.end(),
                   [affected_nodes_during_split](Edge e)
    {return (affected_nodes_during_split.find(e.nds[0])!=affected_nodes_during_split.end() &&
                affected_nodes_during_split.find(e.nds[1])!=affected_nodes_during_split.end());}),
            boundaryEdges.end());

    for(auto kvp : correctly_inferred_edges)
    {
        Edge &existing_edge = kvp.second;
        icy::Element *elem_of_edge0 = existing_edge.elems[0];
        icy::Element *elem_of_edge1 = existing_edge.elems[1];
        short idx0 = existing_edge.edge_in_elem_idx[0];
        short idx1 = existing_edge.edge_in_elem_idx[1];

        if(elem_of_edge0 == nullptr && elem_of_edge1 == nullptr) throw std::runtime_error("disconnected edge?");
        existing_edge.isBoundary = (existing_edge.elems[0] == nullptr || existing_edge.elems[1] == nullptr);

        if(elem_of_edge0 != nullptr) elem_of_edge0->edges[idx0] = existing_edge;
        if(elem_of_edge1 != nullptr) elem_of_edge1->edges[idx1] = existing_edge;

        if(!existing_edge.isBoundary)
        {
            elem_of_edge0->adj_elems[idx0] = elem_of_edge1;
            elem_of_edge1->adj_elems[idx1] = elem_of_edge0;
        }

        if(existing_edge.isBoundary) boundaryEdges.push_back(existing_edge);
    }


    for(Element *elem : expanded_set_elems2) {
        elem->AssertEdges();
    }

//    try {
//    std::cout << "UpdateEdges invoking PrepareFan2: ";
    for(Node *nd : affected_nodes_during_split)
    {
//        std::cout << nd->locId << ", ";
        nd->PrepareFan2();
    }
//    } catch(std::runtime_error &e)
//    {
//        qDebug() << "reverting to createedges";
//        CreateEdges2();
//    }
//    std::cout << '\n' <<  std::endl;


}
