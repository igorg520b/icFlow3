#if !defined(Q_MOC_RUN) // MOC has a glitch when parsing TBB headers

#ifndef FL333_H
#define FL333_H

#include <gmsh.h>

#include <vector>
#include <memory>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <set>
#include <utility>
#include <algorithm>
#include <chrono>
#include <tuple>

#include <QObject>
#include <QString>
#include <QDebug>

#include "linearsystem.h"
#include "element.h"
#include "edge.h"
#include "SimpleObjectPool.h"

#include <set>
#include <concurrent_unordered_map.h>

namespace icy { class Geometry; class Node; class Element; class Edge;}

class icy::Geometry
{
public:
    Geometry();

    // nodes, elements and edges
    std::unique_ptr<std::vector<icy::Node*>> nodes = std::make_unique<std::vector<icy::Node*>>();
    std::unique_ptr<std::vector<icy::Node*>> nodes2 = std::make_unique<std::vector<icy::Node*>>();
    std::unique_ptr<std::vector<icy::Element*>> elems = std::make_unique<std::vector<icy::Element*>>();
    std::unique_ptr<std::vector<icy::Element*>> elems2 = std::make_unique<std::vector<icy::Element*>>();

    double length, width, area;

    // at the "setting up the scene" stage - remesh 2d floe if needed
    void Reset();
    void ImportFloePatch(QString fileName, double CharacteristicLengthMax);
    void Remesh(double CharacteristicLengthMax);

    void PrecomputePersistentVariables(SimParams &prms);
    void AssignLsIds();
    void CreateEdges2();     // from the list of elements, infer inner edges and boundary
    long IdentifyDisconnectedRegions(); // used to deal with small fragments
//    long RemoveDegenerateFragments();
    std::vector<std::tuple<unsigned, double, unsigned>> regions; // region#, area, element count

    unsigned getElemCount() {return elems->size();}
    unsigned getNodeCount() {return nodes->size();}
    unsigned getFreeNodeCount() {
        return std::count_if(nodes->begin(), nodes->end(), [](Node* &nd){return !nd->prescribed;}); }

    void EvaluateStresses(SimParams &prms);     // needed for ComputeFractureDirections
    void DistributeStresses();                  // needed for visualization
    long ComputeFractureDirections(SimParams &prms, double timeStep = 0, bool startingFracture = false); // sets maxNode to breakable node

    long SplitNodeAlt(SimParams &prms);

    // save/load
    void WriteToSerializationBuffers();
    void RestoreFromSerializationBuffers();
    std::vector<double> node_buffer;
    std::vector<int> elems_buffer;

    // fracture
    std::vector<Node*> breakable_range, neighbors_of_crack_tip, local_support;
    std::vector<Element*> local_elems; // elems corresponding to breakable_range;
    std::unordered_set<Element*> local_elems_set; // for computing local_elems
    icy::Node *maxNode = nullptr;

    Edge getEdgeByNodalIdx(int idx1, int idx2);
    tbb::concurrent_unordered_map<uint64_t, Edge> edges_map2;
    std::vector<Edge> allEdges;
    std::vector<Edge> boundaryEdges;

private:
    void RecomputeElasticityMatrix(SimParams &prms);
    // matrices for elements, recomputed when rho/Y/nu change
    Eigen::Matrix3d elasticityMatrix;        // this has to be pre-computed whenever Y and nu change
    Eigen::Matrix2d D_mats;

    void SwapCurrentAndTmp();
    void ResizeNodes(std::size_t newSize);
    void ResizeElems(std::size_t newSize);

    icy::Node* AddNode(icy::Node *otherNd=nullptr);
    icy::Element* AddElement(); // makes a new element

    std::set<Element*> affected_elements_during_split; // a list of elements that were affected by SplitNode
    std::set<Node*> affected_nodes_during_split;
    void UpdateEdges();

    void EstablishSplittingEdge(Edge &splitEdge, Node* nd,
                                const float phi, const float theta, const float fracture_epsilon,
                                const Edge e0, const Edge e1, const Edge e_opposite, Element *elem, SimParams &prms);
    void Fix_X_Topology(Node *nd);
    // preserve boundaries and orientation
    void CarefulSplitBoundaryElem(Element *originalElem, Node *nd, Node *nd0, Node *nd1, float where, Edge &insertedEdge, SimParams &prms);
    void CarefulSplitNonBoundaryElem(Element *originalElem, Element *adjElem, Node *nd,
                                     Node *nd0, Node *nd1, float where, Edge &insertedEdge, SimParams &prms);

    void MeshingStepTwo(double CharacteristicLengthMax);

    icy::SimpleObjectPool<Node> s_pool_nodes;
    icy::SimpleObjectPool<Element> s_pool_elems;

    void CreateSupportRange(int neighborLevel, std::vector<Node*> &initial_set);
    std::unique_ptr<std::unordered_set<Node*>> tmp_range0 = std::make_unique<std::unordered_set<Node*>>();
    std::unique_ptr<std::unordered_set<Node*>> tmp_range1 = std::make_unique<std::unordered_set<Node*>>();

    std::vector<Element*> wave; // used by IdentifyDisconnectedRegions()
//    void RemoveRegion(unsigned idx);

};
#endif
#endif // Q_MOC_RUN
