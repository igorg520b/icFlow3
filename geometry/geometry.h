namespace icy { class Geometry; class Node; class Element; class Edge;}
#if !defined(Q_MOC_RUN) // MOC has a glitch when parsing TBB headers

#ifndef FL333_H
#define FL333_H

#include <gmsh.h>
#include <boost/pool/object_pool.hpp>

#include <vector>
#include <memory>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <set>
#include <utility>
#include <algorithm>
#include <chrono>

#include <QObject>
#include <QString>
#include <QDebug>

#include "linearsystem.h"
#include "element.h"
#include "edge.h"

#include <concurrent_unordered_map.h>

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
    long CreateEdges2();     // from the list of elements, infer inner edges and boundary
    long IdentifyDisconnectedRegions(); // used to deal with small fragments
    std::vector<std::pair<double, unsigned>> regions; // first is area, second is number of elements

    unsigned getElemCount() {return elems->size();}
    unsigned getNodeCount() {return nodes->size();}
    unsigned getFreeNodeCount() {
        return std::count_if(nodes->begin(), nodes->end(), [](Node* &nd){return !nd->prescribed;}); }

    void EvaluateStresses(SimParams &prms);     // needed for ComputeFractureDirections
    void DistributeStresses();                  // needed for visualization
    long ComputeFractureDirections(SimParams &prms, double timeStep = 0, bool startingFracture = false); // sets maxNode to breakable node
    long SplitNode(SimParams &prms);   // split the node with the highest normal traction

    // save/load
    void WriteToSerializationBuffers();
    void RestoreFromSerializationBuffers();
    std::vector<double> node_buffer;
    std::vector<int> elems_buffer;

    // fracture
    std::vector<Node*> breakable_range, neighbors_of_crack_tip, local_support;
    icy::Node *maxNode = nullptr;

    icy::Edge getEdgeByNodalIdx(int idx1, int idx2);
    tbb::concurrent_unordered_map<uint64_t, icy::Edge> edges_map2;

private:
    void RecomputeElasticityMatrix(SimParams &prms);
    // matrices for elements, recomputed when rho/Y/nu change
    Eigen::Matrix3d elasticityMatrix;        // this has to be pre-computed whenever Y and nu change
    Eigen::Matrix2d D_mats;

    void SwapCurrentAndTmp();
    void ResizeNodes(std::size_t newSize);
    void ResizeElems(std::size_t newSize);

    icy::Node* AddNode(icy::Node *otherNd=nullptr);
    icy::Element* AddElement();
    // returns newly inserted face for the fan
    void SplitEdge(Edge edge, double where,
                   Node *centerNode, Node* &splitNode, bool forwardDirection, SimParams &prms,
                   Eigen::Vector2f dir);

    void SplitAlongExistingEdge(Edge edge, Node *centerNode, Node* &splitNode,
                                int oppositeNodeIdx, bool forwardDirection, Eigen::Vector2f dir);

    void MeshingStepTwo(double CharacteristicLengthMax);

    boost::object_pool<icy::Node> pool_nodes{10000, 0};
    boost::object_pool<icy::Element> pool_elems{10000, 0};

    void CreateSupportRange(int neighborLevel, std::vector<Node*> &initial_set);
    std::unique_ptr<std::unordered_set<Node*>> tmp_range0 = std::make_unique<std::unordered_set<Node*>>();
    std::unique_ptr<std::unordered_set<Node*>> tmp_range1 = std::make_unique<std::unordered_set<Node*>>();

    std::vector<Element*> wave; // used by IdentifyDisconnectedRegions()

};
#endif
#endif // Q_MOC_RUN
