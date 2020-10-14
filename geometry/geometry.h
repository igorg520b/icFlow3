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

namespace icy { class Geometry; class Node; class Element; class Edge;}

class icy::Geometry : public QObject
{
    Q_OBJECT
    Q_PROPERTY(int in_Elems READ getElemCount)
    Q_PROPERTY(int in_Nodes READ getNodeCount)
    Q_PROPERTY(int in_FreeNds READ getFreeNodeCount)
    Q_PROPERTY(double in_length MEMBER length NOTIFY propertyChanged)
    Q_PROPERTY(double in_width MEMBER width NOTIFY propertyChanged)
    Q_PROPERTY(double in_area MEMBER area NOTIFY propertyChanged)

public:
    Geometry();

    // nodes, elements and edges
    std::unique_ptr<std::vector<icy::Node*>> nodes = std::make_unique<std::vector<icy::Node*>>();
    std::unique_ptr<std::vector<icy::Node*>> nodes2 = std::make_unique<std::vector<icy::Node*>>();
    std::unique_ptr<std::vector<icy::Element*>> elems = std::make_unique<std::vector<icy::Element*>>();
    std::unique_ptr<std::vector<icy::Element*>> elems2 = std::make_unique<std::vector<icy::Element*>>();
    std::unique_ptr<std::vector<icy::Edge*>> edges = std::make_unique<std::vector<icy::Edge*>>();

    double length, width, area;

    // at the "setting up the scene" stage - remesh 2d floe if needed
    void Reset();
    void ImportFloePatch(QString fileName, double CharacteristicLengthMax);
    void Remesh(double CharacteristicLengthMax);

    void PrecomputePersistentVariables(SimParams &prms);
    void AssignLsIds();
    long CreateEdges();     // from the list of elements, infer inner edges and boundary

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

    // mapping (int,int) -> Edge
    struct EdgeElementPairing
    {
        icy::Element *element0 = nullptr, *element1 = nullptr;
        icy::Edge *edge = nullptr;
        EdgeElementPairing(icy::Element *e0, icy::Element *e1) : element0(e0), element1(e1) {}
        EdgeElementPairing() {}
    };

    // fracture
    std::vector<Node*> breakable_range, neighbors_of_crack_tip, local_support;

    icy::Node *maxNode = nullptr;

    icy::Edge* getEdgeByNodalIdx(int idx1, int idx2);
    //    std::map<std::pair<int,int>, EdgeElementPairing> edges_map;
    std::unordered_map<uint64_t, EdgeElementPairing> edges_map;

private:
    void RecomputeElasticityMatrix(SimParams &prms);
    // matrices for elements, recomputed when rho/Y/nu change
    Eigen::Matrix3d elasticityMatrix;        // this has to be pre-computed whenever Y and nu change
    Eigen::Matrix2d D_mats;

    void SwapCurrentAndTmp();
    void ResizeNodes(std::size_t newSize);
    void ResizeElems(std::size_t newSize);
    void ResizeEdges(std::size_t newSize);

    icy::Node* AddNode(icy::Node *otherNd=nullptr);
    icy::Element* AddElement();
    // returns newly inserted face for the fan
    void SplitEdge(icy::Edge *edge, double where,
                   icy::Node *centerNode, icy::Node* &splitNode, bool forwardDirection, SimParams &prms,
                   Eigen::Vector2f dir);

    void SplitAlongExistingEdge(Edge *edge, Node *centerNode, Node* &splitNode,
                                int oppositeNodeIdx, bool forwardDirection, Eigen::Vector2f dir);

    void MeshingStepTwo(double CharacteristicLengthMax);

    boost::object_pool<icy::Node> pool_nodes{10000, 0};
    boost::object_pool<icy::Element> pool_elems{10000, 0};
    boost::object_pool<icy::Edge> pool_edges{10000, 0};

    void CreateSupportRange(int neighborLevel, std::vector<Node*> &initial_set);
    std::unique_ptr<std::unordered_set<Node*>> tmp_range0 = std::make_unique<std::unordered_set<Node*>>();
    std::unique_ptr<std::unordered_set<Node*>> tmp_range1 = std::make_unique<std::unordered_set<Node*>>();


signals:
    void propertyChanged();
};
#endif
