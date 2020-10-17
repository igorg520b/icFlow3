#ifndef NODE_H
#define NODE_H

#include <vector>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <cmath>
#include <Eigen/Core>
#include "linearsystem.h"
#include "edge.h"

namespace icy { class Node; class SimParams; class Edge; class Element; }

class icy::Node
{
public:
    Node();
    void Reset();
    void InitializeFromAdjacent(icy::Node *nd0, icy::Node *nd1, double f);
    void InitializeFromAnother(icy::Node *nd);

    static const int NumberOfSerializedFields = 19; // number of double values saved to file

    int lsId = -1;      // squential number in the system of equations (-1 if prescribed)
    int locId;          // sequential number in a given floe
    bool prescribed = false;    // motion of prescribed nodes is known up front
    bool isBoundary;
    double area;        // mass that the node "represents", for applying various forces
    double vertical_force; // for testing

    // TODO: combine adjacent_nodes and adjacent_edges ?
    struct Sector
    {
        Sector(icy::Element *elem, icy::Node *nd);
        icy::Element *face;
        float centerAngle; // angle from the node to the center of the adjacent element
        icy::Node* nd[2];
        icy::Edge e[2]; // begins with CW boundary; ends with CCW boundary
        float angle0, angle1, angle_span;
        Eigen::Vector2f u_normalized, v_normalized, u_p, v_p;
        Eigen::Vector2f t0_top, t1_top, t0_bottom, t1_bottom;
    };

    std::vector<icy::Node*> adjacent_nodes;
    std::unordered_map<int, icy::Edge> adjacent_edges_map;
    std::vector<icy::Node::Sector> fan;

    // initial configuration
    Eigen::Matrix<double,DOFS,1> x_initial;

    // at step n: displacement, position, velocity, acceleration
    Eigen::Matrix<double,DOFS,1> un, xn, vn, an;
    Eigen::Matrix<double,DOFS,1> ut, xt, vt, at; // at step n+1

    // forces per node (gravity and any test forces)
    void ComputeElasticForce(SimParams &prms, double timeStep, double totalTime);
    void Assemble(icy::LinearSystem &ls) const;   // distributed F and dF into linear system
    void AcceptTentativeValues();
    Eigen::Matrix<double,DOFS,1> F;
    Eigen::Matrix<double,DOFS,DOFS> dF;

    // visualized values, distributed from elements in Element::DistributeStresses()
    Eigen::Vector3d str_b, str_m, str_b_top, str_b_bottom;
    Eigen::Vector2d str_s, str_b_top_principal, str_b_bottom_principal;



    Eigen::Vector3d normal_n;   // averaged normal of the surrounding elements

    // set the size and initialize with adjacent elements
    void PrepareFan();  // performed when topology changes
    void InitializeFan(); // performed when tentative displacements and stress distribution change
    float fan_angle_span;  // assigned in InitializeFan();

    // separation stress
    struct SepStressResult
    {
        float angle_fwd, angle_bwd;
        icy::Element* faces[2];
        Eigen::Vector2f traction_top[2], traction_bottom[2];
        Eigen::Vector2f tn, tn_p;
        float phi[2], theta[2];
        float trac_normal_top, trac_tangential_top, trac_normal_bottom, trac_tangential_bottom, trac_normal_max;
        icy::Edge e[4];
    };

    void evaluate_tractions(float angle_fwd, SepStressResult &ssr, const float weakening_coeff) const;
    float normal_traction(float angle_fwd, float weakening_coeff) const;

    void ComputeFanVariablesAlt(SimParams &prms);     // compute tractions - alt version
    SepStressResult result_with_max_traction;
    Eigen::Vector2f dir;
    Eigen::Vector2f weakening_direction;    // only used if crack_tip==true
    float max_normal_traction;

//    void ComputeFanVariables(SimParams &prms);     // compute tractions

    // additional fracture parameters
    bool crack_tip, core_node, support_node;
    double timeLoadedAboveThreshold;

    static double OceanWave(double x, double y, double t);
    static double Smoothstep(double edge0, double edge1, double x);
    static double SmoothstepDeriv(double edge0, double edge1, double x);
    static double WaterLine(double x, double y, double t, SimParams &prms);
    static double WaterLineDt(double x, double y, double t, SimParams &prms); // derivative with respect to time

private:
    int idxSepStressResult;
    static const int num_disc = 200;
    SepStressResult sep_stress_results[num_disc];
};

#endif // NODE_H
