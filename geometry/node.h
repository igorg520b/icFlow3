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

    std::vector<icy::Node*> adjacent_nodes;
    std::unordered_map<int, icy::Edge*> adjacent_edges_map;

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

    struct Sector
    {
        Sector(icy::Element *elem, icy::Node *nd);
        icy::Element *face;
        double centerAngle; // angle from the node to the center of the adjacent element
        icy::Node* nd[2];
        icy::Edge* e[2]; // begins with CW boundary; ends with CCW boundary
        double angle0, angle1, angle_span;
        Eigen::Vector3d u_normalized, v_normalized, u_p, v_p;
        Eigen::Vector3d t0_top, t1_top, t0_bottom, t1_bottom;
    };

    std::vector<Sector> fan;
    Eigen::Vector3d normal_n;   // averaged normal of the surrounding elements

    // set the size and initialize with adjacent elements
    void PrepareFan();  // performed when topology changes
    void InitializeFan(); // performed when tentative displacements and stress distribution change
    double fan_angle_span;  // assigned in InitializeFan();

    // separation stress
    struct SepStressResult
    {
        double angle_fwd, angle_bwd;
        icy::Element* faces[2];
        Eigen::Vector3d traction_top[2], traction_bottom[2];
        double phi[2];
        double theta[2];
        double trac_normal_top, trac_tangential_top;
        double trac_normal_bottom, trac_tangential_bottom;
        double trac_normal_max;
        Eigen::Vector3d tn, tn_p;
        icy::Edge* e[4];
    };

    void evaluate_tractions(double angle_fwd, SepStressResult &ssr, const double weakening_coeff) const;
    double normal_traction(double angle_fwd, double weakening_coeff) const;

    void ComputeFanVariablesAlt(SimParams &prms);     // compute tractions - alt version
    SepStressResult result_with_max_traction;
    Eigen::Vector3d dir;
    double max_normal_traction;

    void ComputeFanVariables(SimParams &prms);     // compute tractions

    static const int num_disc = 200;
    SepStressResult sep_stress_results[num_disc];
    int idxSepStressResult;

    // additional fracture parameters
    bool crack_tip, core_node, support_node;
    double timeLoadedAboveThreshold;
    Eigen::Vector3d weakening_direction;

    static double OceanWave(double x, double y, double t);
    static double Smoothstep(double edge0, double edge1, double x);
    static double SmoothstepDeriv(double edge0, double edge1, double x);
    static double WaterLine(double x, double y, double t, SimParams &prms);
    static double WaterLineDt(double x, double y, double t, SimParams &prms); // derivative with respect to time
};

#endif // NODE_H
