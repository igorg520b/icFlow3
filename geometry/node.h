#if !defined(Q_MOC_RUN) // MOC has a glitch when parsing TBB headers
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

#include <concurrent_vector.h>
#include <concurrent_unordered_map.h>

namespace icy { class Node; class SimParams; class Edge; class Element; class Geometry; }

class icy::Node
{
public:
    Node();
    void Reset();
    void InitializeFromAdjacent(const Node *nd0, const Node *nd1, double f);
    void InitializeFromAnother(const Node *nd);

    static const int NumberOfSerializedFields = 19; // number of double values saved to file

    int lsId = -1;      // squential number in the system of equations (-1 if prescribed)
    int locId;          // sequential number in a given floe
    bool isBoundary;
    double area;        // mass that the node "represents", for applying various forces
    float vertical_force; // for testing

    struct Sector
    {
        // Sector(icy::Element *elem, icy::Node *nd);
        icy::Element *face;
        float centerAngle; // angle from the node to the center of the adjacent element
        icy::Node* nd[2];
        float angle0, angle1;
        Eigen::Vector2f u_normalized, v_normalized, u_p, v_p;
        Eigen::Vector2f t0_top, t1_top, t0_bottom, t1_bottom;
    };

    tbb::concurrent_vector<icy::Element*> adjacent_elems;
    std::vector<icy::Node::Sector> fan;
    void PrintoutFan(); // for testing

    // initial configuration
    Eigen::Matrix<double,3,1> x_initial;

    // at step n: displacement, position, velocity, acceleration
    Eigen::Matrix<double,DOFS,1> un, xn, vn, an;
    Eigen::Matrix<double,DOFS,1> ut, xt, vt, at; // at step n+1

    // forces per node (gravity and any test forces)
    void ComputeElasticForce(LinearSystem &ls, SimParams &prms, double timeStep, double totalTime);
    void AcceptTentativeValues();

    // visualized stress values, distributed from elements in Element::DistributeStresses()
    float str_b[3], str_m[3], str_b_top[3], str_b_bottom[3];
    float str_s[2];

    Eigen::Vector3d normal_n;   // averaged normal of the surrounding elements

    // set the size and initialize with adjacent elements
    void PrepareFan2();  // performed when topology changes
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
        icy::Edge e_opposite[2]; // edges that lie opposite of the center node
    };

    void evaluate_tractions(float angle_fwd, SepStressResult &ssr, const float weakening_coeff) const;
    float normal_traction(float angle_fwd, float weakening_coeff) const;

    void ComputeFanVariablesAlt(SimParams &prms);     // compute tractions - alt version
    SepStressResult result_with_max_traction;
    Eigen::Vector2f dir;
    Eigen::Vector2f weakening_direction;    // only used if crack_tip==true
    float max_normal_traction;

    // additional fracture parameters
    bool crack_tip, reset_timing, support_node, potentially_can_fracture;
    double timeLoadedAboveThreshold;

    static double Smoothstep(double edge0, double edge1, double x);
    static double SmoothstepDeriv(double edge0, double edge1, double x);
    static double WaterLine(double x, double y, double t, SimParams &prms);
    static double WaterLineDt(double x, double y, double t, SimParams &prms); // derivative with respect to time

    static double BellShapedPolynomial(double x);
    static double BellShapedPolynomialDx(double x);
};

#endif // NODE_H
#endif // Q_MOC_RUN
