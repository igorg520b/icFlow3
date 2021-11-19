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

#include <boost/container/small_vector.hpp>

namespace icy { class Node; class SimParams; class Edge; class Element; class Geometry; }

class icy::Node
{
public:
    Node() { Reset(); };
    ~Node() = default;
    Node& operator=(Node&) = delete;

    void Reset();
    void Initialize(double x, double y);
    void Initialize(const Node *other);
    void InitializeLERP(const Node *nd0, const Node *nd1, double f);    // linear interpolaiton between two other nodes

    constexpr static int NumberOfSerializedFields = 19; // number of double values saved to file

    int lsId;      // squential number in the system of equations (-1 if prescribed)
    int locId;          // sequential number in a given floe
    double area;        // mass that the node "represents", for applying various forces
    double vertical_force; // for testing
    bool isBoundary;

    struct Sector
    {
        icy::Element *face;
        icy::Node* nd[2];
        double centerAngle; // angle from the node to the center of the adjacent element
        double angle0, angle1;
        Eigen::Vector2d u_normalized, v_normalized, u_p, v_p; // u is CW, v is CCW
        Eigen::Vector2d t0_top, t1_top, t0_bottom, t1_bottom;
        bool operator<(const Sector& other) const {return centerAngle < other.centerAngle;}
    };

    // separation stress
    struct SepStressResult
    {
        double angle_fwd, angle_bwd;
        icy::Element* faces[2];
        Eigen::Vector2d traction_top[2], traction_bottom[2];
        Eigen::Vector2d tn, tn_p;
        double phi[2], theta[2];
        double trac_normal_top, trac_tangential_top, trac_normal_bottom, trac_tangential_bottom, trac_normal_max;
        double angle0[2], angle1[2];
        double sectorSpan(const int idx) const {return theta[idx]+phi[idx];}
    };

    boost::container::small_vector<icy::Element*, 8> adj_elems;
    boost::container::small_vector<icy::Node::Sector,8> fan;

    double fan_angle_span;  // assigned in UpdateFan();
    bool isCrackTip;
    SepStressResult result_with_max_traction;
    Eigen::Vector2d dir;
    Eigen::Vector2d weakening_direction;    // used when isCrackTip==true
    double max_normal_traction;
    double time_loaded_above_threshold;
    bool potentially_can_fracture;

    void CreateUnrotatedFan();
    void PrepareFan();  // performed when topology changes
    void PrintoutFan(); // for testing
    void ComputeFanVariables(const SimParams &prms);
    void ReplaceAdjacentElement(Element *originalElem, Element *replacement);

    static uint64_t make_key(Node *nd0, Node *nd1); // return unique id for a segment defined by two nodes


    // initial configuration
    Eigen::Matrix<double,3,1> x_initial;

    // at step n: displacement, position, velocity, acceleration
    Eigen::Matrix<double,LinearSystem::DOFS,1> un, xn, vn, an;
    Eigen::Matrix<double,LinearSystem::DOFS,1> ut, xt, vt, at; // at step n+1

    Eigen::Vector3d xn3() { return xn.block(0,0,3,1); }
    Eigen::Vector2d xt2() { return xt.block(0,0,2,1); }

    // forces per node (gravity and any test forces)
    void ComputeElasticForce(LinearSystem &ls, SimParams &prms, double timeStep, double totalTime);
    void AcceptTentativeValues();

    // visualized stress values, distributed from elements in Element::DistributeStresses()
    double str_b[3], str_m[3], str_b_top[3], str_b_bottom[3];
    double str_s[2];

    Eigen::Vector3d normal_n;   // averaged normal of the surrounding elements


private:
    void UpdateFan();   // performed when tentative displacements and stress distribution change; invoked from ComputeFanVariables()
    double NormalTraction(double angle_fwd, double weakening_coeff) const;
    void EvaluateTractions(double angle_fwd, SepStressResult &ssr, const double weakening_coeff) const;



    static double Smoothstep(double edge0, double edge1, double x);
    static double SmoothstepDeriv(double edge0, double edge1, double x);
    static double WaterLine(double x, double y, double t, SimParams &prms);
    static double WaterLineDt(double x, double y, double t, SimParams &prms); // derivative with respect to time

    static double BellShapedPolynomial(double x);
    static double BellShapedPolynomialDx(double x);
};

#endif // NODE_H
#endif // Q_MOC_RUN
