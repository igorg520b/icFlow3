#include "node.h"
#include "model.h"
#include "parameters_sim.h"
#include "element.h"

#include <cmath>
#include <cfloat>
#include <algorithm>

#include <spdlog/spdlog.h>
#include <boost/math/tools/minima.hpp>


void icy::Node::Reset()
{
    ut.setZero();
    xt.setZero();
    vt.setZero();
    at.setZero();
    un.setZero();
    xn.setZero();
    vn.setZero();
    an.setZero();
    x_initial.setZero();
    normal_n.setZero();

    lsId = locId = -1;
    area = 0;
    vertical_force = 0;

    adj_elems.clear();
    fan.clear();

    isCrackTip = false;
    max_normal_traction = 0;
    dir.setZero();
    weakening_direction.setZero();
    time_loaded_above_threshold = 0;
}


void icy::Node::Initialize(double x, double y)
{
    x_initial << x,y,0;
    xn3() = x_initial;
}

void icy::Node::Initialize(const Node *nd)
{
    x_initial = nd->x_initial;
    ut = nd->ut;
    un = nd->un;
    xt = nd->xt;
    xn = nd->xn;
    vt = nd->vt;
    vn = nd->vn;
    at = nd->at;
    an = nd->an;
}

void icy::Node::InitializeLERP(const Node *nd0, const Node *nd1, double f)
{
    x_initial = nd0->x_initial*f + nd1->x_initial*(1-f);
    ut = nd0->ut*f + nd1->ut*(1-f);
    un = nd0->un*f + nd1->un*(1-f);
    xt = nd0->xt*f + nd1->xt*(1-f);
    xn = nd0->xn*f + nd1->xn*(1-f);
    vt = nd0->vt*f + nd1->vt*(1-f);
    vn = nd0->vn*f + nd1->vn*(1-f);
    at = nd0->at*f + nd1->at*(1-f);
    an = nd0->an*f + nd1->an*(1-f);
}






void icy::Node::CreateUnrotatedFan()
{
    if(adj_elems.size() == 0) throw std::runtime_error("Node::CreateUnorderedFan: disconnected node");
    fan.clear();
    area = 0;
    for(Element *elem : adj_elems)
    {
        area += elem->area_initial/3;
        Node::Sector s;
        s.face = elem;
        Eigen::Vector3d tcv = elem->getCenter() - x_initial;
        s.centerAngle = atan2(tcv.y(), tcv.x());
        std::tie(s.nd[0],s.nd[1]) = elem->CW_CCW_Node(this);
        fan.push_back(s);
    }
    std::sort(fan.begin(), fan.end());  // order by angle of the element, counter-clockwise
}


void icy::Node::PrepareFan()
{
    CreateUnrotatedFan();
    isBoundary = std::any_of(fan.begin(),fan.end(),[this](Sector &s){return s.face->isOnBoundary(this);});

    // if boundary, then ensure that sectors start with a boundary element and end with a boundary element
    if(isBoundary)
    {
        // find the fan element with the border on the CW direction
        auto cw_boundary = std::find_if(fan.begin(), fan.end(), [this](const Sector &f){return f.face->isCWBoundary(this);});
        if(cw_boundary == fan.end()) { PrintoutFan(); throw std::runtime_error("icy::Node::PrepareFan(): cw boundary not found"); }
        std::rotate(fan.begin(), cw_boundary, fan.end());
        if(!fan.back().face->isCCWBoundary(this)) { PrintoutFan(); throw std::runtime_error("PrepareFan(): no CCW boundary"); }
    }

    // assert that the nodes of the fan connect
    for(std::size_t i = 0;i<fan.size()-1;i++)
    {
        if(fan[i].nd[1] != fan[i+1].nd[0])
        {
            spdlog::critical("fan nodes are not contiguous for node {}", locId);
            PrintoutFan();
            throw std::runtime_error("fan nodes are not contiguous");
        }
    }
}



void icy::Node::PrintoutFan()
{
    spdlog::info("Printing fan for node {}; isCrackTip: {}; isBoundary: {}", locId, isCrackTip, isBoundary);
    spdlog::info("fan.size {}; adj_elems.size {}", fan.size(), adj_elems.size());
    spdlog::info("fan_angle_span: {}", fan_angle_span);

    spdlog::info("┌ {0: ^9} ┬ {1: ^14} ┬ {2: ^14} ┬ {3: ^6} ┬ {4: ^10} ┬ {5: ^8} ┬ {6: ^8}",
                 "nd1-nd2", "nds 0-1-2", "angle0-angle1", "cAngle", "area", "CWB", "CCWB");

    for(Sector &s : fan)
        spdlog::info("│ {0:>4}-{1:<4} │ {2: >4}-{3: >4}-{4: <4} │ {5:>6.4f}-{6:<6.4f} │ {7: ^6.4f} | {8:^6.4e} | {9:^8} | {10:^8}",
                     s.nd[0]->locId, s.nd[1]->locId,
                     s.face->nds[0]->locId,s.face->nds[1]->locId,s.face->nds[2]->locId,
                     s.angle0, s.angle1, s.centerAngle, s.face->area_initial, s.face->isCWBoundary(this), s.face->isCCWBoundary(this));
}


void icy::Node::UpdateFan()
{
    // assumes that u and v are normalized
    auto get_angle = [](Eigen::Vector2d u, Eigen::Vector2d v)
    { return acos(std::clamp((double)u.dot(v),-1.0,1.0)); };

    fan_angle_span = 0;

    for(Sector &f : fan)
    {
        f.u_normalized = (f.nd[0]->xt2() - this->xt2()).normalized();
        f.v_normalized = (f.nd[1]->xt2() - this->xt2()).normalized();

        f.angle0 = fan_angle_span;
        f.angle1 = fan_angle_span += get_angle(f.u_normalized,f.v_normalized);

        f.u_p << -f.u_normalized[1], f.u_normalized[0];
        f.v_p << -f.v_normalized[1], f.v_normalized[0];

        f.t0_top << f.face->str_top * f.u_p;
        f.t1_top << f.face->str_top * f.v_p;
        f.t0_bottom << f.face->str_bottom * f.u_p;
        f.t1_bottom << f.face->str_bottom * f.v_p;
    }
}





void icy::Node::ComputeFanVariables(const SimParams &prms)
{
    double dont_split_nearly_degenerate_elems = prms.FractureAngleThreshold*M_PI/180;

    if(fan.size()==0) throw std::runtime_error("invoking ComputeFanVariables on a Node without elements");
    dir.setZero();
    max_normal_traction = 0;
    UpdateFan();
    unsigned nFan = fan.size();
    if(nFan==1 || fan_angle_span < dont_split_nearly_degenerate_elems) return;

    double weakening_coeff = prms.FractureWeakeningCoeff;
    unsigned gridPts = isBoundary ? nFan+1 : nFan;

    double grid_results[gridPts];
    for(unsigned i=0; i<nFan; i++)
    {
        grid_results[i] = NormalTraction(fan[i].angle0, weakening_coeff);
        if(std::isnan(grid_results[i])) throw std::runtime_error("ComputeFanVariables: traction is NaN");
    }
    if(isBoundary)
    {
        grid_results[nFan] = NormalTraction(fan[nFan-1].angle1, weakening_coeff);
        if(std::isnan(grid_results[nFan])) throw std::runtime_error("ComputeFanVariables: traction is NaN");
    }

    double *highest_grid_pt = std::max_element(grid_results, &grid_results[gridPts]);
    unsigned idx = std::distance(grid_results, highest_grid_pt);

    // reject if the grid max is low
    if(*highest_grid_pt < prms.FractureTractionThreshold*0.4) return;

    // sectors
    int sector1, sector2;

    if(isBoundary && (idx==0 || idx==gridPts-1))
    {
        sector1 = idx==0 ? 0 : gridPts-2;
        sector2 = -1;
    }
    else
    {
        sector1 = idx;
        sector2 = (idx-1+nFan)%nFan;
    }

    int bits = std::numeric_limits<float>::digits/2;

    boost::uintmax_t max_iter = 15;
    auto [fracture_angle, max1] = boost::math::tools::brent_find_minima(
        [&](double x){return -NormalTraction(x, weakening_coeff);},
        fan[sector1].angle0, fan[sector1].angle1, bits, max_iter);
    max_normal_traction = -max1;

    if(sector2 > -1)
    {
        max_iter = 15;
        auto [fracture_angle2, max2] = boost::math::tools::brent_find_minima(
            [&](double x){return -NormalTraction(x, weakening_coeff);},
            fan[sector2].angle0, fan[sector2].angle1, bits, max_iter);
        max2 = -max2;
        if(max2 > max_normal_traction) fracture_angle = fracture_angle2;
    }

    EvaluateTractions(fracture_angle, result_with_max_traction, weakening_coeff);

    if(result_with_max_traction.faces[0]==result_with_max_traction.faces[1] || result_with_max_traction.faces[0]==nullptr)
    {
        spdlog::critical("evaluate_tractions: face0=={}; face1=={}",
                         (void*)result_with_max_traction.faces[0],
                         (void*)result_with_max_traction.faces[1]);
        spdlog::critical("fracture_angle: {}; fan_angle_span {}",fracture_angle, fan_angle_span);
        PrintoutFan();
        throw std::runtime_error("evaluate_tractions: face0==face1");
    }

    if(!result_with_max_traction.faces[0]->containsNode(this))
    {
        spdlog::critical("ComputeFanVariables: mesh topology error 0");
        throw std::runtime_error("ComputeFanVariables: mesh topology error 0");
    }

    if(result_with_max_traction.faces[1] != nullptr && !result_with_max_traction.faces[1]->containsNode(this))
    {
        spdlog::critical("ComputeFanVariables: mesh topology error 1");
        throw std::runtime_error("ComputeFanVariables: mesh topology error 1");
    }

    max_normal_traction = result_with_max_traction.trac_normal_max;
    dir = result_with_max_traction.tn;

    // don't break nearly-degenerate sectors
    double span0 = result_with_max_traction.sectorSpan(0);
    if(result_with_max_traction.phi[0] < dont_split_nearly_degenerate_elems ||
        (result_with_max_traction.faces[0]->area_initial < prms.FractureAreaThreshold &&
         result_with_max_traction.phi[0] < result_with_max_traction.theta[0]))
    {
        result_with_max_traction.phi[0] = 0;
        result_with_max_traction.theta[0] = span0;
        fracture_angle = result_with_max_traction.angle0[0];
    }
    else if(result_with_max_traction.theta[0] < dont_split_nearly_degenerate_elems ||
             (result_with_max_traction.faces[0]->area_initial < prms.FractureAreaThreshold &&
              result_with_max_traction.phi[0] > result_with_max_traction.theta[0]))
    {
        result_with_max_traction.theta[0] = 0;
        result_with_max_traction.phi[0] = span0;
        fracture_angle = result_with_max_traction.angle1[0];
    }

    const double threshold_angle = dont_split_nearly_degenerate_elems;
    if(isBoundary && (fracture_angle < threshold_angle ||
                       fracture_angle > fan_angle_span-threshold_angle))
    {max_normal_traction=0; return;}

    if(!isBoundary)
    {
        double fracture_angle_bwd;
        // similar snap on the other side
        // TODO: there is a narrow case that may result in invalid topology
        double span1 = result_with_max_traction.sectorSpan(1);
        if(result_with_max_traction.phi[1] < dont_split_nearly_degenerate_elems ||
            (result_with_max_traction.faces[1]->area_initial < prms.FractureAreaThreshold &&
             result_with_max_traction.phi[1] < result_with_max_traction.theta[1]))
        {
            result_with_max_traction.phi[1] = 0;
            result_with_max_traction.theta[1] = span1;
            fracture_angle_bwd = result_with_max_traction.angle0[1];
            if(fracture_angle_bwd == fracture_angle) {max_normal_traction=0; return;}
        }
        else if(result_with_max_traction.theta[1] < dont_split_nearly_degenerate_elems ||
                 (result_with_max_traction.faces[1]->area_initial < prms.FractureAreaThreshold &&
                  result_with_max_traction.phi[1] > result_with_max_traction.theta[1]))
        {
            result_with_max_traction.theta[1] = 0;
            result_with_max_traction.phi[1] = span1;
            fracture_angle_bwd = result_with_max_traction.angle1[1];
            if(fracture_angle_bwd == fracture_angle) {max_normal_traction=0; return;}
        }
    }
}



double icy::Node::NormalTraction(double angle_fwd, double weakening_coeff) const
{
    SepStressResult tmpSsr;
    EvaluateTractions(angle_fwd, tmpSsr, weakening_coeff);
    return tmpSsr.trac_normal_max;
}


void icy::Node::EvaluateTractions(double angle_fwd, SepStressResult &ssr, const double weakening_coeff) const
{
    ssr.traction_top[0].setZero();
    ssr.traction_top[1].setZero();
    ssr.traction_bottom[0].setZero();
    ssr.traction_bottom[1].setZero();
    ssr.faces[0] = ssr.faces[1] = nullptr;

    if(angle_fwd == fan_angle_span) angle_fwd = std::max(0.0,angle_fwd-1e-11);
    ssr.angle_fwd = angle_fwd;

    double angle_bwd = angle_fwd+fan_angle_span/2;
    if (angle_bwd >= fan_angle_span) angle_bwd -= fan_angle_span;
    ssr.angle_bwd = angle_bwd;

    // integrate traction
    int sector = (isBoundary || angle_fwd < angle_bwd) ? 0 : 1;

    std::size_t nFans = fan.size();

    for (std::size_t f=0; f < nFans; f++)
    {
        const Sector &fp = fan[f];

        if (angle_fwd >= fp.angle0 && angle_fwd < fp.angle1)
        {
            ssr.faces[0] = fp.face;

            double phi = ssr.phi[0] = angle_fwd - fp.angle0;
            ssr.theta[0] = fp.angle1 - angle_fwd;
            ssr.angle0[0] = fp.angle0;
            ssr.angle1[0] = fp.angle1;

            double ratio = phi/(fp.angle1-fp.angle0);
            ssr.tn = (fp.u_normalized*(1-ratio) + fp.v_normalized*ratio).normalized();
            ssr.tn_p = (fp.u_p*(1-ratio) + fp.v_p*ratio).normalized(); // perpendicular to tn
            Eigen::Vector2d tmult_top = fp.face->str_top * ssr.tn_p;
            Eigen::Vector2d tmult_bottom = fp.face->str_bottom * ssr.tn_p;

            ssr.traction_top[sector] += tmult_top - fp.t0_top;
            ssr.traction_bottom[sector] += tmult_bottom - fp.t0_bottom;
            sector = 1-sector;
            ssr.traction_top[sector] += fp.t1_top - tmult_top;
            ssr.traction_bottom[sector] += fp.t1_bottom - tmult_bottom;
        }
        else if (!isBoundary && angle_bwd >= fp.angle0 && angle_bwd < fp.angle1)
        {
            ssr.faces[1] = fp.face;

            double phi = ssr.phi[1] = angle_bwd - fp.angle0;
            ssr.theta[1] = fp.angle1 - angle_bwd;
            ssr.angle0[1] = fp.angle0;
            ssr.angle1[1] = fp.angle1;

            double ratio = phi/(fp.angle1-fp.angle0);
            Eigen::Vector2d tn_p = (fp.u_p*(1-ratio) + fp.v_p*ratio).normalized(); // perpendicular to tn

            Eigen::Vector2d tmult_top = fp.face->str_top * ssr.tn_p;
            Eigen::Vector2d tmult_bottom = fp.face->str_bottom * ssr.tn_p;

            ssr.traction_top[sector] += tmult_top - fp.t0_top;
            ssr.traction_bottom[sector] += tmult_bottom - fp.t0_bottom;
            sector = 1-sector;
            ssr.traction_top[sector] += fp.t1_top - tmult_top;
            ssr.traction_bottom[sector] += fp.t1_bottom - tmult_bottom;
        }
        else
        {
            ssr.traction_top[sector] += fp.t1_top - fp.t0_top;
            ssr.traction_bottom[sector] += fp.t1_bottom - fp.t0_bottom;
        }
    }   // nFans

    double t0_tangential_top = ssr.traction_top[0].dot(ssr.tn);
    double t1_tangential_top = ssr.traction_top[1].dot(ssr.tn);
    double t0_normal_top = ssr.tn_p.dot(ssr.traction_top[0]);
    double t1_normal_top = -ssr.tn_p.dot(ssr.traction_top[1]);
    ssr.trac_normal_top = t0_normal_top + t1_normal_top;
    ssr.trac_tangential_top = t0_tangential_top - t1_tangential_top;

    double t0_tangential_bottom = ssr.traction_bottom[0].dot(ssr.tn);
    double t1_tangential_bottom = ssr.traction_bottom[1].dot(ssr.tn);
    double t0_normal_bottom = ssr.tn_p.dot(ssr.traction_bottom[0]);
    double t1_normal_bottom = -ssr.tn_p.dot(ssr.traction_bottom[1]);
    ssr.trac_normal_bottom = t0_normal_bottom + t1_normal_bottom;
    ssr.trac_tangential_bottom = t0_tangential_bottom - t1_tangential_bottom;


    if(!isBoundary)
    {
        ssr.trac_normal_top /= 2;
        ssr.trac_tangential_top /= 2;
        ssr.trac_normal_bottom /= 2;
        ssr.trac_tangential_bottom /= 2;
    }

    if(isCrackTip)
    {
        double coeff = ((1-weakening_coeff)+(weakening_coeff)*pow((weakening_direction.dot(ssr.tn)+1)/2, 5));
        ssr.trac_normal_top *= coeff;
        ssr.trac_normal_bottom *= coeff;
    }

    ssr.trac_normal_max = std::max(ssr.trac_normal_top, ssr.trac_normal_bottom);
}









void icy::Node::ReplaceAdjacentElement(Element *originalElem, Element *replacement)
{
    auto iter = std::find(adj_elems.begin(),adj_elems.end(),originalElem);
    if(iter == adj_elems.end())
        throw std::runtime_error("icy::Node::ReplaceAdjacentElement: can't find the original element to replace");
    else *iter = replacement;
}

uint64_t icy::Node::make_key(Node *nd0, Node *nd1)
{
    int nd0idx = nd0->locId;
    int nd1idx = nd1->locId;
    if(nd0idx > nd1idx) std::swap(nd0idx, nd1idx);
    return ((uint64_t)nd0idx << 32) | nd1idx;
}



















void icy::Node::ComputeElasticForce(LinearSystem &ls, SimParams &prms, double timeStep, double totalTime)
{
    if(lsId < 0) return;
    double beta = prms.NewmarkBeta;
    double alpha = prms.HHTalpha;
    double mass = area*prms.Thickness*prms.IceDensity;
    if(mass <= 0) throw std::runtime_error("zero nodal mass");

    Eigen::Matrix<double,LinearSystem::DOFS,1> F;
    Eigen::Matrix<double,LinearSystem::DOFS,LinearSystem::DOFS> dF;

    F = at*mass;
    //F(2) -= prms.gravity*mass;
    F(3) = 0;
    F(4) = 0;

    dF = Eigen::Matrix<double,LinearSystem::DOFS,LinearSystem::DOFS>::Identity()*(mass/(beta*timeStep*timeStep));
    dF(3,3) = 0;
    dF(4,4) = 0;


    // loading with surface waves

    vertical_force = 0;
    double spring = area*prms.WaterDensity*prms.gravity;
    double water_line = WaterLine(x_initial(0), x_initial(1), totalTime, prms);

    double disp_t = xt(2)-water_line;
    double disp_n = xn(2)-water_line;
    std::clamp(disp_t, -prms.Thickness*0.1, prms.Thickness*0.9);
    std::clamp(disp_n, -prms.Thickness*0.1, prms.Thickness*0.9);

    F(2) += disp_t*spring*(1-alpha);
    F(2) += disp_n*spring*alpha;
    dF(2,2) += spring*(1-alpha);
    vertical_force = (disp_t*spring*(1-alpha) + disp_n*spring*alpha)/area;

    // damping force
    double vert_velocity = WaterLineDt(x_initial(0), x_initial(1), totalTime, prms);
    double velocity_difference = vt.z()-vert_velocity;
    F(2) += prms.Damping*mass*(velocity_difference)/timeStep;
    dF(2,2) += prms.Damping*mass*prms.NewmarkGamma/(prms.NewmarkBeta*timeStep*timeStep);

    if(prms.loadType == icy::Model::LoadOpt::stretch_x)
    {
        // horizontal split in the middle

        double attenuation = totalTime < 1 ? totalTime : 1;
        double disp = x_initial.x() > 0 ? attenuation : -attenuation;
        disp*=prms.wave_height;
        double dispx_t = ut.x()-disp;
        double dispx_n = un.x()-disp;

        F(0) += dispx_t*spring*(1-alpha);
        F(0) += dispx_n*spring*alpha;

        dF(0,0) += spring*(1-alpha);

        F(1) += ut.y()*spring*(1-alpha);
        F(1) += un.y()*spring*alpha;
        dF(1,1) += spring*(1-alpha);
    }
    else if(prms.loadType == icy::Model::LoadOpt::stretch_xy)
    {

        if(totalTime < 3) spring+=200*spring;
//        else if(totalTime < 10) spring+=200*spring*((10-totalTime)/5);
        // radial stretch in all directions
        double attenuation10 = totalTime < 20 ? totalTime/20 : 1;
        Eigen::Vector2d vec(x_initial.x(), x_initial.y());
        vec*=(1+attenuation10*0.1);

        Eigen::Vector2d disp_t = xt.block(0,0,2,1)-vec;
        Eigen::Vector2d disp_n = xn.block(0,0,2,1)-vec;

        F(0) += disp_t.x()*spring*(1-alpha);
        F(0) += disp_n.x()*spring*alpha;
        F(1) += disp_t.y()*spring*(1-alpha);
        F(1) += disp_n.y()*spring*alpha;

        dF(0,0) += spring*(1-alpha);
        dF(1,1) += spring*(1-alpha);

        F(0) += prms.Damping*mass*(vt.x())/timeStep;
        dF(0,0) += prms.Damping*mass*prms.NewmarkGamma/(prms.NewmarkBeta*timeStep*timeStep);
        F(1) += prms.Damping*mass*(vt.y())/timeStep;
        dF(1,1) += prms.Damping*mass*prms.NewmarkGamma/(prms.NewmarkBeta*timeStep*timeStep);
    }
    else if(prms.loadType == icy::Model::LoadOpt::indentation)
    {
        // center indentation
        const double ind_radius = 1;
        const double ind_rate = 2.0/100;
        double rsq = x_initial.x() * x_initial.x() + x_initial.y() * x_initial.y();
        double r = sqrt(rsq);
        if(r < ind_radius)
        {
            double sphere_z = sqrt(ind_radius*ind_radius - rsq) - ind_radius + totalTime*ind_rate;
            if(sphere_z > 0)
            {
                double indented_position = -sphere_z;
                double disp_t = xt(2)-indented_position;
                double disp_n = xn(2)-indented_position;
                double spring2 = 0;//spring;//*100;
                if(disp_t>0) { spring2=spring*(1+disp_t*100);
                F(2) += disp_t*spring2*(1-alpha);
                F(2) += disp_n*spring2*alpha;
//                dF(2,2) += spring2*(1-alpha);
                dF(2,2) += spring*(1+disp_t*200)*(1-alpha);
                vertical_force += (disp_t)*spring2*(1-alpha) + (xn(2)-indented_position)*spring2*alpha;
                }
            }
        }
    }
    else if(prms.loadType == icy::Model::LoadOpt::waves_wind)
    {
        if(totalTime>5) {
            // add wind
            double attenuation10 = (totalTime-5) < 5 ? (totalTime-5)/5 : 1;
            double spring2 = 0.01*spring*attenuation10;

            Eigen::Vector2d vec(x_initial.x(), x_initial.y());
            vec.x()+=vec.x()*0.3*attenuation10*((x_initial.y()+50)/50)*((x_initial.y()+50)/50)*(1+abs(vec.x()/50));
            vec.y()+=((vec.y()+50)*(vec.y()+50)/2500)*50*attenuation10*0.5;

            Eigen::Vector2d disp_t = xt.block(0,0,2,1)-vec;
            Eigen::Vector2d disp_n = xn.block(0,0,2,1)-vec;

            F(0) += disp_t.x()*spring2*(1-alpha);
            F(0) += disp_n.x()*spring2*alpha;
            F(1) += disp_t.y()*spring2*(1-alpha);
            F(1) += disp_n.y()*spring2*alpha;

            dF(0,0) += spring2*(1-alpha);
            dF(1,1) += spring2*(1-alpha);

            F(0) += prms.Damping*mass*(vt.x())/timeStep;
            dF(0,0) += prms.Damping*mass*prms.NewmarkGamma/(prms.NewmarkBeta*timeStep*timeStep);
            F(1) += prms.Damping*mass*(vt.y())/timeStep;
            dF(1,1) += prms.Damping*mass*prms.NewmarkGamma/(prms.NewmarkBeta*timeStep*timeStep);
        }
    }

    for(int i=0;i<LinearSystem::DOFS;i++) if(std::isnan(dF(i))) throw std::runtime_error("node; dF contains NaN");

    // assemble
    ls.AddToEquation(F.data(),dF.data(),{lsId});
}

double icy::Node::Smoothstep(double edge0, double edge1, double x)
{
    x = (x-edge0)/(edge1-edge0);
    if(x>1.0) return 1;
    else if(x<0) return 0;
    else return x * x * (3 - 2 * x);
}

double icy::Node::SmoothstepDeriv(double edge0, double edge1, double x)
{
    x = (x-edge0)/(edge1-edge0);
    if(x>1.0 || x<0) return 0;
    else return 6*x*(1-x)/(edge1-edge0);
}

double icy::Node::WaterLine(double x, double y, double t, SimParams &prms)
{
    if(prms.loadType == icy::Model::LoadOpt::waterfall)
    {
        return -prms.wave_height*Smoothstep(0, 1.0, x+t-prms.wave_start_location);
    }
    else if(prms.loadType == icy::Model::LoadOpt::waves_x)
    {
        double result = prms.wave_height*cos(x*2*M_PI/3.99)*sin(t * 2 * M_PI / 1.6);
        if(t < 2) result *= t/2;
        return result;
    }
    else if(prms.loadType == icy::Model::LoadOpt::waves_xy)
    {
        double wave1 = prms.wave_height*cos(x*2*M_PI/3.99)*sin(t * 2 * M_PI / 1.6);
        double wave2 = prms.wave_height*0.8*cos(y*2*M_PI/5)*sin(2+t * 2 * M_PI / 1.2);
        double wave3 = prms.wave_height*0.5*cos(y*2*M_PI/10)*sin(1+t * 2 * M_PI / 2);
        if(t < 2) wave1 *= t/2;
        if(t < 4) wave2 *= t/4;
        if(t < 6) wave3 *= t/6;
        return wave1+wave2+wave3;
    }
    else if(prms.loadType == icy::Model::LoadOpt::waves_diag ||
            prms.loadType == icy::Model::LoadOpt::waves_wind)
    {
        Eigen::Vector2d dir1(1,1);
        Eigen::Vector2d dir2(1,-1);
        dir1.normalize();
        dir2.normalize();
        Eigen::Vector2d dir(x,y);
        double wavelength1=4;
        double wavelength2=6;
        double velocity1 = 1;
        double velocity2 = 1.3;
        double A = prms.wave_height;
        double wave1 = A*sin(M_PI*2*(dir1.dot(dir)-t*velocity1)/wavelength1);
        double wave2 = A*sin(M_PI*2*(dir2.dot(dir)-t*velocity2)/wavelength2);
        if(t < 2) wave1 *= t/2;
        if(t < 4) wave2 *= t/4;
        double total=wave1+wave2;
        double coeff = (y+50)*(y+50)/(50*50);
        if(coeff>1) coeff=1;
        total*=exp(-t/10)*coeff;
        return total;
    }
    else return 0;
}

double icy::Node::WaterLineDt(double x, double y, double t, SimParams &prms)
{
    if(prms.loadType == icy::Model::LoadOpt::waterfall)
    {
        return -prms.wave_height*SmoothstepDeriv(0, 1.0, x+t-prms.wave_start_location);
    }
    else if(prms.loadType == icy::Model::LoadOpt::waves_x)
    {
        double result = prms.wave_height*cos(x*2*M_PI/3.99)*cos(t * 2 * M_PI / 1.6)* 2 * M_PI / 1.6;
        if(t < 2) result *= t/2;
        return result;
    }
    else if(prms.loadType == icy::Model::LoadOpt::waves_xy)
    {
        double wave1 = prms.wave_height*cos(x*2*M_PI/3.99)*cos(t * 2 * M_PI / 1.6)* 2 * M_PI / 1.6;
        double wave2 = prms.wave_height*0.5*cos(y*2*M_PI/5)*cos(2+t * 2 * M_PI / 1.2)* 2 * M_PI / 1.2;
        double wave3 = prms.wave_height*0.5*cos(y*2*M_PI/10)*cos(1+t * 2 * M_PI / 2)* 2 * M_PI / 2;
        if(t < 2) wave1 *= t/2;
        if(t < 4) wave2 *= t/4;
        if(t < 6) wave3 *= t/6;
        return wave1+wave2+wave3;
    }
    else if(prms.loadType == icy::Model::LoadOpt::waves_diag ||
            prms.loadType == icy::Model::LoadOpt::waves_wind)
    {
        Eigen::Vector2d dir1(1,1);
        Eigen::Vector2d dir2(1,-1);
        dir1.normalize();
        dir2.normalize();
        Eigen::Vector2d dir(x,y);
        double wavelength1=4;
        double wavelength2=6;
        double velocity1 = 1;
        double velocity2 = 1.3;
        double A = prms.wave_height;
        double wave1 = A*cos(M_PI*2*(dir1.dot(dir)-t*velocity1)/wavelength1)*M_PI*2*(-velocity1)/wavelength1;
        double wave2 = A*cos(M_PI*2*(dir2.dot(dir)-t*velocity2)/wavelength2)*M_PI*2*(-velocity2)/wavelength2;
        if(t < 2) wave1 *= t/2;
        if(t < 4) wave2 *= t/4;
        double total=wave1+wave2;
        double coeff = (y+50)*(y+50)/(50*50);
        if(coeff>1) coeff=1;
        total*=exp(-t/10)*coeff;
        return total;
    }
    else return 0;
}

double icy::Node::BellShapedPolynomial(double x)
{
    if(x<0) x=-x;
    if(x>2) return 0;
    if(x<1) return 0.25*(4-6*x*x+3*x*x*x);
    return 0.25*(2-x)*(2-x)*(2-x);
}

double icy::Node::BellShapedPolynomialDx(double x)
{
    const double k=(3.0/4.0);
    if(x>=1 && x<2) return -k*(x-2)*(x-2);
    else if(x>-2 && x<=-1) return k*(2+x)*(2+x);
    else if(x>=0 && x<1) return k*x*(3*x-4);
    else if(x>-1 && x <0) return -k*x*(4+3*x);
    else return 0;
}




void icy::Node::AcceptTentativeValues()
{
    un = ut;
    xn = xt;
    vn = vt;
    an = at;
}




