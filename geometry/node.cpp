#include <cmath>
#include <cfloat>
#include <algorithm>
#include "node.h"
#include "model.h"
#include "parameters_sim.h"
#include "element.h"
#include "boost/math/tools/minima.hpp"

icy::Node::Node() { Reset(); }

void icy::Node::ComputeElasticForce(LinearSystem &ls, SimParams &prms, double timeStep, double totalTime)
{
    if(lsId < 0) return;
    double beta = prms.NewmarkBeta;
    double alpha = prms.HHTalpha;
    double mass = area*prms.Thickness*prms.IceDensity;
    if(mass <= 0) throw std::runtime_error("zero nodal mass");

    Eigen::Matrix<double,DOFS,1> F;
    Eigen::Matrix<double,DOFS,DOFS> dF;

    F = Eigen::Matrix<double,DOFS,1>::Zero();
    F = at*mass;
    //F(2) -= prms.gravity*mass;
    F(3) = 0;
    F(4) = 0;

    dF = Eigen::Matrix<double,DOFS,DOFS>::Identity()*(mass/(beta*timeStep*timeStep));
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

#ifdef QT_DEBUG
    for(int i=0;i<DOFS;i++) if(std::isnan(dF(i))) throw std::runtime_error("node; dF contains NaN");
#endif

    // assemble
    ls.SubtractRHS(lsId, F);
    ls.AddLHS(lsId, lsId, dF);
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

// ======================================================

void icy::Node::Reset()
{
    ut=xt=vt=at=un=xn=vn=an=Eigen::Matrix<double,DOFS,1>::Zero();
    x_initial=Eigen::Matrix<double,3,1>::Zero();
    normal_n = Eigen::Vector3d::Zero();
    area = 0;
    adjacent_elems.reserve(7);
    adjacent_elems.clear();
    crack_tip = support_node = reset_timing = false;
    timeLoadedAboveThreshold = 0;
    max_normal_traction = 0;
    dir=weakening_direction = Eigen::Vector2f::Zero();
}

void icy::Node::InitializeFromAdjacent(const Node *nd0, const Node *nd1, double f)
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

void icy::Node::InitializeFromAnother(const Node *nd)
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


void icy::Node::AcceptTentativeValues()
{
    un = ut;
    xn = xt;
    vn = vt;
    an = at;
}

void icy::Node::InitializeFan()
{
    auto get_angle = [](Eigen::Vector2f u, Eigen::Vector2f v)
    {
        double dot = u.dot(v);
//        double dot = u.dot(v)/(u.norm()*v.norm());
        if(dot > 1) dot = 1.0;
        else if(dot < -1.0) dot = -1.0;
        return acos(dot);
    };

    fan_angle_span = 0;

    for(Sector &f : fan)
    {
        // TODO: Simplify these two statements
        f.u_normalized = f.face->CWEdge(this).getVec(this).normalized();
        f.v_normalized = f.face->CCWEdge(this).getVec(this).normalized();

        f.angle0 = fan_angle_span;
        fan_angle_span += get_angle(f.u_normalized,f.v_normalized);
        f.angle1 = fan_angle_span;

        f.u_p << -f.u_normalized.y(), f.u_normalized.x();
        f.v_p << -f.v_normalized.y(), f.v_normalized.x();

        f.t0_top << f.face->str_top * f.u_p;
        f.t1_top << f.face->str_top * f.v_p;
        f.t0_bottom << f.face->str_bottom * f.u_p;
        f.t1_bottom << f.face->str_bottom * f.v_p;
    }
}


void icy::Node::evaluate_tractions(float angle_fwd, SepStressResult &ssr, const float weakening_coeff) const
{
    ssr.traction_top[0] = ssr.traction_top[1] = Eigen::Vector2f::Zero();
    ssr.traction_bottom[0] = ssr.traction_bottom[1] = Eigen::Vector2f::Zero();
    ssr.faces[0] = ssr.faces[1] = nullptr;

    if(angle_fwd == fan_angle_span) angle_fwd -= 1e-4;
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
            ssr.e[0] = fp.face->CWEdge(this);
            ssr.e[1] = fp.face->CCWEdge(this);
            ssr.e_opposite[0] = fp.face->OppositeEdge(this);

            float phi = ssr.phi[0] = angle_fwd - fp.angle0;
            ssr.theta[0] = fp.angle1 - angle_fwd;

            float ratio = phi/(fp.angle1-fp.angle0);
            ssr.tn = (fp.u_normalized*(1-ratio) + fp.v_normalized*ratio).normalized();
            ssr.tn_p = (fp.u_p*(1-ratio) + fp.v_p*ratio).normalized(); // perpendicular to tn
            //ssr.tn_p = normal_n.cross(ssr.tn).normalized();
            Eigen::Vector2f tmult_top = fp.face->str_top * ssr.tn_p;
            Eigen::Vector2f tmult_bottom = fp.face->str_bottom * ssr.tn_p;

            ssr.traction_top[sector] += tmult_top - fp.t0_top;
            ssr.traction_bottom[sector] += tmult_bottom - fp.t0_bottom;
            sector = 1-sector;
            ssr.traction_top[sector] += fp.t1_top - tmult_top;
            ssr.traction_bottom[sector] += fp.t1_bottom - tmult_bottom;
        }
        else if (!isBoundary && angle_bwd >= fp.angle0 && angle_bwd < fp.angle1)
        {
            ssr.faces[1] = fp.face;
            ssr.e[2] = fp.face->CWEdge(this);
            ssr.e[3] = fp.face->CCWEdge(this);
            ssr.e_opposite[1] = fp.face->OppositeEdge(this);

            float phi = ssr.phi[1] = angle_bwd - fp.angle0;
            ssr.theta[1] = fp.angle1 - angle_bwd;

            float ratio = phi/(fp.angle1-fp.angle0);
            Eigen::Vector2f tn_p = (fp.u_p*(1-ratio) + fp.v_p*ratio).normalized(); // perpendicular to tn

            Eigen::Vector2f tmult_top = fp.face->str_top * tn_p;
            Eigen::Vector2f tmult_bottom = fp.face->str_bottom * tn_p;

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

    float t0_tangential_top = ssr.traction_top[0].dot(ssr.tn);
    float t1_tangential_top = ssr.traction_top[1].dot(ssr.tn);
    float t0_normal_top = ssr.tn_p.dot(ssr.traction_top[0]);
    float t1_normal_top = -ssr.tn_p.dot(ssr.traction_top[1]);
    ssr.trac_normal_top = t0_normal_top + t1_normal_top;
    ssr.trac_tangential_top = t0_tangential_top - t1_tangential_top;

    float t0_tangential_bottom = ssr.traction_bottom[0].dot(ssr.tn);
    float t1_tangential_bottom = ssr.traction_bottom[1].dot(ssr.tn);
    float t0_normal_bottom = ssr.tn_p.dot(ssr.traction_bottom[0]);
    float t1_normal_bottom = -ssr.tn_p.dot(ssr.traction_bottom[1]);
    ssr.trac_normal_bottom = t0_normal_bottom + t1_normal_bottom;
    ssr.trac_tangential_bottom = t0_tangential_bottom - t1_tangential_bottom;

    if(!isBoundary)
    {
        ssr.trac_normal_bottom /= 2;
        ssr.trac_tangential_bottom /= 2;
        ssr.trac_normal_top /= 2;
        ssr.trac_tangential_top /= 2;
    }

    if(crack_tip)
    {
        // TODO: replace pow(...) with something simpler to compute
        float coeff = ((1-weakening_coeff)+(weakening_coeff)*pow((weakening_direction.dot(ssr.tn)+1)/2, 5));
        ssr.trac_normal_bottom*=coeff;
        ssr.trac_normal_top*=coeff;
    }

    ssr.trac_normal_max = std::max(ssr.trac_normal_top, ssr.trac_normal_bottom);
}

float icy::Node::normal_traction(float angle_fwd, float weakening_coeff) const
{
    SepStressResult tmpSsr;
    evaluate_tractions(angle_fwd, tmpSsr, weakening_coeff);
    return tmpSsr.trac_normal_max;
}

void icy::Node::ComputeFanVariablesAlt(SimParams &prms)
{
    InitializeFan();
    dir=Eigen::Vector2f::Zero();
    max_normal_traction = 0;
    unsigned nFan = fan.size();

    double weakening_coeff = prms.weakening_coeff;

    unsigned gridPts = isBoundary ? nFan+1 : nFan;

    float grid_results[gridPts];
    for(unsigned i=0; i<nFan; i++)
    {
        grid_results[i] = normal_traction(fan[i].angle0, weakening_coeff);
        if(std::isnan(grid_results[i])) throw std::runtime_error("traction is nan");
    }
    if(isBoundary) {
        grid_results[nFan] = normal_traction(fan[nFan-1].angle1, weakening_coeff);
        if(std::isnan(grid_results[nFan])) throw std::runtime_error("traction is nan");
    }

    float *highest_grid_pt = std::max_element(grid_results, &grid_results[gridPts]);
    unsigned idx = std::distance(grid_results, highest_grid_pt);

    // reject if the grid max is low
    if(*highest_grid_pt < prms.normal_traction_threshold*prms.cutoff_coefficient) return;

    // sectors
    int sector1, sector2;

    if(isBoundary && (idx == 0 || idx==gridPts-1))
    {
        sector1 = idx == 0 ? 0 : gridPts-2;
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
                    [=](double x){return -normal_traction(x, weakening_coeff);},
        fan[sector1].angle0, fan[sector1].angle1, bits, max_iter);
    max_normal_traction = -max1;

    if(sector2 > -1)
    {
        max_iter = 15;
        auto [fracture_angle2, max2] = boost::math::tools::brent_find_minima(
                        [=](double x){return -normal_traction(x, weakening_coeff);},
            fan[sector2].angle0, fan[sector2].angle1, bits, max_iter);
        max2 = -max2;
        if(max2 > max_normal_traction) fracture_angle = fracture_angle2;
    }

    evaluate_tractions(fracture_angle, result_with_max_traction, weakening_coeff);
    if(result_with_max_traction.faces[0]==result_with_max_traction.faces[1])
        throw std::runtime_error("evaluate_tractions: face0==face1");
    if(!result_with_max_traction.faces[0]->ContainsNode(this))
        throw std::runtime_error("ComputeFanVariablesAlt: mesh topology error 0");
    if(result_with_max_traction.faces[1]!= nullptr && !result_with_max_traction.faces[1]->ContainsNode(this))
        throw std::runtime_error("ComputeFanVariablesAlt: mesh topology error 1");
    max_normal_traction = result_with_max_traction.trac_normal_max;
    dir = result_with_max_traction.tn;

    const float threshold_angle = fan_angle_span*0.1;
    if(isBoundary && (fracture_angle < threshold_angle ||
                      fracture_angle > fan_angle_span-threshold_angle || fan_angle_span < M_PI/2))
    {max_normal_traction=0; return;}
}


void icy::Node::PrepareFan2()
{
    Eigen::Vector3d nd_vec = x_initial;

    unsigned nElems = adjacent_elems.size();
    if(nElems == 0)
    {
        qDebug() << "node " << this->locId;
        throw std::runtime_error("disconnected node");
    }
    fan.clear();
    fan.reserve(nElems);
    area = 0;
    for(unsigned k=0;k<nElems;k++)
    {
        icy::Element *elem = adjacent_elems[k];
        if(!elem->ContainsNode(this))  throw std::runtime_error("PrepareFan2: mesh topology error 0");
        area += elem->area_initial/3;

        Sector s;
        s.face = elem;
        Eigen::Vector3d tcv = elem->getCenter() - nd_vec;
        s.centerAngle = atan2(tcv.y(), tcv.x());

        short thisIdx, CWIdx, CCWIdx;
        elem->getIdxs(this, thisIdx, CWIdx, CCWIdx);

        s.nd[0] = elem->nds[CWIdx];
        s.nd[1] = elem->nds[CCWIdx];

        // note that the indices are swapped
        Edge e0 = elem->CWEdge(this);
        Edge e1 = elem->CCWEdge(this);
        fan.push_back(s);

        if(e0.isBoundary || e1.isBoundary) isBoundary = true;
    }

    std::sort(fan.begin(), fan.end(),
              [](const Sector &f0, const Sector &f1)
    {return f0.centerAngle < f1.centerAngle; });

    if(isBoundary) // assert means that PrepareFan2 is not called from Fix_X
    {
        // find the fan element with the border on the CW direction
        auto cw_boundary = std::find_if(fan.begin(), fan.end(), [this](const Sector &f){return f.face->CWEdge(this).isBoundary;});
        if(cw_boundary == fan.end())
        {
            PrintoutFan();
            throw std::runtime_error("cw boundary not found");
        }
        else
        {
        std::rotate(fan.begin(), cw_boundary, fan.end());
        }
    }

    // assert that the nodes of the fan connect
    for(std::size_t i = 0;i<fan.size()-1;i++)
    {
        if(fan[i].nd[1] != fan[i+1].nd[0])
        {
            std::cout << "\n\n\nfan nodes are not contiguous " << locId << std::endl;
            PrintoutFan();
            throw std::runtime_error("fan nodes are not contiguous");
        }
//        if(fan[i].e[1].nds[0] != fan[i+1].e[0].nds[0] || fan[i].e[1].nds[1] != fan[i+1].e[0].nds[1])
//            throw std::runtime_error("edges not shared");
    }
}

void icy::Node::PrintoutFan()
{
    std::cout << "Printing fan for node " << locId << (crack_tip ? " crack_tip" : " ") << std::endl;
    std::cout << "fan size " << fan.size() << "; isBoundary " << isBoundary << std::endl;
    std::cout << "adj elems size " << adjacent_elems.size() << std::endl;

    for(Sector &s : fan)
    {
        std::cout << s.nd[0]->locId << "-" << s.nd[1]->locId;
        std::cout << " ; " << s.face->nds[0]->locId << "-" << s.face->nds[1]->locId << "-"<< s.face->nds[2]->locId;
//        std::cout << " ; C " << s.e[0].nds[0]->locId << "-" << s.e[0].nds[1]->locId << (s.e[0].isBoundary ? " b " : " nb ");
//        std::cout << (s.e[0].toSplit ? "* " : " ");
        //std::cout << s.e[0].elems[0] << " " << s.e[0].elems[1];
//        std::cout << "; CC " << s.e[1].nds[0]->locId << "-" << s.e[1].nds[1]->locId << (s.e[1].isBoundary ? " b " : " nb ");
//        std::cout << (s.e[1].toSplit ? "* " : " ");
        //std::cout << s.e[1].elems[0] << " " << s.e[1].elems[1];
        std::cout << std::endl;
    }
    std::cout << "--------------------------------\n";
    std::cout << std::endl;
    for(int i=0;i<100;i++)
    std::cout << std::flush;
}


