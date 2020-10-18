#include <cmath>
#include <cfloat>
#include <algorithm>
#include "node.h"
#include "model.h"
#include "parameters_sim.h"
#include "element.h"
#include "boost/math/tools/minima.hpp"

icy::Node::Node()
{
    Reset();
    adjacent_nodes.reserve(10);
    fan.reserve(10);
}

void icy::Node::ComputeElasticForce(SimParams &prms, double timeStep, double totalTime)
{
    if(prescribed || lsId < 0) return;
    double beta = prms.NewmarkBeta;
    double alpha = prms.HHTalpha;
    double mass = area*prms.Thickness*prms.IceDensity;
    if(mass <= 0) throw std::runtime_error("zero nodal mass");

    F = Eigen::Matrix<double,DOFS,1>::Zero();
    F = at;
    F(2) -= prms.gravity*mass;
    F(3) = 0;
    F(4) = 0;

    dF = Eigen::Matrix<double,DOFS,DOFS>::Identity()*(mass/(beta*timeStep*timeStep));
    dF(3,3) = 0;
    dF(4,4) = 0;

    if(prms.loadType < 10)
    {
        // loading with surface waves

        vertical_force = 0;
        double spring = area*prms.WaterDensity*std::abs(prms.gravity);
        double water_line = WaterLine(x_initial(0), x_initial(1), totalTime, prms);

        double disp_t = xt(2)-water_line;
        double disp_n = xn(2)-water_line;
        std::clamp(disp_t, -prms.Thickness*2, prms.Thickness*2);
        std::clamp(disp_n, -prms.Thickness*2, prms.Thickness*2);

        F(2) += disp_t*spring*(1-alpha);
        F(2) += disp_n*spring*alpha;
        dF(2,2) += spring*(1-alpha);
        vertical_force = disp_t*spring*(1-alpha) + disp_n*spring*alpha;

        // damping force
        double vert_velocity = WaterLineDt(x_initial(0), x_initial(1), totalTime, prms);
        //    double max_velocity_difference = 3;
        double velocity_difference = vt.z()-vert_velocity;
        //    if(velocity_difference > max_velocity_difference) velocity_difference = max_velocity_difference;
        //    else if(velocity_difference < -max_velocity_difference) velocity_difference = -max_velocity_difference;

        F(2) += prms.Damping*mass*(velocity_difference)/timeStep;
        dF(2,2) += prms.Damping*mass*prms.NewmarkGamma/(prms.NewmarkBeta*timeStep*timeStep);
    }
    else
    {
        // other types of loading

    }



/*
    else if(prms.loadType == 2)
    {
        // center indentation
        const double ind_radius = 0.50;
        const double ind_rate = 1.0/120;
        double spring = area*prms.WaterDensity*std::abs(prms.gravity);
        double rsq = x_initial.x() * x_initial.x() + x_initial.y() * x_initial.y();
        double r = sqrt(rsq);
        if(r < ind_radius)
        {
            double sphere_z = sqrt(ind_radius*ind_radius - rsq) - ind_radius + totalTime*ind_rate;
            if(sphere_z > 0)
            {
                double indented_position = -sphere_z;
                double spring2 = spring*100;
                F(2) += (xt(2)-indented_position)*spring2*(1-alpha);
                F(2) += (xn(2)-indented_position)*spring2*alpha;
                dF(2,2) += spring2*(1-alpha);
                vertical_force = (xt(2)-indented_position)*spring2*(1-alpha) + (xn(2)-indented_position)*spring2*alpha;
            }
        }
        // add normal buoyancy "spring"
        F(2) += xt(2)*spring*(1-alpha);
        F(2) += xn(2)*spring*alpha;
        dF(2,2) += spring*(1-alpha);
    }
    */

    // assert
//    for(int i=0;i<DOFS;i++) if(std::isnan(dF(i))) throw std::runtime_error("node; dF contains NaN");

}

double icy::Node::OceanWave(double x, double y, double t)
{
    const int len = 4;
    double amp[] = {0.2, 0.0, 0.0, 0.5};
    double dir[][2] = {{1, 0},
                        {0.0995037, 0.995037},
                        {0.980581, 0.196116},
                        {0., 1}};
    double sp[] = {1, 3, 2, 3};
    double wavelen[] = {25.5, 10, 7, 38};
    double phase[] = {1, 0.1, 0.5, 0.1};

    amp[3] = 0;
    if(t<10) amp[0] *= t/10;
    else if(t>60) amp[0] *= exp(-(t-60)/10);


    double result = 0;
    for(int k=0;k<len;k++)
        result += amp[k]*sin(phase[k]+2*M_PI*((dir[k][0]*x+dir[k][1]*y-sp[k]*t)/wavelen[k]));
    return result;
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
    if(prms.loadType == 5)
    {
        return -prms.wave_height*Smoothstep(0, 1.0, x+t-25.3);
    }
    else if(prms.loadType == 6)
    {
        double result = prms.wave_height*cos(x*2*M_PI/3.99)*sin(t * 2 * M_PI / 1.6);
        if(t < 2) result *= t/2;
        return result;
    }
    else if(prms.loadType == 7)
    {
        double wave1 = prms.wave_height*cos(x*2*M_PI/3.99)*sin(t * 2 * M_PI / 1.6);
        double wave2 = prms.wave_height*0.1*cos(y*2*M_PI/5)*sin(2+t * 2 * M_PI / 1.2);
        double wave3 = prms.wave_height*0.1*cos(y*2*M_PI/10)*sin(1+t * 2 * M_PI / 2);
        if(t < 2) wave1 *= t/2;
        if(t < 4) wave2 *= t/4;
        if(t < 6) wave3 *= t/6;
        return wave1+wave2+wave3;
    }
    else return 0;
}

double icy::Node::WaterLineDt(double x, double y, double t, SimParams &prms)
{
    if(prms.loadType == 5)
    {
        return -prms.wave_height*SmoothstepDeriv(0, 1.0, x+t-25.3);
    }
    else if(prms.loadType == 6)
    {
        double result = prms.wave_height*cos(x*2*M_PI/3.99)*cos(t * 2 * M_PI / 1.6)* 2 * M_PI / 1.6;
        if(t < 2) result *= t/2;
        return result;
    }
    else if(prms.loadType == 7)
    {
        double wave1 = prms.wave_height*cos(x*2*M_PI/3.99)*cos(t * 2 * M_PI / 1.6)* 2 * M_PI / 1.6;
        double wave2 = prms.wave_height*0.1*cos(y*2*M_PI/5)*cos(2+t * 2 * M_PI / 1.2)* 2 * M_PI / 1.2;
        double wave3 = prms.wave_height*0.1*cos(y*2*M_PI/10)*cos(1+t * 2 * M_PI / 2)* 2 * M_PI / 2;
        if(t < 2) wave1 *= t/2;
        if(t < 4) wave2 *= t/4;
        if(t < 6) wave3 *= t/6;
        return wave1+wave2+wave3;
    }
    else return 0;
}


void icy::Node::Reset()
{
    x_initial=ut=xt=vt=at=un=xn=vn=an=Eigen::Matrix<double,DOFS,1>::Zero();
    weakening_direction = Eigen::Vector2f::Zero();
    normal_n = Eigen::Vector3d::Zero();
    prescribed = false;
    area = 0;
    adjacent_nodes.clear();
    //adjacent_edges_map.clear();
    crack_tip = support_node = core_node = false;
    timeLoadedAboveThreshold = 0;
}

void icy::Node::InitializeFromAdjacent(icy::Node *nd0, icy::Node *nd1, double f)
{
    x_initial = nd0->x_initial*f + nd1->x_initial*(1-f);
    ut = nd0->ut*f + nd1->ut*(1-f);
    un = nd0->un*f + nd1->un*(1-f);
    xt = x_initial + ut;
    xn = x_initial + un;
    vt = nd0->vt*f + nd1->vt*(1-f);
    vn = nd0->vn*f + nd1->vn*(1-f);
    at = nd0->at*f + nd1->at*(1-f);
    an = nd0->an*f + nd1->an*(1-f);
}

void icy::Node::InitializeFromAnother(icy::Node *nd)
{
    x_initial = nd->x_initial;
    ut = nd->ut;
    un = nd->un;
    xt = x_initial + ut;
    xn = x_initial + un;
    vt = nd->vt;
    vn = nd->vn;
    at = nd->at;
    an = nd->an;
}

void icy::Node::Assemble(icy::LinearSystem &ls) const
{
    if(prescribed || lsId < 0) return;
    ls.SubtractRHS(lsId, F);
    ls.AddLHS(lsId, lsId, dF);
}

void icy::Node::AcceptTentativeValues()
{
//    if(prescribed) return;
    un = ut;
    xn = xt;
    vn = vt;
    an = at;
}
/*
void icy::Node::PrepareFan(icy::Geometry *geom)
{
    std::sort(fan.begin(), fan.end(),
              [](const Sector &f0, const Sector &f1)
    {return f0.centerAngle < f1.centerAngle; });

    // TODO: find a more efficient way to infer these edges (from the element ?)
    for(Sector &f : fan)
    {
//        f.e[0] = adjacent_edges_map.at(f.nd[0]->locId);
//        f.e[1] = adjacent_edges_map.at(f.nd[1]->locId);
        f.e[0] = geom->getEdgeByNodalIdx(this->locId, f.nd[0]->locId);
        f.e[1] = geom->getEdgeByNodalIdx(this->locId, f.nd[1]->locId);
    }

    if(isBoundary)
    {
        // find the fan element with the border on the CW direction
        auto cw_boundary = std::find_if(fan.begin(), fan.end(), [](const Sector &f){return f.e[0].isBoundary;});
        if(cw_boundary == fan.end()) throw std::runtime_error("cw boundary not found");
        std::rotate(fan.begin(), cw_boundary, fan.end());
    }

    // assert that the nodes of the fan connect
    for(std::size_t i = 0;i<fan.size()-1;i++)
    {
        if(fan[i].nd[1] != fan[i+1].nd[0])
        {
            qDebug() << "nd: " << locId << "; fan nodes not contiguous";
            qDebug() << "fan size " << fan.size() << "; isBoundary " << isBoundary;
            qDebug() << "between idx " << i << " and " << i+1;
            qDebug() << "nodes ids: " << fan[i].nd[1]->locId << " and " << fan[i+1].nd[0]->locId;
            qDebug() << "----";
            for(std::size_t j = 0;j<fan.size();j++)
                qDebug() << "j="<<j<<"; nd_ids: " << fan[j].nd[0]->locId << " -- " << fan[j].nd[1]->locId;
            throw std::runtime_error("fan nodes are not contiguous");
        }
        if(fan[i].e[1].nds[0] != fan[i+1].e[0].nds[0] || fan[i].e[1].nds[1] != fan[i+1].e[0].nds[1])
            throw std::runtime_error("edges not shared");
    }
}
*/
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

    // TODO: make sure that normal_n is already computed at this stage
//    normal_n = Eigen::Vector3d::Zero();
//    for(Sector &f : fan) normal_n += f.face->normal_n;
//    normal_n.normalize();

    fan_angle_span = 0;

    for(Sector &f : fan)
    {
        f.u_normalized = f.e[0].getVec(this).normalized();
        f.v_normalized = f.e[1].getVec(this).normalized();

        f.angle0 = fan_angle_span;
        f.angle_span = get_angle(f.u_normalized,f.v_normalized);
        fan_angle_span += f.angle_span;
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

    if(angle_fwd == fan_angle_span) angle_fwd -= 1e-10;
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
            ssr.e[0] = fp.e[0];
            ssr.e[1] = fp.e[1];

            float phi = ssr.phi[0] = angle_fwd - fp.angle0;
            ssr.theta[0] = fp.angle1 - angle_fwd;

            float ratio = phi/fp.angle_span;
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
            ssr.e[2] = fp.e[0];
            ssr.e[3] = fp.e[1];

            float phi = ssr.phi[1] = angle_bwd - fp.angle0;
            ssr.theta[1] = fp.angle1 - angle_bwd;

            float ratio = phi/fp.angle_span;
            //Eigen::Vector2f tn_bwd = fp.u_normalized*(1-ratio) + fp.v_normalized*ratio;
            Eigen::Vector2f tn_p = (fp.u_p*(1-ratio) + fp.v_p*ratio).normalized(); // perpendicular to tn

//            Eigen::Vector3f tn_p = normal_n.cross(tn_bwd).normalized();
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
        // TODO: replace pow(...) with a something simpler to compute
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
    max_normal_traction = 0;
    unsigned nFan = fan.size();

    double weakening_coeff = prms.weakening_coeff;

    unsigned gridPts = isBoundary ? nFan+1 : nFan;

    float grid_results[gridPts];
    for(unsigned i=0; i<nFan; i++) grid_results[i] = normal_traction(fan[i].angle0, weakening_coeff);
    if(isBoundary) grid_results[nFan] = normal_traction(fan[nFan-1].angle1, weakening_coeff);

    float *highest_grid_pt = std::max_element(grid_results, &grid_results[gridPts]);
    unsigned idx = std::distance(grid_results, highest_grid_pt);

    // reject if the grid max is low
    if(*highest_grid_pt < prms.normal_traction_threshold/3) return;

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

    int bits = std::numeric_limits<float>::digits;

    boost::uintmax_t max_iter = 20;
    auto [fracture_angle, max1] = boost::math::tools::brent_find_minima(
                    [=](double x){return -normal_traction(x, weakening_coeff);},
        fan[sector1].angle0, fan[sector1].angle1, bits, max_iter);
    max_normal_traction = -max1;

    if(sector2 > -1)
    {
        max_iter = 20;
        auto [fracture_angle2, max2] = boost::math::tools::brent_find_minima(
                        [=](double x){return -normal_traction(x, weakening_coeff);},
            fan[sector2].angle0, fan[sector2].angle1, bits, max_iter);
        max2 = -max2;
        if(max2 > max_normal_traction) fracture_angle = fracture_angle2;
    }

    evaluate_tractions(fracture_angle, result_with_max_traction, weakening_coeff);
    max_normal_traction = result_with_max_traction.trac_normal_max;
    dir = result_with_max_traction.tn;

    const float threshold_angle = fan_angle_span*0.1;
    if(isBoundary && (fracture_angle < threshold_angle ||
                      fracture_angle > fan_angle_span-threshold_angle || fan_angle_span < M_PI/2))
    {max_normal_traction=0; return;}
}
/*
icy::Node::Sector::Sector(icy::Element *elem, icy::Node *ndd) : face(elem)
{
    // (the list of edges has not been built yet)
    Eigen::Vector3d nd_vec = ndd->x_initial.block(0,0,3,1);
    Eigen::Vector3d tcv = face->getCenter() - nd_vec;
    centerAngle = atan2(tcv.y(), tcv.x());

    nd[0] = face->getCWNode(ndd);
    nd[1] = face->getCCWNode(ndd);
}
*/

void icy::Node::PrepareFan2()
{
    Eigen::Vector3d nd_vec = x_initial.block(0,0,3,1);

    unsigned nElems = adjacent_elems.size();
    if(nElems == 0) throw std::runtime_error("disconnected node");
    fan.clear();
    for(unsigned k=0;k<nElems;k++)
    {
        icy::Element *elem = adjacent_elems[k];
        area += elem->area_initial/3;

        Sector s;
        s.face = elem;
        Eigen::Vector3d tcv = elem->getCenter() - nd_vec;
        s.centerAngle = atan2(tcv.y(), tcv.x());

        short CWIdx = elem->getCWIdx(this);
        short CCWIdx = elem->getCCWIdx(this);

        s.nd[0] = elem->nds[CWIdx];
        s.nd[1] = elem->nds[CCWIdx];

        // note that the indices are swapped
        s.e[0] = elem->edges[CCWIdx];
        s.e[1] = elem->edges[CWIdx];
        fan.push_back(s);

        if(s.e[0].isBoundary || s.e[1].isBoundary) isBoundary = true;
    }

    std::sort(fan.begin(), fan.end(),
              [](const Sector &f0, const Sector &f1)
    {return f0.centerAngle < f1.centerAngle; });

    if(isBoundary)
    {
        // find the fan element with the border on the CW direction
        auto cw_boundary = std::find_if(fan.begin(), fan.end(), [](const Sector &f){return f.e[0].isBoundary;});
        if(cw_boundary == fan.end())
        {
            PrintoutFan();
            throw std::runtime_error("cw boundary not found");
        }
        std::rotate(fan.begin(), cw_boundary, fan.end());
    }

    // assert that the nodes of the fan connect
    for(std::size_t i = 0;i<fan.size()-1;i++)
    {
        if(fan[i].nd[1] != fan[i+1].nd[0])
        {
            PrintoutFan();
            throw std::runtime_error("fan nodes are not contiguous");
        }
        if(fan[i].e[1].nds[0] != fan[i+1].e[0].nds[0] || fan[i].e[1].nds[1] != fan[i+1].e[0].nds[1])
            throw std::runtime_error("edges not shared");
    }
}

void icy::Node::PrintoutFan()
{
    std::cout << "Printing fan for node " << locId << std::endl;
    std::cout << "fan size " << fan.size() << "; isBoundary " << isBoundary << std::endl;
    std::cout << "adj elems size " << adjacent_elems.size() << std::endl;

    for(Sector &s : fan)
    {
        std::cout << s.nd[0]->locId << "-" << s.nd[1]->locId;
        std::cout << " ; eCW " << s.e[0].nds[0]->locId << "-" << s.e[0].nds[1]->locId << (s.e[0].isBoundary ? " b; " : " nb;");
        std::cout << " ; eCCW " << s.e[1].nds[0]->locId << "-" << s.e[1].nds[1]->locId << (s.e[1].isBoundary ? " b; " : " nb;");
        std::cout << std::endl;
    }
    std::cout << "--------------------------------\n";
    std::cout << std::endl;

}
