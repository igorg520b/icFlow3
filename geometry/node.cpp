#include <cmath>
#include <cfloat>
#include "node.h"
#include "model.h"
#include "parameters_sim.h"
#include "element.h"

icy::Node::Node()
{
    Reset();
    adjacent_nodes.reserve(10);
    fan.reserve(10);
}




void icy::Node::ComputeElasticForce(SimParams &prms, double timeStep, double totalTime, double added_mass_coeff)
{
    if(prescribed || lsId < 0) return;
    double beta = prms.NewmarkBeta;
    double alpha = prms.HHTalpha;
    double mass = area*prms.Thickness*prms.IceDensity;
    if(mass <= 0) throw std::runtime_error("zero nodal mass");

    F = Eigen::Matrix<double,DOFS,1>::Zero();
    F = at*mass*added_mass_coeff;
//    F(2) -= gravity*mass;
    F(3) = 0;
    F(4) = 0;

    dF = Eigen::Matrix<double,DOFS,DOFS>::Identity()*(mass*added_mass_coeff/(beta*timeStep*timeStep));
    dF(3,3) = 0;
    dF(4,4) = 0;


    vertical_force = 0;

    if(prms.loadType == 1)
    {
        // river rapids
        double spring = area*prms.WaterDensity*std::abs(prms.gravity);
        double riverRapids = 0.5*RiverRapids(x_initial.x(), totalTime*0.3-4);
        double a = prms.Thickness*prms.IceDensity/prms.WaterDensity;
        double c = prms.Thickness - a;

        double disp_t = xt(2)-riverRapids;
        if(disp_t > -c && disp_t < a) dF(2,2) += spring*(1-alpha);
        else if(disp_t < -c) disp_t = -c;
        else if(disp_t > a) disp_t = a;

        double disp_n = xn(2)-riverRapids;
        if(disp_n < -c) disp_n = -c;
        else if(disp_n > a) disp_n = a;

        F(2) += disp_t*spring*(1-alpha);
        F(2) += disp_n*spring*alpha;
        vertical_force = disp_t*spring*(1-alpha) + disp_n*spring*alpha;
    }
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
    else if(prms.loadType == 3)
    {
        // standing wave(s)
        double spring = area*prms.WaterDensity*std::abs(prms.gravity);
//        double rest_position = OceanWave(x_initial.x(), x_initial.y(), totalTime);
        double rest_position = sin(x_initial.x()*1.5);
        rest_position *= 0.07*sin(totalTime * 2 * M_PI / 10);
        if(totalTime < 10) rest_position *= totalTime/10;
        double a = prms.Thickness*prms.IceDensity/prms.WaterDensity;
        double c = prms.Thickness - a;

        double disp_t = xt(2)-rest_position;
        if(disp_t > -c && disp_t < a) dF(2,2) += spring*(1-alpha);
        else if(disp_t < -c) disp_t = -c;
        else if(disp_t > a) disp_t = a;

        double disp_n = xn(2)-rest_position;
        if(disp_n < -c) disp_n = -c;
        else if(disp_n > a) disp_n = a;

        F(2) += disp_t*spring*(1-alpha);
        F(2) += disp_n*spring*alpha;
        vertical_force = disp_t*spring*(1-alpha) + disp_n*spring*alpha;
    }
    else if(prms.loadType == 4)
    {
        // standing wave(s)
        double spring = area*prms.WaterDensity*std::abs(prms.gravity);
//        double rest_position = OceanWave(x_initial.x(), x_initial.y(), totalTime);
        double rest_position = sin(x_initial.x()*1)+sin(x_initial.y()*1+1);
        rest_position *= 0.05*sin(totalTime * 2 * M_PI / 10);
        if(totalTime < 10) rest_position *= totalTime/10;
        double a = prms.Thickness*prms.IceDensity/prms.WaterDensity;
        double c = prms.Thickness - a;

        double disp_t = xt(2)-rest_position;
        if(disp_t > -c && disp_t < a) dF(2,2) += spring*(1-alpha);
        else if(disp_t < -c) disp_t = -c;
        else if(disp_t > a) disp_t = a;

        double disp_n = xn(2)-rest_position;
        if(disp_n < -c) disp_n = -c;
        else if(disp_n > a) disp_n = a;

        F(2) += disp_t*spring*(1-alpha);
        F(2) += disp_n*spring*alpha;
        vertical_force = disp_t*spring*(1-alpha) + disp_n*spring*alpha;
    }
    else if(prms.loadType == 5)
    {
        // river rapids
        double spring = area*prms.WaterDensity*std::abs(prms.gravity);
        double riverRapids = 0.5*RiverRapids(x_initial.x(), totalTime*0.3-25.3);
        double a = prms.Thickness*prms.IceDensity/prms.WaterDensity;
        double c = prms.Thickness - a;

        double disp_t = xt(2)-riverRapids;
        if(disp_t > -c && disp_t < a) dF(2,2) += spring*(1-alpha);
        else if(disp_t < -c) disp_t = -c;
        else if(disp_t > a) disp_t = a;

        double disp_n = xn(2)-riverRapids;
        if(disp_n < -c) disp_n = -c;
        else if(disp_n > a) disp_n = a;

        F(2) += disp_t*spring*(1-alpha);
        F(2) += disp_n*spring*alpha;
        vertical_force = disp_t*spring*(1-alpha) + disp_n*spring*alpha;
    }
    else if(prms.loadType == 6)
    {
        // standing wave(s)
        double spring = area*prms.WaterDensity*std::abs(prms.gravity)*1500;
        double rest_position = prms.wave_height*cos(x_initial.x()*2*M_PI/3.99)*sin(totalTime * 2 * M_PI / 1.6);
        if(totalTime < 2) rest_position *= totalTime/2;
        double a = prms.Thickness*prms.IceDensity/prms.WaterDensity;
        double c = prms.Thickness - a;

        double disp_t = xt(2)-rest_position;
//        if(disp_t > -c && disp_t < a) dF(2,2) += spring*(1-alpha);
//        else if(disp_t < -c) disp_t = -c;
//        else if(disp_t > a) disp_t = a;
        dF(2,2) += spring*(1-alpha);

        double disp_n = xn(2)-rest_position;
//        if(disp_n < -c) disp_n = -c;
//        else if(disp_n > a) disp_n = a;

        F(2) += disp_t*spring*(1-alpha);
        F(2) += disp_n*spring*alpha;
        vertical_force = disp_t*spring*(1-alpha) + disp_n*spring*alpha;
    }
    else if(prms.loadType == 7)
    {
        // standing wave(s)
        double spring = area*prms.WaterDensity*std::abs(prms.gravity)*15;
        double rest_position = prms.wave_height*cos(x_initial.x()*2*M_PI/3.99)*sin(totalTime * 2 * M_PI / 1.6);
        rest_position += prms.wave_height*0.1*cos(x_initial.y()*2*M_PI/5)*sin(2+totalTime * 2 * M_PI / 1.2)*(totalTime < 6 ? totalTime/6 : 1);
        rest_position += prms.wave_height*0.1*cos(x_initial.y()*2*M_PI/10)*sin(1+totalTime * 2 * M_PI / 2)*(totalTime < 4 ? totalTime/4 : 1);
//        rest_position += 0.1*RiverRapids(x_initial.x(), totalTime*0.3-36.5);

        if(totalTime < 2) rest_position *= totalTime/2;
        double a = prms.Thickness*prms.IceDensity/prms.WaterDensity;
        double c = prms.Thickness - a;

        double disp_t = xt(2)-rest_position;
//        if(disp_t > -c && disp_t < a) dF(2,2) += spring*(1-alpha);
//        else if(disp_t < -c) disp_t = -c;
//        else if(disp_t > a) disp_t = a;
        dF(2,2) += spring*(1-alpha);

        double disp_n = xn(2)-rest_position;
//        if(disp_n < -c) disp_n = -c;
//        else if(disp_n > a) disp_n = a;

        F(2) += disp_t*spring*(1-alpha);
        F(2) += disp_n*spring*alpha;
        vertical_force = disp_t*spring*(1-alpha) + disp_n*spring*alpha;
    }


    // damping force
    F(2) += prms.Damping*mass*vt.z()/timeStep;
    dF(2,2) += prms.Damping*mass*prms.NewmarkGamma/(prms.NewmarkBeta*timeStep*timeStep);


    // assert
//    for(int i=0;i<DOFS;i++) if(std::isnan(dF(i))) throw std::runtime_error("node; dF contains NaN");

}

void icy::Node::Reset()
{
    x_initial=ut=xt=vt=at=un=xn=vn=an=Eigen::Matrix<double,DOFS,1>::Zero();
    weakening_direction = normal_n = Eigen::Vector3d::Zero();
    prescribed = false;
    area = 0;
    adjacent_nodes.clear();
    adjacent_edges_map.clear();
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
    if(!prescribed) {
        un = ut;
        xn = xt;
        vn = vt;
        an = at;
    }
}


void icy::Node::PrepareFan()
{
    std::sort(fan.begin(), fan.end(),
              [](icy::Node::FanPrecomp &f0, icy::Node::FanPrecomp &f1)
    {return f0.centerAngle < f1.centerAngle; });

    for(FanPrecomp &f : fan)
    {
        f.e[0] = adjacent_edges_map.at(f.nd[0]->locId);
        f.e[1] = adjacent_edges_map.at(f.nd[1]->locId);
    }

    if(isBoundary)
    {
        // find the fan element with the border on the CW direction
        auto cw_boundary = std::find_if(fan.begin(), fan.end(),
                                        [](FanPrecomp &f){return f.e[0]->isBoundary;});
        if(cw_boundary == fan.end()) throw std::runtime_error("cw boundary not found");
        std::rotate(fan.begin(), cw_boundary, fan.end());
    }

    // assert that nodes of the fan connect
    for(std::size_t i = 0;i<fan.size()-1;i++)
    {
        if(fan[i].nd[1] != fan[i+1].nd[0]) {
            qDebug() << "nd: " << locId << "; fan nodes not contiguous";
            qDebug() << "fan size " << fan.size() << "; isBoundary " << isBoundary;
            qDebug() << "between idx " << i << " and " << i+1;
            qDebug() << "nodes ids: " << fan[i].nd[1]->locId << " and " << fan[i+1].nd[0]->locId;
            qDebug() << "----";
            for(std::size_t j = 0;j<fan.size();j++)
                qDebug() << "j="<<j<<"; nd_ids: " << fan[j].nd[0]->locId << " -- " << fan[j].nd[1]->locId;
            throw std::runtime_error("fan nodes are not contiguous");
        }
        if(fan[i].e[1] != fan[i+1].e[0]) throw std::runtime_error("edges not shared");
    }
}

void icy::Node::ComputeFanVariables(SimParams &prms)
{

    idxSepStressResult = -1;
    max_normal_traction = -DBL_MAX;
    if(fan.size()==0) return;
    auto get_angle = [](Eigen::Vector3d u, Eigen::Vector3d v)
    {
        double dot = u.dot(v)/(u.norm()*v.norm());
        if(dot > 1) dot = 1.0;
        else if(dot < -1.0) dot = -1.0;
        return acos(dot);
    };

    // precompute fan
    double end_angle = 0;
    // str_bottom_avg = Eigen::Matrix3d::Zero();
    normal_n = Eigen::Vector3d::Zero();

    for(FanPrecomp &f : fan) normal_n += f.face->normal_n;
    normal_n.normalize();

    for(FanPrecomp &f : fan)
    {
        Eigen::Vector3d u = f.e[0]->getVec(this);
        Eigen::Vector3d v = f.e[1]->getVec(this);
        f.u_normalized = u.normalized();
        f.v_normalized = v.normalized();

        f.angle0 = end_angle;
        f.angle_span = get_angle(u,v);
        end_angle += f.angle_span;
        f.angle1 = end_angle;

        f.u_p = normal_n.cross(u).normalized();
        f.v_p = normal_n.cross(v).normalized();

        f.t0_top << f.face->str_top * f.u_p;
        f.t1_top << f.face->str_top * f.v_p;
        f.t0_bottom << f.face->str_bottom * f.u_p;
        f.t1_bottom << f.face->str_bottom * f.v_p;

        if(!f.face) throw std::runtime_error("fan elem without face");
//        str_bottom_avg += f.face->str_bottom*f.face->area_initial;
    }
//    str_bottom_avg /= area;

    std::size_t nFans = fan.size();

    dir = Eigen::Vector3d::Zero();

    // discretize (CCW)
    for (std::size_t i=0; i<num_disc; i++) {
        SepStressResult &ssr = sep_stress_results[i];
        ssr.trac_avg_normal = ssr.trac_avg_tangential = 0;
        ssr.tn = Eigen::Vector3d::Zero();
        ssr.traction_top[0] = ssr.traction_top[1] = Eigen::Vector3d::Zero();
        ssr.traction_bottom[0] = ssr.traction_bottom[1] = Eigen::Vector3d::Zero();

        ssr.faces[0] = ssr.faces[1] = nullptr;
        ssr.sep_stress_top = ssr.sep_stress_bottom = 0;

        double angle_fwd = (double)i*end_angle/num_disc;
        ssr.angle_fwd = angle_fwd;

        double half_turn = end_angle/2;
        double angle_bwd = angle_fwd+half_turn;
        if (angle_bwd >= end_angle) angle_bwd -= end_angle;
        ssr.angle_bwd = angle_bwd;

        // integrate traction
        int sector = (isBoundary || angle_fwd < angle_bwd) ? 0 : 1;

        for (std::size_t f=0; f < nFans; f++) {
            FanPrecomp &fp = fan[f];

            if (angle_fwd >= fp.angle0 && angle_fwd < fp.angle1)
            {
                ssr.faces[0] = fp.face;
                ssr.e[0] = fp.e[0];
                ssr.e[1] = fp.e[1];

                double phi = ssr.phi[0] = angle_fwd - fp.angle0;
                ssr.theta[0] = fp.angle1 - angle_fwd;

                double ratio = phi/fp.angle_span;
                ssr.tn = fp.u_normalized*(1-ratio) + fp.v_normalized*ratio;
                ssr.tn.normalize();
                Eigen::Vector3d tn_perp = normal_n.cross(ssr.tn).normalized();
                ssr.tn_perp = tn_perp;
                Eigen::Vector3d tmult_top = fp.face->str_top * tn_perp;
                Eigen::Vector3d tmult_bottom = fp.face->str_bottom * tn_perp;

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

                double phi = ssr.phi[1] = angle_bwd - fp.angle0;
                ssr.theta[1] = fp.angle1 - angle_bwd;

                double ratio = phi/fp.angle_span;
                Eigen::Vector3d tn_bwd = fp.u_normalized*(1-ratio) + fp.v_normalized*ratio;
                Eigen::Vector3d tn_perp = normal_n.cross(tn_bwd).normalized();
                Eigen::Vector3d tmult_top = fp.face->str_top * tn_perp;
                Eigen::Vector3d tmult_bottom = fp.face->str_bottom * tn_perp;

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
/*
        double sx = str_b_top.coeff(0);
        double sy = str_b_top.coeff(1);
        double txy = str_b_top.coeff(2);
        Eigen::Matrix3d str_top;
        str_top << sx, txy, 0, txy, sy, 0, 0, 0, 0;

        sx = str_b_bottom.coeff(0);
        sy = str_b_bottom.coeff(1);
        txy = str_b_bottom.coeff(2);
        Eigen::Matrix3d str_bottom;
        str_bottom << sx, txy, 0, txy, sy, 0, 0, 0, 0;

        Eigen::Vector3d traction_avg_bottom = str_bottom*ssr.tn_perp;
        Eigen::Vector3d traction_avg_top = str_top*ssr.tn_perp;

        double avg_trac_normal_bottom = traction_avg_bottom.dot(ssr.tn);
        double avg_trac_tangential_bottom = traction_avg_bottom.dot(ssr.tn_perp);
        double avg_trac_normal_top = traction_avg_top.dot(ssr.tn);
        double avg_trac_tangential_top = traction_avg_top.dot(ssr.tn_perp);

        ssr.trac_avg_tangential = avg_trac_tangential_bottom;
        ssr.trac_avg_normal = avg_trac_normal_bottom;
*/
        ssr.t0_tangential_top = ssr.traction_top[0].dot(ssr.tn);
        ssr.t1_tangential_top = ssr.traction_top[1].dot(ssr.tn);
        ssr.t0_normal_top = ssr.tn_perp.dot(ssr.traction_top[0]);
        ssr.t1_normal_top = -ssr.tn_perp.dot(ssr.traction_top[1]);
        ssr.trac_normal_top = ssr.t0_normal_top + ssr.t1_normal_top;
        ssr.trac_tangential_top = ssr.t0_tangential_top - ssr.t1_tangential_top;

        ssr.t0_tangential_bottom = ssr.traction_bottom[0].dot(ssr.tn);
        ssr.t1_tangential_bottom = ssr.traction_bottom[1].dot(ssr.tn);
        ssr.t0_normal_bottom = ssr.tn_perp.dot(ssr.traction_bottom[0]);
        ssr.t1_normal_bottom = -ssr.tn_perp.dot(ssr.traction_bottom[1]);
        ssr.trac_normal_bottom = ssr.t0_normal_bottom + ssr.t1_normal_bottom;
        ssr.trac_tangential_bottom = ssr.t0_tangential_bottom - ssr.t1_tangential_bottom;

        if(!isBoundary)
        {
            ssr.trac_normal_bottom /= 2;
            ssr.trac_tangential_bottom /= 2;
            ssr.trac_normal_top /= 2;
            ssr.trac_tangential_top /= 2;
        }

        if(crack_tip)
        {
            double coeff = ((1-prms.weakening_coeff)+(prms.weakening_coeff)*pow((weakening_direction.dot(ssr.tn)+1)/2, 5));
            ssr.trac_normal_bottom*=coeff;
            ssr.trac_normal_top*=coeff;
        }

        if(max_normal_traction < ssr.trac_normal_top) {
            max_normal_traction = ssr.trac_normal_top;
            idxSepStressResult = i;
            dir = ssr.tn;
        }

        if(max_normal_traction < ssr.trac_normal_bottom) {
            max_normal_traction = ssr.trac_normal_bottom;
            idxSepStressResult = i;
            dir = ssr.tn;
        }

    } // num_disc

    // exclude small-angle cuts at the boundary
    if(isBoundary && (idxSepStressResult < num_disc/10 || idxSepStressResult > num_disc*9/10))
        max_normal_traction = 0;
}

icy::Node::FanPrecomp::FanPrecomp(icy::Element *elem, icy::Node *ndd) : face(elem)
{
    // (the list of edges has not been built yet)

    Eigen::Vector3d nd_vec = ndd->x_initial.block(0,0,3,1);
    Eigen::Vector3d tcv = face->getCenter() - nd_vec;
    centerAngle = atan2(tcv.y(), tcv.x());

    nd[0] = face->getCWNode(ndd);
    nd[1] = face->getCCWNode(ndd);
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
    if(x>1.0) x=1.0;
    else if(x<0) x=0;
    return x * x * (3 - 2 * x);
}

double icy::Node::RiverRapids(double x, double t)
{
    return (-Smoothstep(-t, -t+1, x));
}

