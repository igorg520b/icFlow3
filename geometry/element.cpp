#include "element.h"
#include "edge.h"
#include "model.h"
#include "node.h"

#include <cstdlib>
#include <algorithm>

// barycentric coords of Gauss points
double icy::Element::N[3][3] = {
       {2.0/3.0, 1.0/6.0, 1.0/6.0},
       {1.0/6.0, 2.0/3.0, 1.0/6.0},
       {1.0/6.0, 1.0/6.0, 2.0/3.0}};

void icy::Element::InitializePersistentVariables()
{
    x_initial << nds[0]->x_initial, nds[1]->x_initial, nds[2]->x_initial;
    // translate the element
    Eigen::Vector3d p1, p2;
    p1 = nds[1]->x_initial.block(0,0,3,1) - nds[0]->x_initial.block(0,0,3,1);
    p2 = nds[2]->x_initial.block(0,0,3,1) - nds[0]->x_initial.block(0,0,3,1);

    rotationMatrix(p1, p2, R0, area_initial, normal_initial);
    for(int j=0;j<3;j++) nds[j]->area += area_initial/3; // distribute area to adjacent nodes

    R0t = R0.transpose();
    pr1_initial = p1;
    pr2_initial = p2;

    // (use 1-based indices, yij = yi-yj)
    double x1 = 0;
    double x2 = pr1_initial.x();
    double x3 = pr2_initial.x();
    double y1 = 0;
    double y2 = pr1_initial.y();
    double y3 = pr2_initial.y();

    double A2 = 2*area_initial;
    double y23, y31, y12, x32, x13, x21;
    y23 = y2-y3;
    y31 = y3-y1;
    y12 = y1-y2;
    x32 = x3-x2;
    x13 = x1-x3;
    x21 = x2-x1;

    // derivatives of shape functions
    double dN1dx = y23/A2;
    double dN2dx = y31/A2;
    double dN3dx = y12/A2;
    double dN1dy = x32/A2;
    double dN2dy = x13/A2;
    double dN3dy = x21/A2;

    // bending strain-displacement matrix
    bmat_b <<
         0,0,0,dN1dx,0,       0,0,0,dN2dx,0,        0,0,0,dN3dx,0,
         0,0,0,0,dN1dy,       0,0,0,0,dN2dy,        0,0,0,0,dN3dy,
         0,0,0,dN1dy,dN1dx,  0,0,0,dN2dy,dN2dx,   0,0,0,dN3dy,dN3dx;

    // membrane strain-displacement matrix
    bmat_m <<
         dN1dx,0,0,0,0,       dN2dx,0,0,0,0,      dN3dx,0,0,0,0,
         0,dN1dy,0,0,0,       0,dN2dy,0,0,0,      0,dN3dy,0,0,0,
         dN1dy,dN1dx,0,0,0,   dN2dy,dN2dx,0,0,0,  dN3dy,dN3dx,0,0,0;

    // shear
    for(int i=0;i<3;i++) {
        bmat_s[i] <<
           0,0,dN1dx,-N[i][0],0,   0,0,dN2dx,-N[i][1],0,     0,0,dN3dx,-N[i][2], 0,
           0,0,dN1dy,0,-N[i][0],   0,0,dN2dy,0,-N[i][1],     0,0,dN3dy,0,-N[i][2];
    }
}

void icy::Element::PrecomputeStiffnessMatrix(icy::SimParams &prms,
                                             Eigen::Matrix3d &elasticityMatrix,
                                             Eigen::Matrix2d &D_mats)
{
    double thickness = prms.Thickness;
    // K and M depend on rho, Young's modulus and Poisson's ratio,
    // therefore they are computed after these parameters are set

    Eigen::Matrix<double,DOFS*3,DOFS*3> K_b, K_m, K_s;
    K = Eigen::Matrix<double,DOFS*3,DOFS*3>::Zero();
    K += bmat_m.transpose()*elasticityMatrix*bmat_m*(area_initial*thickness);

    double coeff = area_initial*thickness*thickness*thickness/12.0;
    K += bmat_b.transpose()*elasticityMatrix*bmat_b*coeff;

    for(int i=0;i<3;i++) K += bmat_s[i].transpose()*D_mats*bmat_s[i]*(thickness*area_initial/3.0);
}

void icy::Element::FdF(
        Eigen::Matrix<double,DOFS*3,1> &x,
        Eigen::Matrix<double,DOFS*3,1> &Fo,
        Eigen::Matrix<double,DOFS*3,DOFS*3> *dFo)
{
    Eigen::Matrix<double,DOFS*3,1> u = x-x_initial;
    Fo = K*u;
    if(dFo != nullptr) *dFo = K;

    /*
    Eigen::Vector3d xc0, xc1, xc2, normal;
    xc0 = x.block(0, 0, 3, 0);
    xc1 = x.block(5, 0, 3, 0);
    xc2 = x.block(10, 0, 3, 0);
    Eigen::Vector3d  p1, p2;
    p1 = xc1-xc0;
    p2 = xc2-xc0;
    Eigen::Matrix3d R, Rt;
    double area;

    rotationMatrix(p1, p2, R, area, normal);
    R = R*R0t;
//    R = Eigen::Matrix3d::Identity();
    Rt = R.transpose();

    Eigen::Vector3d pr1 = Rt*p1;
    Eigen::Vector3d pr2 = Rt*p2;

    Eigen::Vector3d u1, u2;
    u1 = pr1-pr1_initial;
    u2 = pr2-pr2_initial;

    Eigen::Matrix<double,DOFS*3,1> u = x;
    u.block(0,0,3,0) = Eigen::Vector3d::Zero();
    u.block(5,0,3,0) = u1;
    u.block(10,0,3,0) = u2;

    Eigen::Matrix<double,DOFS*3,DOFS*3> R_elem = Eigen::Matrix<double,DOFS*3,DOFS*3>::Identity();
    R_elem.block(0,0,3,3) = R;
    R_elem.block(5,5,3,3) = R;
    R_elem.block(10,10,3,3) = R;

    Fo = R_elem*K*u;
//    Fo = K*u;

    if(dFo != nullptr) {
        *dFo = R_elem*K*R_elem.transpose();
//        *dFo = K;
    }
    */
}

void icy::Element::ComputeElasticForce(icy::LinearSystem &ls, icy::SimParams &prms, double)
{
    if(nds[0]->lsId < 0 && nds[1]->lsId < 0 && nds[2]->lsId < 0) return;

    // reserve non-zero entries in the sparse matrix structure
    ls.AddElementToStructure(nds[0]->lsId, nds[1]->lsId);
    ls.AddElementToStructure(nds[0]->lsId, nds[2]->lsId);
    ls.AddElementToStructure(nds[1]->lsId, nds[2]->lsId);

/*    Eigen::Matrix<double,DOFS*3,1> un, x_initial;
    Eigen::Matrix<double,DOFS*3,1> ut;
    un << nds[0]->un, nds[1]->un, nds[2]->un;
    ut << nds[0]->ut, nds[1]->ut, nds[2]->ut;
    x_initial << nds[0]->x_initial,nds[1]->x_initial,nds[2]->x_initial;
*/
    Eigen::Matrix<double,DOFS*3,1> xn;
    Eigen::Matrix<double,DOFS*3,1> xt;
    xn << nds[0]->xn, nds[1]->xn, nds[2]->xn;
    xt << nds[0]->xt, nds[1]->xt, nds[2]->xt;

    // calculate elastic forces and Hessian at step n+1

    // absolute position is not important, the origin is set at node 0
    Eigen::Matrix<double,DOFS*3,1> Fn, Fnp1;     // internal force at steps n and n+1
    Eigen::Matrix<double,DOFS*3,DOFS*3> dFnp1;   // Hessian of the force function

    FdF(xn, Fn, nullptr);
    FdF(xt, Fnp1, &dFnp1);

    // combine the results into a linearized equation of motion with HHT-alpha integration scheme
    double alpha = prms.HHTalpha;
    F = Fn*alpha + Fnp1*(1-alpha);
    dF= dFnp1*(1-alpha);

    //    double dampingMass = prms.DampingMass;
    //    double dampingStiffness = prms.DampingStiffness;
    //    double beta = prms.NewmarkBeta;
    //    Eigen::Matrix<double,9,9> Cnp1 = M*dampingMass+dFnp1*dampingStiffness;
//    F = M*a + Fn*alpha + Fnp1*(1-alpha) + Cn*vn*alpha + Cnp1*vt*(1-alpha);
//    dF= M/(beta*timeStep*timeStep) + dFnp1*(1-alpha) + Cnp1*((1-alpha)*gamma/(timeStep*beta));

    // assert
    for(int i=0;i<DOFS*3;i++)
        for(int j=0;j<DOFS*3;j++)
            if(std::isnan(dF(i,j)))
                throw std::runtime_error("elem.ComputeElasticForce: dF contains NaN");
}

void icy::Element::Assemble(icy::LinearSystem &ls) const
{
    if(nds[0]->lsId < 0 && nds[1]->lsId < 0 && nds[2]->lsId < 0) return;
    for(int i=0;i<3;i++) {
        int row = nds[i]->lsId;
        Eigen::Matrix<double,DOFS,1> locF = F.block(i*DOFS,0,DOFS,1);
        ls.SubtractRHS(row, locF);
        for(int j=0;j<3;j++) {
            int col = nds[j]->lsId;
            Eigen::Matrix<double,DOFS,DOFS> loc_dF = dF.block(i*DOFS,j*DOFS,DOFS,DOFS);
            ls.AddLHS(row, col, loc_dF);
        }
    }
}

void icy::Element::EvaluateStresses(icy::SimParams &prms,
                                    Eigen::Matrix3d &elasticityMatrix,
                                    Eigen::Matrix2d &D_mats)
{
    Eigen::Matrix<double,DOFS*3,1> u;
    u << nds[0]->ut, nds[1]->ut, nds[2]->ut;

    double thickness = prms.Thickness;
    double coeff = thickness*thickness*thickness/12.0;
    str_b = elasticityMatrix*bmat_b*u*coeff;
    str_m = -elasticityMatrix*bmat_m*u*thickness;
    for(int i=0;i<3;i++) str_s[i] = D_mats*bmat_s[i]*u*(thickness/3.0);

    // stress on top/bottom of the plate and corresponding principal stresses
    str_b_top = elasticityMatrix*bmat_b*u*area_initial*(-thickness/2);

    float sx = (float)str_b_top.coeff(0);
    float sy = str_b_top.coeff(1);
    float txy = str_b_top.coeff(2);
    str_top << sx, txy, txy, sy;

//    double coeff1 = sqrt((sx-sy)*(sx-sy)+txy*txy*4.0);
//    double s1 = 0.5*(sx+sy+coeff1);
//    double s2 = 0.5*(sx+sy-coeff1);
//    str_b_top_principal << s1, s2;

    str_b_bottom = elasticityMatrix*bmat_b*u*area_initial*(thickness/2);

    sx = str_b_bottom.coeff(0) + str_m(0);
    sy = str_b_bottom.coeff(1) + str_m(1);
    txy = str_b_bottom.coeff(2) + str_m(2);
    str_bottom << sx, txy, txy, sy;

//    coeff1 = sqrt((sx-sy)*(sx-sy)+txy*txy*4.0);
//    s1 = 0.5*(sx+sy+coeff1);
//    s2 = 0.5*(sx+sy-coeff1);
//    str_b_bottom_principal << s1, s2;

    ComputeNormal();
}

void icy::Element::DistributeStresses()
{
    // distribute the values from elements to nodes
    for(int i=0;i<3;i++) {
        Node *nd = nds[i];
        double coeff1 = area_initial/nd->area;
        double coeff2 = area_initial/(3*nd->area);

        nd->str_b += str_b*coeff1;
        nd->str_m += str_m*coeff1;

        for(int j=0;j<3;j++) nd->str_s += str_s[j]*coeff2*N[i][j];

        nd->str_b_top += str_b_top*coeff1;
        nd->str_b_bottom += str_b_bottom*coeff1;
        nd->str_b_top_principal += str_b_top_principal*coeff1;
        nd->str_b_bottom_principal += str_b_bottom_principal*coeff1;

        nd->normal_n += normal_n;
    }
}

void icy::Element::rotationMatrix(Eigen::Vector3d &p1, Eigen::Vector3d &p2, Eigen::Matrix3d &result,
                                   double &area, Eigen::Vector3d &normal)
{
    Eigen::Vector3d r1 = p1.normalized();
    normal = p1.cross(p2); area=normal.norm()/2; normal.normalize();
    Eigen::Vector3d r2 = normal.cross(r1).normalized();
    result << r1, r2, normal;
}

void icy::Element::rotationMatrix_alt(Eigen::Vector3d &p12, Eigen::Vector3d &p13,
                        Eigen::Matrix3d &result,
                        double &area, Eigen::Vector3d &normal)
{
    normal = p12.cross(p13); area=normal.norm()/2; normal.normalize();
    double nx = normal.coeff(0);
    double ny = normal.coeff(0);
    double nz = normal.coeff(2);

    Eigen::Vector3d r1;
    if(nx == 0 && ny == 0)
    {
        r1 << 1, 0, 0;
    }
    else if(nx != 0 && nz != 0)
    {
        double c1 = nx/nz;
        double c1sq = c1*c1;
        r1 << 1.0/sqrt(1.0+c1sq), 0, -1.0/sqrt(1.0+1.0/c1sq);
        r1.normalize();
    }
    else {
        r1 = p12.normalized();
    }

    Eigen::Vector3d r2 = normal.cross(r1).normalized();
    result << r1, r2, normal;
}


void icy::Element::ComputeNormal()
{
    Eigen::Vector3d u = (nds[1]->xn - nds[0]->xn).block(0,0,3,1);
    Eigen::Vector3d v = (nds[2]->xn - nds[0]->xn).block(0,0,3,1);
    normal_n = u.cross(v);
}

/*
template <typename T> T sqr (const T& x) { return x*x; }

int icy::Element::dsyevc3(const Eigen::Matrix3d &A, Eigen::Vector3d &w) {
  double m, c1, c0;

  // Determine coefficients of characteristic poynomial. We write
  //       | a   d   f  |
  //  A =  | d*  b   e  |
  //       | f*  e*  c  |
  double de = A(0,1) * A(1,2);                                    // d * e
  double dd = sqr(A(0,1));                                         // d^2
  double ee = sqr(A(1,2));                                         // e^2
  double ff = sqr(A(0,2));                                         // f^2
  m  = A(0,0) + A(1,1) + A(2,2);
  c1 = (A(0,0)*A(1,1) + A(0,0)*A(2,2) + A(1,1)*A(2,2))        // a*b + a*c + b*c - d^2 - e^2 - f^2
          - (dd + ee + ff);
  c0 = A(2,2)*dd + A(0,0)*ee + A(1,1)*ff - A(0,0)*A(1,1)*A(2,2)
            - 2.0 * A(0,2)*de;                                     // c*d^2 + a*e^2 + b*f^2 - a*b*c - 2*f*d*e)

  double p, sqrt_p, q, c, s, phi;
  p = sqr(m) - 3.0*c1;
  q = m*(p - (3.0/2.0)*c1) - (27.0/2.0)*c0;
  sqrt_p = sqrt(fabs(p));

  phi = 27.0 * ( 0.25*sqr(c1)*(p - c1) + c0*(q + 27.0/4.0*c0));
  phi = (1.0/3.0) * atan2(sqrt(fabs(phi)), q);

  c = sqrt_p*cos(phi);
  s = (1.0/1.73205080756887729352744634151)*sqrt_p*sin(phi);

  w[1]  = (1.0/3.0)*(m - c);
  w[2]  = w[1] + s;
  w[0]  = w[1] + c;
  w[1] -= s;

  return 0;
}


void icy::Element::swap_columns(Eigen::Matrix3d &A, int col1, int col2)
{
    Eigen::Vector3d cv1 = A.col(col1);
    Eigen::Vector3d cv2 = A.col(col2);
    A.col(col1) = cv2;
    A.col(col2) = cv1;
}

void icy::Element::swap_rows(Eigen::Vector3d &w, int col1, int col2)
{
    double cv1 = w[col1];
    double cv2 = w[col2];
    w[col1] = cv2;
    w[col2] = cv1;
}

// http://www.mpi-hd.mpg.de/personalhomes/globes/3x3
void icy::Element::eigen_decomposition(const Eigen::Matrix3d &B,
                         Eigen::Matrix3d &Q, Eigen::Vector3d &l)
{
  Eigen::Matrix3d A = B;

  double norm;          // Squared norm or inverse norm of current eigenvector
  double n0, n1;        // Norm of first and second columns of A
  double n0tmp, n1tmp;  // "Templates" for the calculation of n0/n1 - saves a few FLOPS
  double thresh;        // Small number used as threshold for floating point comparisons
  double error;         // Estimated maximum roundoff error in some steps
  double wmax;          // The eigenvalue of maximum modulus
  double f, t;          // Intermediate storage
  int i, j;             // Loop counters

  // Calculate eigenvalues
  dsyevc3(A, l);

  wmax = fabs(l[0]);
  if ((t=fabs(l[1])) > wmax)
    wmax = t;
  if ((t=fabs(l[2])) > wmax)
    wmax = t;
  thresh = sqr(8.0 * DBL_EPSILON * wmax);

  // Prepare calculation of eigenvectors
  n0tmp   = sqr(A(0,1)) + sqr(A(0,2));
  n1tmp   = sqr(A(0,1)) + sqr(A(1,2));
  Q(0,1) = A(0,1)*A(1,2) - A(0,2)*A(1,1);
  Q(1,1) = A(0,2)*A(0,1) - A(1,2)*A(0,0);
  Q(2,1) = sqr(A(0,1));

  // Calculate first eigenvector by the formula
  //   v[0] = (A - e.l[0]).e1 x (A - e.l[0]).e2
  A(0,0) -= l[0];
  A(1,1) -= l[0];
  Q(0,0) = Q(0,1) + A(0,2)*l[0];
  Q(1,0) = Q(1,1) + A(1,2)*l[0];
  Q(2,0) = A(0,0)*A(1,1) - Q(2,1);
  norm    = sqr(Q(0,0)) + sqr(Q(1,0)) + sqr(Q(2,0));
  n0      = n0tmp + sqr(A(0,0));
  n1      = n1tmp + sqr(A(1,1));
  error   = n0 * n1;

  if (n0 <= thresh)         // If the first column is zero, then (1,0,0) is an eigenvector
  {
    Q(0,0) = 1.0;
    Q(1,0) = 0.0;
    Q(2,0) = 0.0;
  }
  else if (n1 <= thresh)    // If the second column is zero, then (0,1,0) is an eigenvector
  {
    Q(0,0) = 0.0;
    Q(1,0) = 1.0;
    Q(2,0) = 0.0;
  }
  else if (norm < sqr(64.0 * DBL_EPSILON) * error)
  {                         // If angle between A[0] and A[1] is too small, don't use
    t = sqr(A(0,1));       // cross product, but calculate v ~ (1, -A0/A1, 0)
    f = -A(0,0) / A(0,1);
    if (sqr(A(1,1)) > t)
    {
      t = sqr(A(1,1));
      f = -A(0,1) / A(1,1);
    }
    if (sqr(A(1,2)) > t)
      f = -A(0,2) / A(1,2);
    norm    = 1.0/sqrt(1 + sqr(f));
    Q(0,0) = norm;
    Q(1,0) = f * norm;
    Q(2,0) = 0.0;
  }
  else                      // This is the standard branch
  {
    norm = sqrt(1.0 / norm);
    for (j=0; j < 3; j++)
      Q(j,0) = Q(j,0) * norm;
  }


  // Prepare calculation of second eigenvector
  t = l[0] - l[1];
  if (fabs(t) > 8.0 * DBL_EPSILON * wmax)
  {
    // For non-degenerate eigenvalue, calculate second eigenvector by the formula
    //   v[1] = (A - e.l[1]).e1 x (A - e.l[1]).e2
    A(0,0) += t;
    A(1,1) += t;
    Q(0,1)  = Q(0,1) + A(0,2)*l[1];
    Q(1,1)  = Q(1,1) + A(1,2)*l[1];
    Q(2,1)  = A(0,0)*A(1,1) - Q(2,1);
    norm     = sqr(Q(0,1)) + sqr(Q(1,1)) + sqr(Q(2,1));
    n0       = n0tmp + sqr(A(0,0));
    n1       = n1tmp + sqr(A(1,1));
    error    = n0 * n1;

    if (n0 <= thresh)       // If the first column is zero, then (1,0,0) is an eigenvector
    {
      Q(0,1) = 1.0;
      Q(1,1) = 0.0;
      Q(2,1) = 0.0;
    }
    else if (n1 <= thresh)  // If the second column is zero, then (0,1,0) is an eigenvector
    {
      Q(0,1) = 0.0;
      Q(1,1) = 1.0;
      Q(2,1) = 0.0;
    }
    else if (norm < sqr(64.0 * DBL_EPSILON) * error)
    {                       // If angle between A[0] and A[1] is too small, don't use
      t = sqr(A(0,1));     // cross product, but calculate v ~ (1, -A0/A1, 0)
      f = -A(0,0) / A(0,1);
      if (sqr(A(1,1)) > t)
      {
        t = sqr(A(1,1));
        f = -A(0,1) / A(1,1);
      }
      if (sqr(A(1,2)) > t)
        f = -A(0,2) / A(1,2);
      norm    = 1.0/sqrt(1 + sqr(f));
      Q(0,1) = norm;
      Q(1,1) = f * norm;
      Q(2,1) = 0.0;
    }
    else
    {
      norm = sqrt(1.0 / norm);
      for (j=0; j < 3; j++)
        Q(j,1) = Q(j,1) * norm;
    }
  }
  else
  {
    // For degenerate eigenvalue, calculate second eigenvector according to
    //   v[1] = v[0] x (A - e.l[1]).e[i]
    //
    // This would really get to complicated if we could not assume all of A to
    // contain meaningful values.
    A(1,0)  = A(0,1);
    A(2,0)  = A(0,2);
    A(2,1)  = A(1,2);
    A(0,0) += l[0];
    A(1,1) += l[0];
    for (i=0; i < 3; i++)
    {
      A(i,i) -= l[1];
      n0       = sqr(A(0,i)) + sqr(A(1,i)) + sqr(A(2,i));
      if (n0 > thresh)
      {
        Q(0,1)  = Q(1,0)*A(2,i) - Q(2,0)*A(1,i);
        Q(1,1)  = Q(2,0)*A(0,i) - Q(0,0)*A(2,i);
        Q(2,1)  = Q(0,0)*A(1,i) - Q(1,0)*A(0,i);
        norm     = sqr(Q(0,1)) + sqr(Q(1,1)) + sqr(Q(2,1));
        if (norm > sqr(256.0 * DBL_EPSILON) * n0) // Accept cross product only if the angle between
        {                                         // the two vectors was not too small
          norm = sqrt(1.0 / norm);
          for (j=0; j < 3; j++)
            Q(j,1) = Q(j,1) * norm;
          break;
        }
      }
    }

    if (i == 3)    // This means that any vector orthogonal to v[0] is an EV.
    {
      for (j=0; j < 3; j++)
        if (Q(j,0) != 0.0)                                   // Find nonzero element of v[0] ...
        {                                                     // ... and swap it with the next one
          norm          = 1.0 / sqrt(sqr(Q(j,0)) + sqr(Q((j+1)%3,0)));
          Q(j,1)       = Q((j+1)%3,0) * norm;
          Q((j+1)%3,1) = -Q(j,0) * norm;
          Q((j+2)%3,1) = 0.0;
          break;
        }
    }
  }

  // Calculate third eigenvector according to
  //   v[2] = v[0] x v[1]
  Q(0,2) = Q(1,0)*Q(2,1) - Q(2,0)*Q(1,1);
  Q(1,2) = Q(2,0)*Q(0,1) - Q(0,0)*Q(2,1);
  Q(2,2) = Q(0,0)*Q(1,1) - Q(1,0)*Q(0,1);

  // sort eigenvectors
  if (l[1] > l[0]) { swap_columns(Q, 0, 1); swap_rows(l, 0, 1); }
  if (l[2] > l[0]) { swap_columns(Q, 0, 2); swap_rows(l, 0, 2); }
  if (l[2] > l[1]) { swap_columns(Q, 1, 2); swap_rows(l, 1, 2); }
}

Eigen::Matrix3d icy::Element::make_positive(const Eigen::Matrix3d B)
{
    Eigen::Matrix3d Q;
    Eigen::Vector3d l;

    eigen_decomposition(B, Q, l);

    for(int i=0;i<3;i++) if(l[i] < 0) l[i] = 0;
    Eigen::DiagonalMatrix<double,3> diag (l);
    return Q*diag*Q.transpose();
}
*/

//=============== fracture helpers
icy::Node* icy::Element::getOppositeNode(icy::Edge *edge)
{
    icy::Node *nd0 = edge->nds[0];
    icy::Node *nd1 = edge->nds[1];

    for(int i=0;i<3;i++)
    {
        int idx_next = (i+1)%3;
        if((nds[i] == nd0 && nds[idx_next] == nd1)||
                (nds[i] == nd1 && nds[idx_next] == nd0))
            return nds[(i+2)%3];
    }
    throw std::runtime_error("opposite node not found");
}

std::pair<int,int> icy::Element::getOppositeEdge(icy::Node *nd)
{
    int nd_idx;
    if(nd == nds[0]) nd_idx = 0;
    else if(nd == nds[1]) nd_idx = 1;
    else if(nd == nds[2]) nd_idx = 2;
    else throw std::runtime_error("node not found");

    int edge_idx0 = nds[(nd_idx+1)%3]->locId;
    int edge_idx1 = nds[(nd_idx+2)%3]->locId;
    if(edge_idx0<edge_idx1) return std::make_pair(edge_idx0, edge_idx1);
    else return std::make_pair(edge_idx1, edge_idx0);
}

void icy::Element::ReplaceNode(icy::Node *replaceWhat, icy::Node *replaceWith)
{
    if(nds[0] == replaceWhat) { nds[0] = replaceWith; return; }
    if(nds[1] == replaceWhat) { nds[1] = replaceWith; return; }
    if(nds[2] == replaceWhat) { nds[2] = replaceWith; return; }
    throw std::runtime_error("replaceWhat is not in nds[]");
}

void icy::Element::Initialize(icy::Node* nd0, icy::Node* nd1, icy::Node* nd2, bool orientation)
{
    nds[0] = nd0;
    if(orientation) {
        nds[1] = nd1;
        nds[2] = nd2;
    } else {
        nds[2] = nd1;
        nds[1] = nd2;
    }
    InitializePersistentVariables();
}

bool icy::Element::Orientation(icy::Node* nd0, icy::Node* nd1)
{
    if(nds[0] == nd0 && nds[1] == nd1) return true;
    if(nds[1] == nd0 && nds[0] == nd1) return false;

    if(nds[1] == nd0 && nds[2] == nd1) return true;
    if(nds[2] == nd0 && nds[1] == nd1) return false;

    if(nds[2] == nd0 && nds[0] == nd1) return true;
    if(nds[0] == nd0 && nds[2] == nd1) return false;

    throw std::runtime_error("nodes do not belong to this element");
}

Eigen::Vector3d icy::Element::getCenter()
{
    return (nds[0]->x_initial + nds[1]->x_initial + nds[2]->x_initial).block(0,0,3,1)/3.0;
}

icy::Node* icy::Element::getCCWNode(icy::Node* nd)
{
    if(normal_initial.z()<0) {
        if(nd==nds[0]) return nds[1];
        if(nd==nds[1]) return nds[2];
        if(nd==nds[2]) return nds[0];
    } else
    {
        if(nd==nds[0]) return nds[2];
        if(nd==nds[1]) return nds[0];
        if(nd==nds[2]) return nds[1];
    }
    throw std::runtime_error("nd not found");
}

icy::Node* icy::Element::getCWNode(icy::Node* nd)
{
    if(normal_initial.z()<0) {

        if(nd==nds[0]) return nds[2];
        if(nd==nds[1]) return nds[0];
        if(nd==nds[2]) return nds[1];
    } else {
        if(nd==nds[0]) return nds[1];
        if(nd==nds[1]) return nds[2];
        if(nd==nds[2]) return nds[0];
    }
    throw std::runtime_error("nd not found");
}
