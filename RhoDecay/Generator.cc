// Copyright [06/2018] Misha Mikhasenko, mikhail.mikhasenko@gmail.com

#include <iostream>
#include <vector>
#include <functional>
#include <cmath>

//------------------------------------------------------------------------------//
//------------------------------------------------------------------------------//
//---------------------GENERATOR OF ARTIGICIAL EVENTS---------------------------//
//------------------------------------------------------------------------------//
//------------------------------------------------------------------------------//

#ifndef LAMBDA
        #define LAMBDA(x,y,z) ((x)*(x)+(y)*(y)+(z)*(z)-2*(x)*(y)-2*(y)*(z)-2*(z)*(x))
#endif

namespace Generator {

double fourv_prod(const std::vector<double> &p1, const std::vector<double> &p2) {
        if (p1.size() != 4 || p2.size() != 4) std::cerr << "Error: p.size() != 4, something is wrong!\n";
        return p1[3]*p2[3]-p1[0]*p2[0]-p1[1]*p2[1]-p1[2]*p2[2];
}
double inv_masssq(const std::vector<double> &p) {
        return fourv_prod(p,p);
}
double inv_masssq_of_sum(const std::vector<double> &p1, const std::vector<double> &p2) {
        return fourv_prod(p1,p1) + fourv_prod(p2,p2) + 2*fourv_prod(p1,p2);
}

std::vector<std::vector<double> > decay_p(const std::vector<double> &p,
                                          double mu1sq, double mu2sq,
                                          double cosTh, double phi = 0.0) {
        double msq = inv_masssq(p);
        double pstar = sqrt(LAMBDA(msq, mu1sq, mu2sq)/(4*msq));
        double sinTh = sqrt(1.0-cosTh*cosTh);
        std::vector<double> rest_frame_q1 = {pstar*sinTh*cos(phi),
                                             pstar*sinTh*sin(phi),
                                             pstar*cosTh,
                                             sqrt(pstar*pstar+mu1sq)};
        std::vector<double> rest_frame_q2 = {-rest_frame_q1[0],
                                             -rest_frame_q1[1],
                                             -rest_frame_q1[2],
                                             sqrt(pstar*pstar+mu2sq)};
        // boost
        auto boost = [](std::vector<double> &p, double gamma) -> void {
                             double gb = sqrt(gamma*gamma-1.0);
                             double pp = p[2]*gamma+p[3]*gb;
                             double E  = p[3]*gamma+p[2]*gb;
                             p[2] = pp; p[3] = E;
                     };
        boost(rest_frame_q1, p[3]/sqrt(msq));
        boost(rest_frame_q2, p[3]/sqrt(msq));

        // rotate
        auto rot = [](std::vector<double> &p, double theta, double phi) -> void {
                           //theta
                           double px =  p[0]*cos(theta)+p[2]*sin(theta);
                           double pz = -p[0]*sin(theta)+p[2]*cos(theta);
                           p[0] = px; p[2] = pz;
                           // phi
                           px = p[0]*cos(phi)-p[1]*sin(phi);
                           double py = p[0]*sin(phi)+p[1]*cos(phi);
                           p[0] = px; p[1] = py;
                   };
        double vec_psq = p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
        double theta_p1 = vec_psq > 0.0 ? acos(p[2]/sqrt(vec_psq)) : 0.0;
        double phi_p1 = atan2(p[1], p[0]);
        rot(rest_frame_q1, theta_p1, phi_p1);
        rot(rest_frame_q2, theta_p1, phi_p1);

        return {rest_frame_q1, rest_frame_q2};
}

std::vector<std::vector<double> > make_up_some_vectors(double s, double cosTh, double phi,
                                                       double m1sq, double m2sq, double m3sq, double m4sq) {
        std::vector<double> p1 = {0.0, 0.0, -sqrt(LAMBDA(s,m1sq,m2sq))/(2*sqrt(s)), (m2sq-s-m1sq)/(2*sqrt(s))};
        std::vector<double> p2 = {0.0, 0.0, -sqrt(LAMBDA(s,m1sq,m2sq))/(2*sqrt(s)), (s-m1sq+m2sq)/(2*sqrt(s))};
        double pstar = sqrt(LAMBDA(s,m3sq,m4sq))/(2*sqrt(s));
        std::vector<double> tot = {0.0, 0.0, 0.0, sqrt(s)};
        std::vector<std::vector<double> > p3p4 = decay_p(tot, m3sq, m4sq, cosTh, phi);
        return {p1,p2,p3p4[0],p3p4[1]};
}

}
