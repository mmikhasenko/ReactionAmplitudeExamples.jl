

// #include <functional>
#include <iostream>
// #include <complex>
#include <vector>
#include <random>
// #include <map>
// #include <iomanip>

#define _USE_MATH_DEFINES
#include <cmath>

#include "Generator.cc"

#define E_BEAM 190.0
#define PROT_MASS 0.938
#define PI_MASS 0.13957
#define PI0_MASS 0.13497
#define RHO_MASS 0.7755

#define PI_MASS_SQ (PI_MASS*PI_MASS)
#define PI0_MASS_SQ (PI0_MASS*PI0_MASS)
#define RHO_MASS_SQ (RHO_MASS*RHO_MASS)
#define PROT_MASS_SQ (PROT_MASS*PROT_MASS)

int main() {

        std::random_device rd;
        std::mt19937 en(rd());
        std::uniform_real_distribution<> dist_cosT(-1.0, 1.0);
        std::uniform_real_distribution<> dist_phi(-M_PI, M_PI);

        double s0 = PI_MASS_SQ+PROT_MASS_SQ + 2*E_BEAM*PROT_MASS_SQ;
        double t = -0.1;  // GeV^2
        double cosTh_from_t = (2*s0*(t-PI_MASS_SQ - RHO_MASS_SQ) +
                               (s0 + PI_MASS_SQ - PROT_MASS_SQ)*
                               (s0+ RHO_MASS_SQ - PROT_MASS_SQ) ) /
                              sqrt( LAMBDA(s0,  PI_MASS_SQ, PROT_MASS_SQ)*
                                    LAMBDA(s0, RHO_MASS_SQ, PROT_MASS_SQ));

        const int Nevents = 100;
        for (int i = 0; i < Nevents; i++) {
                std::vector<double> piProt {0.0, 0.0, sqrt(E_BEAM*E_BEAM-PI_MASS_SQ), E_BEAM+PROT_MASS};
                auto rhoprot = Generator::decay_p(piProt,RHO_MASS_SQ, PROT_MASS_SQ, cosTh_from_t, dist_phi(en));
                auto rho = rhoprot[0];
                auto pions = Generator::decay_p(rho,PI_MASS_SQ, PI0_MASS_SQ, dist_cosT(en), dist_phi(en));
                auto pi0 = pions[1];
                auto gammas = Generator::decay_p(pi0,0.0,0.0, dist_cosT(en), dist_phi(en));
                // print photons
                std::cout << "photon energies: " << gammas[0][3] <<", " << gammas[1][3] << "\n";

                std::cout << "pi0 mass = " << sqrt(Generator::inv_masssq_of_sum(gammas[0],gammas[1])) << "\n";
        }
        return 0.0;
}
