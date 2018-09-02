#include "vsh_translation.hpp"
#include "special.hpp"
#include "interactions.hpp"
#include "forces.hpp"
#include <iostream>
#include <complex>

using namespace std;

int main() {
    int lmax = 2;
    int rmax = lmax*(lmax+2);

    auto val = vsh_translation(2, 4, -4, 4, 500, M_PI/2, M_PI/2, 2*M_PI/600, vsh_mode::outgoing);
    int n = 4, m=-2, v = 4, u=-4;
    double phi = M_PI/2;
    complex<double> factor = 0.5*pow(-1, m)*sqrt((2*v+1)*(2*n+1)*factorial(v-u)*factorial(n-m)
            /double(v*(v+1)*n*(n+1)*factorial(v+u)*factorial(n+m)))*exp(1i*double(u+m)*phi);
    uint f = (2*v+1)*(2*n+1)*factorial(6)*factorial(8);
    factor = (2*v+1)*(2*n+1)*factorial(v-u)*factorial(n-m)
            /double(v*(v+1)*n*(n+1)*factorial(v+u)*factorial(n+m));
}
