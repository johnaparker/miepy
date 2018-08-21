#include "vsh_translation.hpp"
#include "interactions.hpp"
#include "forces.hpp"
#include <iostream>
#include <complex>

using namespace std;

int main() {
    int lmax = 2;
    int rmax = lmax*(lmax+2);

    ComplexVector p_scat = ComplexVector::Ones(2*rmax); 
    ComplexVector p_inc  = ComplexVector::Ones(2*rmax); 
    for (int i = 0; i < p_scat.size(); i++) {
        p_scat(i) += complex<double>(-.2*i,.4*i);
        p_inc(i) += complex<double>(-.3*i,.3*i);
    }

    vec3 F = force(p_scat, p_inc, 1.0);
}
