#include "vsh_translation.hpp"
#include "interactions.hpp"
#include <iostream>

using namespace std;

int main() {
    int Nparticles = 2;
    int lmax = 1;
    int rmax = lmax*(lmax+2);
    double k = 1;

    position_t pos = position_t::Zero(Nparticles,3); 
    pos(1,0) = 5;

    ComplexMatrix mie = ComplexMatrix::Ones(Nparticles, 2*lmax);
    ComplexMatrix T = sphere_aggregate_tmatrix(pos, mie, k);
}
