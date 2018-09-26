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

    ComplexTensor3 tmatrix(2, 16, 16);
    tmatrix.setZero();

    position_t positions = position_t::Zero(2,3);
    positions(0,0) = -1;
    positions(1,0) = 1;

    double k = 1;
    
    particle_aggregate_tmatrix(positions, tmatrix, k);
}
