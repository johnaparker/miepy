#include "indices.hpp"

int rmax_to_lmax(int rmax) {
    int lmax = -1 + int(sqrt(1+rmax));
    return lmax;
}

int lmax_to_rmax(int lmax) {
    return lmax*(lmax + 2);
}
