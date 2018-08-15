#include "vsh_translation.hpp"
#include <iostream>

using namespace std;

int main() {
    int m =1, n=2, u=2, v=2;
    double rad = 0.3, theta = 0.3, phi = 0.3;
    double k = 1;
    auto mode = vsh_mode::outgoing;
    cout << m << endl;

    // auto A = vsh_translation(m, n, u, v, rad, theta, phi, k, mode);
    auto A = vsh_translation(0, 2, 1, 1, 2*M_PI*4/6, M_PI/4, M_PI/4, 1.0, mode);
    
}
