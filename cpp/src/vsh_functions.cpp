#include "vsh_functions.hpp"
#include "special.hpp"
#include "indices.hpp"

using std::complex;
using namespace std::complex_literals;


vec3 rad_hat(double theta, double phi) {
    vec3 ret;
    ret << sin(theta)*cos(phi),
           sin(theta)*sin(phi),
           cos(theta);
    return ret;
}

vec3 theta_hat(double theta, double phi) {
    vec3 ret;
    ret << cos(theta)*cos(phi),
           cos(theta)*sin(phi),
           -sin(theta);
    return ret;
}

vec3 phi_hat(double phi) {
    vec3 ret;
    ret << -sin(phi),
           cos(phi),
           0;
    return ret;
}

cvec3 vec_sph_to_cart(const Ref<const cvec3>& v, double theta, double phi) {
    cvec3 ret;
    ret << v.dot(rad_hat(theta, phi)),
           v.dot(theta_hat(theta, phi)),
           v.dot(phi_hat(phi));

    return ret;
}

std::function<complex<double>(int, double, bool)> get_zn(vsh_mode mode) {
    switch(mode) {
        case vsh_mode::outgoing:
            return spherical_hn;
        case vsh_mode::ingoing:
            return spherical_hn_2;
        case vsh_mode::incident:
        case vsh_mode::interior:
            return spherical_jn;
        default:
            return spherical_hn_2;
    }
}

complex<double> Emn(int m, int n) {
    return pow(1i, n)*sqrt((2*n+1)*factorial(n-m)/(n*(n+1)*factorial(n+m)));
}

cvec3 vsh_electric(int n, int m, vsh_mode mode, double rad,
        double theta, double phi, double k) {

    auto zn   = get_zn(mode);
    double kr = k*rad;
    complex<double> exp_phi = exp(1i*double(m)*phi);

    complex<double> H  = zn(n, kr, false);
    complex<double> Hp = zn(n, kr, true);
    double Pnm = associated_legendre(n, m, cos(theta));

    complex<double> factor     = (H + kr*Hp)*exp_phi/kr;
    complex<double> r_comp     = n*(n+1)*Pnm*H/(kr)*exp_phi;
    complex<double> theta_comp = factor*tau_func(n, m, theta);
    complex<double> phi_comp   = 1i*factor*pi_func(n, m, theta);

    cvec3 E;
    E << r_comp, theta_comp, phi_comp;
    return E;
}

cvec3 vsh_magnetic(int n, int m, vsh_mode mode, double rad,
        double theta, double phi, double k) {

    auto zn   = get_zn(mode);
    double kr = k*rad;
    complex<double> exp_phi = exp(1i*double(m)*phi);

    complex<double> H  = zn(n, kr, false);

    complex<double> factor     = H*exp_phi;
    complex<double> r_comp     = 0;
    complex<double> theta_comp = 1i*factor*pi_func(n, m, theta);
    complex<double> phi_comp   = -factor*tau_func(n, m, theta);

    cvec3 E;
    E << r_comp, theta_comp, phi_comp;
    return E;
}

E_type expand_E_cluster(const Ref<const position_t>& pos, const Ref<const ComplexMatrix>& p_expand,
        vsh_mode mode, const Ref<const Array>& x, const Ref<const Array>& y,
        const Ref<const Array>& z, double k) {

    int Npts = x.size();
    E_type E_field = E_type::Zero(Npts, 3);

    // shape(p) = Nparticles x 2 x rmax
    int Nparticles = p_expand.rows();
    int rmax = p_expand.cols()/2;
    int lmax = rmax_to_lmax(rmax);

    complex<double> factor = (mode == vsh_mode::outgoing) ? 1i : -1i;

    for (int j = 0; j < Npts; j++) {
        for (int i = 0; i < Nparticles; i++) {
            double xr = x(j) - pos(i,0);
            double yr = y(j) - pos(i,1);
            double zr = z(j) - pos(i,2);

            double radius = sqrt(pow(xr, 2) + pow(yr, 2) + pow(zr, 2));
            double theta  = acos(zr/radius);
            double phi    = atan2(yr, xr);

            for (int n = 1; n <= lmax; n++) {
                for (int m = -n; m <= n; m++) {
                    complex<double> Emn_val = Emn(m, n);

                    int idx = n*(n+2) - n + m - 1;
                    cvec3 E_sph = vsh_electric(n, m, mode, radius, theta, phi, k);
                    cvec3 E_cart = vec_sph_to_cart(E_sph, theta, phi);
                    E_field.row(j) += -factor*Emn_val*p_expand(i,idx)*E_cart;

                    idx += rmax;
                    E_sph = vsh_magnetic(n, m, mode, radius, theta, phi, k);
                    E_cart = vec_sph_to_cart(E_sph, theta, phi);
                    E_field.row(j) += -factor*Emn_val*p_expand(i,idx)*E_cart;
                }
            }
        }
    }

    return E_field;
}
