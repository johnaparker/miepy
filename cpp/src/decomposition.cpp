#include "decomposition.hpp"

using std::complex;


complex<double> trapz(const Ref<const Array>& x, const Ref<const ComplexArray>& y) {
    int n = y.size();
    complex<double> sum = 0;

    for (int i=1; i < n; i++) {
        sum += (y(i) + y(i-1))*(x(i) - x(i-1))/2.0;
    }

    return sum;
}

complex<double> trapz_2d(const Ref<const Array>& x, const Ref<const Array>& y, const Ref<const ComplexMatrix>& f) {
    int n = y.size();
    ComplexArray z(n);

    for (int i=0; i < n; i++) {
        z(i) = trapz(x, f.col(i));
    }

    return trapz(y, z);
}
