#include "decomposition.hpp"

double trapz(const Ref<const Array>& x, const Ref<const Array>& y) {
    int n = y.size();
    double sum = 0;

    for (int i=1; i < n; i++) {
        sum += (y(i) + y(i-1))*(x(i) - x(i-1))/2;
    }

    return sum;
}

double trapz_2d(const Ref<const Array>& x, const Ref<const Array>& y, const Ref<const Matrix>& f) {
    int n = y.size();
    Array z(n);

    for (int i=0; i < n; i++) {
        z(i) = trapz(x, f.row(i));
    }

    return trapz(y, z);
}
