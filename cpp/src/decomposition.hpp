#include "vec.hpp"

std::complex<double> trapz(const Ref<const Array>& x, const Ref<const ComplexArray>& y);
std::complex<double> trapz_2d(const Ref<const Array>& x, const Ref<const Array>& y, const Ref<const ComplexMatrix>& f);
