#ifndef GUARD_flux_h
#define GUARD_flux_h

#include "vec.hpp"
#include <complex>
#include <array>

std::array<Array,3> particle_cross_sections(const Ref<const ComplexVector>& p_scat, const Ref<const ComplexVector>& p_inc,
        const Ref<const ComplexVector>& p_src, double k);

std::array<Array,3> cluster_cross_sections(const Ref<const ComplexVector>& p_cluster,
        const Ref<const ComplexVector>& p_src, double k);

#endif
