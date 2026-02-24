#ifndef GUARD_forces_h
#define GUARD_forces_h

#include "vec.hpp"
#include <complex>

vec3 force(const Ref<const ComplexVector>& p_scat, const Ref<const ComplexVector>& p_inc,
        double k, double eps_b = 1, double mu_b = 1);

vec3 torque(const Ref<const ComplexVector>& p_scat, const Ref<const ComplexVector>& p_inc,
        double k, double eps_b = 1, double mu_b = 1);

Matrix force_all(const Ref<const ComplexMatrix>& p_scat_all, const Ref<const ComplexMatrix>& p_inc_all,
        double k, double eps_b = 1, double mu_b = 1);

Matrix torque_all(const Ref<const ComplexMatrix>& p_scat_all, const Ref<const ComplexMatrix>& p_inc_all,
        double k, double eps_b = 1, double mu_b = 1);

#endif
