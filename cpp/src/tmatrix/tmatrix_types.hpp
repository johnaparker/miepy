#ifndef GUARD_tmatrix_types_h
#define GUARD_tmatrix_types_h

// MIEPY_HAS_QUAD is set by CMake when libquadmath is found and linked.
// Do not auto-detect here to avoid compiling quadmath calls without linking the library.
#ifndef MIEPY_HAS_QUAD
    #define MIEPY_HAS_QUAD 0
#endif

// __float128 math overloads MUST come before <complex> to resolve ambiguity
// in std::abs(std::complex<__float128>) which internally uses unqualified sqrt()
#if MIEPY_HAS_QUAD
#include <quadmath.h>
#include <cmath>
namespace std {
    inline __float128 sqrt(__float128 __x) { return sqrtq(__x); }
}
#endif

#include <complex>
#include <cmath>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>

#if MIEPY_HAS_QUAD
// Eigen NumTraits specializations for __float128
namespace Eigen {

template<> struct NumTraits<__float128> : GenericNumTraits<__float128> {
    typedef __float128 Real;
    typedef __float128 NonInteger;
    typedef __float128 Nested;

    enum {
        IsComplex = 0,
        IsInteger = 0,
        IsSigned = 1,
        RequireInitialization = 1,
        ReadCost = 2,
        AddCost = 2,
        MulCost = 6
    };

    static inline __float128 epsilon() { return FLT128_EPSILON; }
    static inline __float128 dummy_precision() { return FLT128_EPSILON * 100; }
    static inline __float128 highest() { return FLT128_MAX; }
    static inline __float128 lowest() { return -FLT128_MAX; }
    static inline int digits10() { return FLT128_DIG; }
};

template<> struct NumTraits<std::complex<__float128>> : GenericNumTraits<std::complex<__float128>> {
    typedef __float128 Real;
    typedef std::complex<__float128> NonInteger;
    typedef std::complex<__float128> Nested;

    enum {
        IsComplex = 1,
        IsInteger = 0,
        IsSigned = 1,
        RequireInitialization = 1,
        ReadCost = 4,
        AddCost = 4,
        MulCost = 12
    };

    static inline __float128 epsilon() { return FLT128_EPSILON; }
    static inline __float128 dummy_precision() { return FLT128_EPSILON * 100; }
};

} // namespace Eigen
#endif // MIEPY_HAS_QUAD

namespace tmatrix {

// Pi constant with full precision
template<typename Real>
constexpr Real pi() {
    return Real(3.14159265358979323846264338327950288419716939937510L);
}

// Type trait bundle for a given real precision
template<typename Real>
struct Types {
    using real_t = Real;
    using complex_t = std::complex<Real>;
    using Matrix = Eigen::Matrix<complex_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using Vector = Eigen::Matrix<complex_t, Eigen::Dynamic, 1>;
    using RealVector = Eigen::Matrix<Real, Eigen::Dynamic, 1>;
    using RealMatrix = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
};

// Convenience aliases for double precision
using DoubleTypes = Types<double>;

// Math helpers that work for both double and __float128
template<typename Real>
inline Real tm_sqrt(Real x) { return std::sqrt(x); }

template<typename Real>
inline Real tm_abs(Real x) { return std::abs(x); }

template<typename Real>
inline Real tm_cos(Real x) { return std::cos(x); }

template<typename Real>
inline Real tm_sin(Real x) { return std::sin(x); }

template<typename Real>
inline Real tm_atan(Real x) { return std::atan(x); }

template<typename Real>
inline Real tm_atan2(Real y, Real x) { return std::atan2(y, x); }

template<typename Real>
inline std::complex<Real> tm_exp(std::complex<Real> z) { return std::exp(z); }

template<typename Real>
inline std::complex<Real> tm_sqrt(std::complex<Real> z) { return std::sqrt(z); }

template<typename Real>
inline Real tm_abs(std::complex<Real> z) { return std::abs(z); }

template<typename Real>
inline std::complex<Real> tm_sin(std::complex<Real> z) { return std::sin(z); }

template<typename Real>
inline std::complex<Real> tm_cos(std::complex<Real> z) { return std::cos(z); }

// Machine epsilon and smallest positive normal
template<typename Real>
inline Real machine_eps() { return std::numeric_limits<Real>::epsilon(); }

template<typename Real>
inline Real smallest_pos() { return std::numeric_limits<Real>::min(); }

#if MIEPY_HAS_QUAD

// Specialize math for __float128
template<> inline __float128 tm_sqrt(__float128 x) { return sqrtq(x); }
template<> inline __float128 tm_abs(__float128 x) { return fabsq(x); }
template<> inline __float128 tm_cos(__float128 x) { return cosq(x); }
template<> inline __float128 tm_sin(__float128 x) { return sinq(x); }
template<> inline __float128 tm_atan(__float128 x) { return atanq(x); }
template<> inline __float128 tm_atan2(__float128 y, __float128 x) { return atan2q(y, x); }

// Complex __float128 math
template<>
inline std::complex<__float128> tm_exp(std::complex<__float128> z) {
    __float128 r = expq(z.real());
    return std::complex<__float128>(r * cosq(z.imag()), r * sinq(z.imag()));
}

template<>
inline std::complex<__float128> tm_sqrt(std::complex<__float128> z) {
    __float128 r = sqrtq(z.real() * z.real() + z.imag() * z.imag());
    __float128 mag = sqrtq(r);
    __float128 angle = atan2q(z.imag(), z.real()) / 2;
    return std::complex<__float128>(mag * cosq(angle), mag * sinq(angle));
}

template<>
inline __float128 tm_abs(std::complex<__float128> z) {
    return sqrtq(z.real() * z.real() + z.imag() * z.imag());
}

template<>
inline std::complex<__float128> tm_sin(std::complex<__float128> z) {
    return std::complex<__float128>(
        sinq(z.real()) * coshq(z.imag()),
        cosq(z.real()) * sinhq(z.imag())
    );
}

template<>
inline std::complex<__float128> tm_cos(std::complex<__float128> z) {
    return std::complex<__float128>(
        cosq(z.real()) * coshq(z.imag()),
        -sinq(z.real()) * sinhq(z.imag())
    );
}

// Machine epsilon and smallest positive for __float128
template<> inline __float128 machine_eps<__float128>() { return FLT128_EPSILON; }
template<> inline __float128 smallest_pos<__float128>() { return FLT128_MIN; }

using QuadTypes = Types<__float128>;

#endif // MIEPY_HAS_QUAD

} // namespace tmatrix

#endif // GUARD_tmatrix_types_h
