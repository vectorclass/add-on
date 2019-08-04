/***************************  quaternion.h   *********************************
* Author:        Agner Fog
* Date created:  2012-08-01
* Last modified: 2019-07-13
* Version:       2.00
* Project:       Extension to vector class library
* Description:
* Quaternions are used in theoretical algebra
* Classes for quaternions:
* Quaternion1f:  One quaternion consisting of four single precision floats
* Quaternion1d:  One quaternion consisting of four double precision floats
*
* (c) Copyright 2012-2019 Apache License version 2.0 or later
******************************************************************************/


#ifndef QUATERNION_H
#define QUATERNION_H  200

#include "vectorclass.h"
#include <cmath>

#ifdef VCL_NAMESPACE
namespace VCL_NAMESPACE {
#endif

/*****************************************************************************
*
*                     Class Quaternion1f
*         One quaternion consisting of four single precision floats
*
*****************************************************************************/

class Quaternion1f {
protected:
    __m128 xmm; // vector of 4 single precision floats
public:
    // default constructor
    Quaternion1f() {
    }
    // construct from real, no imaginary part
    Quaternion1f(float re) {
        xmm = _mm_load_ss(&re);
    }
    // construct from real and imaginary parts = re + im0*i + im1*j + im2*k
    Quaternion1f(float re, float im0, float im1, float im2) {
        xmm = Vec4f(re, im0, im1, im2);
    }
    // Constructor to convert from type __m128 used in intrinsics:
    Quaternion1f(__m128 const x) {
        xmm = x;
    }
    // Assignment operator to convert from type __m128 used in intrinsics:
    Quaternion1f & operator = (__m128 const x) {
        xmm = x;
        return *this;
    }
    // Constructor to convert from Vec4f
    Quaternion1f(Vec4f const x) {
        xmm = x;
    }
    // Type cast operator to convert to __m128 used in intrinsics
    operator __m128() const {
        return xmm;
    }
    // Member function to convert to vector
    Vec4f to_vector() const {
        return xmm;
    }
    // Member function to load from array
    Quaternion1f & load(float const * p) {
        xmm = Vec4f().load(p);
        return *this;
    }
    // Member function to store into array
    void store(float * p) const {
        Vec4f(xmm).store(p);
    }
    // Member function to extract real part
    float real() const {
        return _mm_cvtss_f32(xmm);
    }
    // Member function to extract imaginary parts, sets real part to 0
    Quaternion1f imag() const {
        return Quaternion1f(permute4<-1,1,2,3>(Vec4f(xmm)));
    }
#ifdef COMPLEXVEC_H  // relations to complexvec1.h
    // construct from two Complex1f = a0 + a1 * j
    Quaternion1f(Complex1f const a0, Complex1f const a1) {
        xmm = _mm_movelh_ps(a0, a1);
    }
    // Member functions to split into two Complex1f:
    // q = q.get_low() + q.get_high()*j
    Complex1f get_low() const {
        return Complex1f(Vec4f(xmm).cutoff(2));
    }
    Complex1f get_high() const {
        __m128 t = _mm_movehl_ps(_mm_setzero_ps(), xmm);
        return Complex1f(t);
    } 
#endif
#ifdef VECTOR3D_H   // relations to vector3d.h
    // Constructor to convert from Vec3f used in geometrics:
    Quaternion1f(Vec3f const x) {
        xmm = permute4<3,0,1,2>(Vec4f(x));  // rotate elements
    }

    // Type cast operator to convert to Vec3f used in geometrics:
    operator Vec3f() const {
        return Vec3f(permute4<1,2,3,0>(Vec4f(xmm)));  // rotate elements
    }
#endif // VECTOR3D_H 
};


/*****************************************************************************
*
*          Operators for Quaternion1f
*
*****************************************************************************/

// operator + : add
static inline Quaternion1f operator + (Quaternion1f const a, Quaternion1f const b) {
    return Quaternion1f(a.to_vector() + b.to_vector());
}

// operator += : add
static inline Quaternion1f & operator += (Quaternion1f & a, Quaternion1f const b) {
    a = a + b;
    return a;
}

// operator - : subtract
static inline Quaternion1f operator - (Quaternion1f const a, Quaternion1f const b) {
    return Quaternion1f(a.to_vector() - b.to_vector());
}

// operator - : unary minus
static inline Quaternion1f operator - (Quaternion1f const a) {
    return Quaternion1f(- a.to_vector());
}

// operator -= : subtract
static inline Quaternion1f & operator -= (Quaternion1f & a, Quaternion1f const b) {
    a = a - b;
    return a;
}

// operator * : quaternion multiply
static inline Quaternion1f operator * (Quaternion1f const a, Quaternion1f const b) {
    __m128 a1123 = _mm_shuffle_ps(a,a,0xE5);
    __m128 a2231 = _mm_shuffle_ps(a,a,0x7A);
    __m128 b1000 = _mm_shuffle_ps(b,b,0x01);
    __m128 b2312 = _mm_shuffle_ps(b,b,0x9E);
    __m128 t1    = _mm_mul_ps(a1123, b1000);
    __m128 t2    = _mm_mul_ps(a2231, b2312);
    __m128 t12   = _mm_add_ps(t1, t2);
    __m128 t12m  = change_sign<1,0,0,0>(Vec4f(t12));
    __m128 a3312 = _mm_shuffle_ps(a,a,0x9F);
    __m128 b3231 = _mm_shuffle_ps(b,b,0x7B);
    __m128 a0000 = _mm_shuffle_ps(a,a,0x00);
    __m128 t3    = _mm_mul_ps(a3312, b3231);
    __m128 t0    = _mm_mul_ps(a0000, b);
    __m128 t03   = _mm_sub_ps(t0, t3);
    return         _mm_add_ps(t03, t12m);
}

// operator *= : multiply
static inline Quaternion1f & operator *= (Quaternion1f & a, Quaternion1f const b) {
    a = a * b;
    return a;
}

// operator ~ : complex conjugate
// ~(a + b*i + c*j + d*k) = (a - b*i - c*j - d*k)
static inline Quaternion1f operator ~ (Quaternion1f const a) {
    return Quaternion1f(change_sign<0,1,1,1>(a.to_vector()));
}

// function reciprocal: multiplicative inverse
static inline Quaternion1f reciprocal (Quaternion1f const a) {
    Vec4f sq  = _mm_mul_ps(a,a);
    float nsq = horizontal_add(sq);
    return Quaternion1f((~a).to_vector() / Vec4f(nsq));
}

// operator / : quaternion divide is defined as
// a / b = a * reciprocal(b)
static inline Quaternion1f operator / (Quaternion1f const a, Quaternion1f const b) {
    return a * reciprocal(b);
}

// operator /= : divide
static inline Quaternion1f & operator /= (Quaternion1f & a, Quaternion1f const b) {
    a = a / b;
    return a;
}

// operator == : returns true if a == b
static inline bool operator == (Quaternion1f const a, Quaternion1f const b) {
    Vec4fb t1 = a.to_vector() == b.to_vector();
    return horizontal_and(t1);
}

// operator != : returns true if a != b
static inline bool operator != (Quaternion1f const a, Quaternion1f const b) {
    Vec4fb t1 = a.to_vector() != b.to_vector();
    return horizontal_or(t1);
}


/*****************************************************************************
*
*          Operators mixing Quaternion1f and float
*
*****************************************************************************/

// operator + : add
static inline Quaternion1f operator + (Quaternion1f const a, float b) {
    return _mm_add_ss(a, _mm_set_ss(b));
}

static inline Quaternion1f operator + (float a, Quaternion1f const b) {
    return b + a;
}

static inline Quaternion1f & operator += (Quaternion1f & a, float & b) {
    a = a + b;
    return a;
}

// operator - : subtract
static inline Quaternion1f operator - (Quaternion1f const a, float b) {
    return _mm_sub_ss(a, _mm_set_ss(b));
}

static inline Quaternion1f operator - (float a, Quaternion1f const b) {
    return _mm_sub_ps(_mm_set_ss(a), b);
}

static inline Quaternion1f & operator -= (Quaternion1f & a, float & b) {
    a = a - b;
    return a;
}

// operator * : multiply
static inline Quaternion1f operator * (Quaternion1f const a, float b) {
    return _mm_mul_ps(a, _mm_set1_ps(b));
}

static inline Quaternion1f operator * (float a, Quaternion1f const b) {
    return b * a;
}

static inline Quaternion1f & operator *= (Quaternion1f & a, float & b) {
    a = a * b;
    return a;
}

// operator / : divide
static inline Quaternion1f operator / (Quaternion1f const a, float b) {
    return _mm_div_ps(a, _mm_set1_ps(b));
}

static inline Quaternion1f operator / (float a, Quaternion1f const b) {
    return reciprocal(b) * a;
}

static inline Quaternion1f & operator /= (Quaternion1f & a, float b) {
    a = a / b;
    return a;
}


/*****************************************************************************
*
*          Functions for Quaternion1f
*
*****************************************************************************/

// function abs: calculate the norm
// abs(a + b*i + c*j + d*k) = sqrt(a*a + b*B + c*c + d*d)
static inline float abs(Quaternion1f const a) {
    Vec4f sq  = _mm_mul_ps(a,a);
    float nsq = horizontal_add(sq);
    return std::sqrt(nsq);
}

// function select
static inline Quaternion1f select (bool s, Quaternion1f const a, Quaternion1f const b) {
    return Quaternion1f(s ? a : b);
}



/*****************************************************************************
*
*                     Class Quaternion1d
*         One quaternion consisting of four double precision floats
*
*****************************************************************************/

class Quaternion1d {
protected:
    Vec4d y; // vector of 4 doubles
public:
    // default constructor
    Quaternion1d() {
    }
    // construct from real and imaginary parts = re + im0*i + im1*j + im2*k
    Quaternion1d(double re, double im0, double im1, double im2) {
        y = Vec4d(re, im0, im1, im2);
    }
    // construct from real, no imaginary part
    Quaternion1d(double re) {
        y = Vec4d(re, 0., 0., 0.);
    }
    // Constructor to convert from type __m256d used in intrinsics or Vec256de used in emulation
#if INSTRSET >= 7  // AVX
    Quaternion1d(__m256d const x) {
#else
    Quaternion1d(Vec256de const x) {
#endif
        y = x;
    }
    // Assignment operator to convert from type __m256d used in intrinsics or Vec256de used in emulation
#if INSTRSET >= 7  // AVX
    Quaternion1d & operator = (__m256d const x) {
#else
    Quaternion1d & operator = (Vec256de const x) {
#endif
        y = x;
        return *this;
    }
    // Constructor to convert from Vec4d
    Quaternion1d(Vec4d const x) {
        y = x;
    }
    // Type cast operator to convert to __m256d used in intrinsics or Vec256de used in emulation
#if INSTRSET >= 7  // AVX
    operator __m256d() const {
#else
    operator Vec256de() const {
#endif
        return y;
    }
    // Member function to convert to vector
    Vec4d to_vector() const {
        return y;
    }
    // Member function to load from array
    Quaternion1d & load(double const * p) {
        y.load(p);
        return *this;
    }
    // Member function to store into array
    void store(double * p) const {
        y.store(p);
    }
#ifdef COMPLEXVEC_H  // relations to complexvec1.h
    // construct from two Complex1d = a0 + a1 * j
    Quaternion1d(Complex1d const a0, Complex1d const a1) {
        y = Vec4d(Vec2d(a0), Vec2d(a1));
    }
    // Member functions to split into two Complex1d:
    // q = q.get_low() + q.get_high()*j
    Complex1d get_low() const {
        return Complex1d(y.get_low());
    }
    Complex1d get_high() const {
        return Complex1d(y.get_high());
    }
#endif
    // Member function to extract real part
    double real() const {
        return y.extract(0);
    }
    // Member function to extract imaginary parts, sets real part to 0
    Quaternion1d imag() const {
        return Quaternion1d(permute4<-1,1,2,3>(Vec4d(y)));
    }
#ifdef VECTOR3D_H
    // Constructor to convert from Vec3d used in geometrics:
    Quaternion1d(Vec3d const x) {
        y = permute4<3,0,1,2>(Vec4d(x));  // rotate elements
    }
    // Type cast operator to convert to Vec3d used in geometrics:
    operator Vec3d() const {
        return Vec3d(permute4<1,2,3,0>(y));  // rotate elements
    }
#endif // VECTOR3D_H 
};


/*****************************************************************************
*
*          Operators for Quaternion1d
*
*****************************************************************************/

// operator + : add
static inline Quaternion1d operator + (Quaternion1d const a, Quaternion1d const b) {
    return Quaternion1d(a.to_vector() + b.to_vector());
}

// operator += : add
static inline Quaternion1d & operator += (Quaternion1d & a, Quaternion1d const b) {
    a = a + b;
    return a;
}

// operator - : subtract
static inline Quaternion1d operator - (Quaternion1d const a, Quaternion1d const b) {
    return Quaternion1d(a.to_vector() - b.to_vector());
}

// operator - : unary minus
static inline Quaternion1d operator - (Quaternion1d const a) {
    return Quaternion1d(- a.to_vector());
}

// operator -= : subtract
static inline Quaternion1d & operator -= (Quaternion1d & a, Quaternion1d const b) {
    a = a - b;
    return a;
}

// operator * : quaternion multiply
static inline Quaternion1d operator * (Quaternion1d const a, Quaternion1d const b) {
    Vec4d a1123 = permute4<1,1,2,3>(a.to_vector());
    Vec4d a2231 = permute4<2,2,3,1>(a.to_vector());
    Vec4d b1000 = permute4<1,0,0,0>(b.to_vector());
    Vec4d b2312 = permute4<2,3,1,2>(b.to_vector());
    Vec4d t1    = a1123 * b1000;
    Vec4d t2    = a2231 * b2312;
    Vec4d t12   = t1 + t2;
    Vec4d t12m  = change_sign<1,0,0,0>(t12);
    Vec4d a3312 = permute4<3,3,1,2>(a.to_vector());
    Vec4d b3231 = permute4<3,2,3,1>(b.to_vector());
    Vec4d a0000 = permute4<0,0,0,0>(a.to_vector());
    Vec4d t3    = a3312 * b3231;
    Vec4d t0    = a0000 * b.to_vector();
    Vec4d t03   = t0  - t3;
    return        t03 + t12m;
}

// operator *= : multiply
static inline Quaternion1d & operator *= (Quaternion1d & a, Quaternion1d const b) {
    a = a * b;
    return a;
}

// operator ~ : complex conjugate
// ~(a + b*i + c*j + d*k) = (a - b*i - c*j - d*k)
static inline Quaternion1d operator ~ (Quaternion1d const a) {
    return Quaternion1d(change_sign<0,1,1,1>(a.to_vector()));
}

// function reciprocal: multiplicative inverse
static inline Quaternion1d reciprocal (Quaternion1d const a) {
    Vec4d sq  = a.to_vector() * a.to_vector();
    double nsq = horizontal_add(sq);
    return Quaternion1d((~a).to_vector() / Vec4d(nsq));
}

// operator / : quaternion divide is defined as
// a / b = a * reciprocal(b)
static inline Quaternion1d operator / (Quaternion1d const a, Quaternion1d const b) {
    return a * reciprocal(b);
}

// operator /= : divide
static inline Quaternion1d & operator /= (Quaternion1d & a, Quaternion1d const b) {
    a = a / b;
    return a;
}

// operator == : returns true if a == b
static inline bool operator == (Quaternion1d const a, Quaternion1d const b) {
    Vec4db t1 = a.to_vector() == b.to_vector();
    return horizontal_and(t1);
}

// operator != : returns true if a != b
static inline bool operator != (Quaternion1d const a, Quaternion1d const b) {
    Vec4db t1 = a.to_vector() != b.to_vector();
    return horizontal_or(t1);
}


/*****************************************************************************
*
*          Operators mixing Quaternion1d and double
*
*****************************************************************************/

// operator + : add
static inline Quaternion1d operator + (Quaternion1d const a, double b) {
    return a + Quaternion1d(b);
}

static inline Quaternion1d operator + (double a, Quaternion1d const b) {
    return b + a;
}

static inline Quaternion1d & operator += (Quaternion1d & a, double & b) {
    a = a + b;
    return a;
}

// operator - : subtract
static inline Quaternion1d operator - (Quaternion1d const a, double b) {
    return a - Quaternion1d(b);
}

static inline Quaternion1d operator - (double a, Quaternion1d const b) {
    return Quaternion1d(a) - b;
}

static inline Quaternion1d & operator -= (Quaternion1d & a, double & b) {
    a = a - b;
    return a;
}

// operator * : multiply
static inline Quaternion1d operator * (Quaternion1d const a, double b) {
    return Quaternion1d(a.to_vector() * b);
}

static inline Quaternion1d operator * (double a, Quaternion1d const b) {
    return b * a;
}

static inline Quaternion1d & operator *= (Quaternion1d & a, double & b) {
    a = a * b;
    return a;
}

// operator / : divide
static inline Quaternion1d operator / (Quaternion1d const a, double b) {
    return Quaternion1d(a.to_vector() / Vec4d(b));
}

static inline Quaternion1d operator / (double a, Quaternion1d const b) {
    return reciprocal(b) * a;
}

static inline Quaternion1d & operator /= (Quaternion1d & a, double b) {
    a = a / b;
    return a;
}


/*****************************************************************************
*
*          Functions for Quaternion1d
*
*****************************************************************************/

// function abs: calculate the norm
// abs(a + b*i + c*j + d*k) = sqrt(a*a + b*B + c*c + d*d)
static inline double abs(Quaternion1d const a) {
    Vec4d sq  = a.to_vector() * a.to_vector();
    double nsq = horizontal_add(sq);
    return std::sqrt(nsq);
}

// function select
static inline Quaternion1d select (bool s, Quaternion1d const a, Quaternion1d const b) {
    return Quaternion1d(s ? a : b);
}


#ifdef VCL_NAMESPACE
}
#endif

#endif  // QUATERNION_H
