/****************************  vector3d.h   ***********************************
* Author:        Agner Fog
* Date created:  2012-08-01
* Last modified: 2023-05-14
* Version:       2.02.00
* Project:       Extension to vector class library
* Description:   Classes for 3-dimensional vectors, including operators and functions
* The following classes are defined:
* Vec3Df:        A vector of 3 single precision floats
* Vec3Dd:        A vector of 3 double precision floats
*
* (c) Copyright 2012-2023 Apache License version 2.0 or later
\*****************************************************************************/

#ifndef VECTOR3D_H
#define VECTOR3D_H  20200

#include "vectorclass.h"
#include <cmath>          // define math library functions

#if VECTORCLASS_H < 20000
#error Incompatible version of vector class library. Must use version 2 or later
#endif

#ifdef VCL_NAMESPACE
namespace VCL_NAMESPACE {
#endif

/*****************************************************************************
*
*               Class Vec3Df: vector of 3 single precision floats
*
*****************************************************************************/

class Vec3Df {
protected:
    __m128 xmm; // Float vector
public:
    // default constructor
    Vec3Df() = default;
    // construct from three coordinates
    Vec3Df(float x, float y, float z) {
        xmm = Vec4f(x, y, z, 0.f);
    }
    // Constructor to convert from Vec4f
    Vec3Df(Vec4f const x) {
        xmm = x;
        // cutoff(3);
    }
    // Constructor to convert from type __m128 used in intrinsics:
    Vec3Df(__m128 const x) {
        xmm = x;
    }
    // Assignment operator to convert from type __m128 used in intrinsics:
    Vec3Df & operator = (__m128 const x) {
        xmm = x;
        return *this;
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
    Vec3Df & load(float const * p) {
        xmm = Vec4f().load_partial(3, p);
        return *this;
    }
    // Member function to store into array
    void store(float * p) const {
        Vec4f(xmm).store_partial(3, p);
    }
    // get x part
    float get_x() const {
        return _mm_cvtss_f32(xmm);
    }
    // get y part
    float get_y() const {
        return Vec4f(xmm).extract(1);
    }
    // get z part
    float get_z() const {
        return Vec4f(xmm).extract(2);
    }
    // Member function to extract one coordinate
    float extract(int index) const {
        return Vec4f(xmm).extract(index);
    }
    // Operator [] to extract one coordinate
    // Operator [] can only read an element, not write.
    float operator [] (uint32_t index) const {
        return extract(index);
    }
    // Insert one coordinate
    Vec3Df & insert (uint32_t index, float x) {
        xmm = Vec4f(xmm).insert(index, x);
        return *this;
    }
    static constexpr int size() {
        return 1;
    }
    static constexpr int elementtype() {
        return 0x210;
    }
};

/*****************************************************************************
*
*          Operators for Vec3Df
*
*****************************************************************************/

// operator + : add
static inline Vec3Df operator + (Vec3Df const a, Vec3Df const b) {
    return Vec3Df(Vec4f(a) + Vec4f(b));
}

// operator += : add
static inline Vec3Df & operator += (Vec3Df & a, Vec3Df const b) {
    a = a + b;
    return a;
}

// operator - : subtract
static inline Vec3Df operator - (Vec3Df const a, Vec3Df const b) {
    return Vec3Df(Vec4f(a) - Vec4f(b));
}

// operator - : unary minus
static inline Vec3Df operator - (Vec3Df const a) {
    return Vec3Df(- Vec4f(a));
}

// operator -= : subtract
static inline Vec3Df & operator -= (Vec3Df & a, Vec3Df const b) {
    a = a - b;
    return a;
}

// operator * : multiply element-by-element
// (see also cross_product and dot_product)
static inline Vec3Df operator * (Vec3Df const a, Vec3Df const b) {
    return Vec3Df(Vec4f(a) * Vec4f(b));
}

// operator *= : multiply element-by-element
static inline Vec3Df & operator *= (Vec3Df & a, Vec3Df const b) {
    a = a * b;
    return a;
}

// operator / : divide element-by-element
static inline Vec3Df operator / (Vec3Df const a, Vec3Df const b) {
    return Vec3Df(Vec4f(a) / Vec4f(b));
}

// operator /= : divide element-by-element
static inline Vec3Df & operator /= (Vec3Df & a, Vec3Df const b) {
    a = a / b;
    return a;
}

// operator == : returns true if a == b
static inline bool operator == (Vec3Df const a, Vec3Df const b) {
    Vec4fb t1 = Vec4f(a) == Vec4f(b);
#if INSTRSET >= 10
    return (uint8_t(t1) & 7) == 7;
#else
    Vec4fb t2 = _mm_shuffle_ps(t1, t1, 0x24);  // ignore unused top element
    return horizontal_and(t2);
#endif
}

// operator != : returns true if a != b
static inline bool operator != (Vec3Df const a, Vec3Df const b) {
    Vec4fb t1 = Vec4f(a) != Vec4f(b);
#if INSTRSET >= 10
    return (uint8_t(t1) & 7) != 0;
#else
    Vec4fb t2 = _mm_shuffle_ps(t1, t1, 0x24);  // ignore unused top element
    return horizontal_or(t2);
#endif
}

/*****************************************************************************
*
*          Operators mixing Vec3Df and float
*
*****************************************************************************/

// operator * : multiply
static inline Vec3Df operator * (Vec3Df const a, float b) {
    return _mm_mul_ps(a, _mm_set1_ps(b));
}
static inline Vec3Df operator * (float a, Vec3Df const b) {
    return b * a;
}
static inline Vec3Df & operator *= (Vec3Df & a, float & b) {
    a = a * b;
    return a;
}

// operator / : divide
static inline Vec3Df operator / (Vec3Df const a, float b) {
    return _mm_div_ps(a, _mm_set1_ps(b));
}

static inline Vec3Df & operator /= (Vec3Df & a, float b) {
    a = a / b;
    return a;
}


/*****************************************************************************
*
*          Functions for Vec3Df
*
*****************************************************************************/

// function cross_product
static inline Vec3Df cross_product (Vec3Df const a, Vec3Df const b) {
    Vec4f a1 = permute4<1,2,0,V_DC>(Vec4f(a));
    Vec4f b1 = permute4<1,2,0,V_DC>(Vec4f(b));
    Vec4f a2 = permute4<2,0,1,V_DC>(Vec4f(a));
    Vec4f b2 = permute4<2,0,1,V_DC>(Vec4f(b));
    Vec4f c  = a1 * b2 - a2 * b1;
    return c.cutoff(3);
}

// function dot_product
static inline float dot_product (Vec3Df const a, Vec3Df const b) {
    Vec4f c = (Vec4f(a) * Vec4f(b)).cutoff(3);
    return horizontal_add(c);
}

// function vector_length
static inline float vector_length (Vec3Df const a) {
    return std::sqrt(dot_product(a,a));
}

// function normalize_vector
static inline Vec3Df normalize_vector (Vec3Df const a) {
    return a / vector_length(a);
}

// function select
static inline Vec3Df select (bool s, Vec3Df const a, Vec3Df const b) {
    return s ? a : b;
}

// function rotate
// The vector a is rotated by multiplying by the matrix defined by the three columns col0, col1, col2
static inline Vec3Df rotate (Vec3Df const col0, Vec3Df const col1, Vec3Df const col2, Vec3Df const a) {
    Vec4f xbroad = permute4<0,0,0,V_DC>(Vec4f(a));  // broadcast x
    Vec4f ybroad = permute4<1,1,1,V_DC>(Vec4f(a));  // broadcast y
    Vec4f zbroad = permute4<2,2,2,V_DC>(Vec4f(a));  // broadcast z
    Vec4f r = col0.to_vector() * xbroad + col1.to_vector() * ybroad + col2.to_vector() * zbroad;
    return r.cutoff(3);
}


/*****************************************************************************
*
*               Class Vec3Dd: vector of 3 double precision floats
*
*****************************************************************************/

class Vec3Dd  {
protected:
    Vec4d yy; // vector of 4 doubles
public:
    // default constructor
    Vec3Dd() = default;
    // construct from three coordinates
    Vec3Dd(double x, double y, double z) {
        yy = Vec4d(x, y, z, 0.);
    }
    // Constructor to convert from Vec4d
    Vec3Dd(Vec4d const x) {
        yy = x;
        // cutoff(3);
    }
    // Constructor to convert from type __m256d used in intrinsics or Vec256de used in emulation
#if INSTRSET >= 7  // AVX
    Vec3Dd(__m256d const x) {
        yy = x;
    }
#else
    Vec3Dd(Vec256de const x) {
        yy = x;
    }
#endif
    // Assignment operator to convert from type __m256d used in intrinsics or Vec256de used in emulation
#if INSTRSET >= 7  // AVX
    Vec3Dd & operator = (__m256d const x) {
#else
    Vec3Dd & operator = (Vec256de const x) {
#endif
        yy = x;
        return *this;
    }
    // Type cast operator to convert to __m256d used in intrinsics or Vec256de used in emulation
#if INSTRSET >= 7  // AVX
    operator __m256d() const {
        return yy;
    }
#endif
    // Member function to load from array
    Vec3Dd & load(double const * p) {
        yy.load_partial(3, p);
        return *this;
    }
    // Member function to store into array
    void store(double * p) const {
        yy.store_partial(3, p);
    }
    // Member function to convert to vector
    Vec4d to_vector() const {
        return yy;
    }
    // get x part
    double get_x() const {
        return _mm_cvtsd_f64(yy.get_low());
    }
    // get y part
    double get_y() const {
        return yy.extract(1);
    }
    // get z part
    double get_z() const {
        return yy.extract(2);
    }
    // Member function to extract one coordinate
    double extract(uint32_t index) const {
        return yy.extract(index);
    }
    // Operator [] to extract one coordinate
    // Operator [] can only read an element, not write.
    double operator [] (uint32_t index) const {
        return extract(index);
    }
    // Insert one coordinate
    Vec3Dd & insert (uint32_t index, double x) {
        yy.insert(index, x);
        return *this;
    }
    static constexpr int size() {
        return 1;
    }
    static constexpr int elementtype() {
        return 0x211;
    }
};

/*****************************************************************************
*
*          Operators for Vec3Dd
*
*****************************************************************************/

// operator + : add
static inline Vec3Dd operator + (Vec3Dd const a, Vec3Dd const b) {
    return Vec3Dd(a.to_vector() + b.to_vector());
}

// operator += : add
static inline Vec3Dd & operator += (Vec3Dd & a, Vec3Dd const b) {
    a = a + b;
    return a;
}

// operator - : subtract
static inline Vec3Dd operator - (Vec3Dd const a, Vec3Dd const b) {
    return Vec3Dd(a.to_vector() - b.to_vector());
}

// operator - : unary minus
static inline Vec3Dd operator - (Vec3Dd const a) {
    return Vec3Dd(- a.to_vector());
}

// operator -= : subtract
static inline Vec3Dd & operator -= (Vec3Dd & a, Vec3Dd const b) {
    a = a - b;
    return a;
}

// operator * : multiply element-by-element
// (see also cross_product and dot_product)
static inline Vec3Dd operator * (Vec3Dd const a, Vec3Dd const b) {
    return Vec3Dd(a.to_vector() * b.to_vector());
}

// operator *= : multiply element-by-element
static inline Vec3Dd & operator *= (Vec3Dd & a, Vec3Dd const b) {
    a = a * b;
    return a;
}

// operator / : divide element-by-element
static inline Vec3Dd operator / (Vec3Dd const a, Vec3Dd const b) {
    return Vec3Dd(a.to_vector() / b.to_vector());
}

// operator /= : divide element-by-element
static inline Vec3Dd & operator /= (Vec3Dd & a, Vec3Dd const b) {
    a = a / b;
    return a;
}

// operator == : returns true if a == b
static inline bool operator == (Vec3Dd const a, Vec3Dd const b) {
    Vec4db t1 = a.to_vector() == b.to_vector();
#if INSTRSET >= 10
    return (uint8_t(t1) & 7) == 7;
#elif INSTRSET >= 7  // AVX
    Vec4db t2 = Vec4db(permute4<0,1,2,2>(Vec4d(t1))); // ignore unused top element
    return horizontal_and(t2);
#else
    Vec2db u0 = t1.get_low();
    Vec2db u1 = t1.get_high();
    u1 = permute2<0,0>(Vec2d(u1));                    // ignore unused top element
    return horizontal_and(u0 & u1);
#endif
}

// operator != : returns true if a != b
static inline bool operator != (Vec3Dd const a, Vec3Dd const b) {
    Vec4db t1 = a.to_vector() != b.to_vector();
#if INSTRSET >= 10
    return (uint8_t(t1) & 7) != 0;
#elif INSTRSET >= 7  // AVX
    Vec4db t2 = Vec4db(permute4<0,1,2,2>(Vec4d(t1))); // ignore unused top element
    return horizontal_and(t2);
#else
    Vec2db u0 = t1.get_low();
    Vec2db u1 = t1.get_high();
    u1 = permute2<0,0>(Vec2d(u1));                    // ignore unused top element
    return horizontal_or(u0 | u1);
#endif
}

/*****************************************************************************
*
*          Operators mixing Vec3Dd and double
*
*****************************************************************************/

// operator * : multiply
static inline Vec3Dd operator * (Vec3Dd const a, double b) {
    return a.to_vector() * Vec4d(b);
}
static inline Vec3Dd operator * (double a, Vec3Dd const b) {
    return b * a;
}
static inline Vec3Dd & operator *= (Vec3Dd & a, double & b) {
    a = a * b;
    return a;
}

// operator / : divide
static inline Vec3Dd operator / (Vec3Dd const a, double b) {
    return a.to_vector() / Vec4d(b);
}

static inline Vec3Dd & operator /= (Vec3Dd & a, double b) {
    a = a / b;
    return a;
}


/*****************************************************************************
*
*          Functions for Vec3Dd
*
*****************************************************************************/

// function cross_product
static inline Vec3Dd cross_product (Vec3Dd const a, Vec3Dd const b) {
    Vec4d a1 = permute4<1,2,0,V_DC>(a.to_vector());
    Vec4d b1 = permute4<1,2,0,V_DC>(b.to_vector());
    Vec4d a2 = permute4<2,0,1,V_DC>(a.to_vector());
    Vec4d b2 = permute4<2,0,1,V_DC>(b.to_vector());
    Vec4d c  = a1 * b2 - a2 * b1;
    return c.cutoff(3);
}

// function dot_product
static inline double dot_product (Vec3Dd const a, Vec3Dd const b) {
    Vec4d c  = (a.to_vector() * b.to_vector()).cutoff(3);
    return horizontal_add(c);
}

// function vector_length
static inline double vector_length (Vec3Dd const a) {
    return std::sqrt(dot_product(a,a));
}

// function normalize_vector
static inline Vec3Dd normalize_vector (Vec3Dd const a) {
    return a / vector_length(a);
}

// function select
static inline Vec3Dd select (bool s, Vec3Dd const a, Vec3Dd const b) {
    return s ? a : b;
}

// function rotate
// The vector a is rotated by multiplying by the matrix defined by the three columns col0, col1, col2
static inline Vec3Dd rotate (Vec3Dd const col0, Vec3Dd const col1, Vec3Dd const col2, Vec3Dd const a) {
    Vec3Dd xbroad = permute4<0,0,0,V_DC>(a.to_vector());  // broadcast x
    Vec3Dd ybroad = permute4<1,1,1,V_DC>(a.to_vector());  // broadcast y
    Vec3Dd zbroad = permute4<2,2,2,V_DC>(a.to_vector());  // broadcast z
    Vec3Dd r = col0 * xbroad + col1 * ybroad + col2 * zbroad;
    return r.to_vector().cutoff(3);
}


/*****************************************************************************
*
*          Conversion functions
*
*****************************************************************************/

// function to_single: convert Vec3Dd to Vec3Df
static inline Vec3Df to_float(Vec3Dd const a) {
#if INSTRSET >= 7  // AVX
    return _mm256_cvtpd_ps(a);
#else
    //return Vec3Df(Vec4f(compress(a.to_vector().get_low(), a.to_vector().get_high())));
    return to_float(a.to_vector());
#endif
}

// function to_double: convert Vec3Df to Vec3Dd
static inline Vec3Dd to_double(Vec3Df const a) {
#if INSTRSET >= 7  // AVX
    return _mm256_cvtps_pd(a);
#else
    return to_double(a.to_vector());
#endif
}

#ifdef VCL_NAMESPACE
}
#endif

#endif  // VECTOR3D_H
