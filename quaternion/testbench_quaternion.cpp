/*************************  testbench_quaternion.cpp   ************************
* Author:        Agner Fog
* Date created:  2019-07-13
* Last modified: 2019-07-13
* Version:       2.00
* Project:       Testbench for quaternion.h using vector class library
* Description:
* Compile and run this program to test operators and functions in quaternion package
* This file contains test cases for general operators and functions.
* Each function or operator is tested with many different combinations of input data.
*
* Instructions:
* The following parameters must be defined on the command line or added in the
* top of this file:
*
* vtype: Vector type to test
* testcase: A number defining a function or operator to test. See the cases in this file.
* seed:  Seed for random number generator. May be any integer
*
* Compile with any compiler supported by VCL.
* Specify the desired instruction set and optimization options as parameters
* to the compiler.
*
* (c) Copyright 2019 Agner Fog.
* Apache license 2.0
******************************************************************************

Test cases:
1:   operator +
2:   operator -
3:   operator *
4:   operator /
5:   vector + real
6:   vector - real
7:   real - vector
8:   vector * real
9:   vector / real
10:  real / vector
11:  reciprocal
12:  operator ~: complex conjugate
13:  abs
14:  operator ==
15:  operator !=
16:  select
20:  constructor from real scalar
21:  constructor from one real and three imaginary parts
22:  constructor from two complex numbers Complex1f
23:  constructor from two complex numbers Complex1d
25:  real(): Get real part as scalar
26:  imag(): Get imaginary part as a quaternion with the real part set to zero
27:  get_low():  Get first complex number
28:  get_high(): Get second complex number

*****************************************************************************/

#include <stdio.h>
#include <cmath>
#if defined (__linux__) && !defined(__LP64__)
#include <fpu_control.h>      // set floating point control word
#endif

#define MAX_VECTOR_SIZE 512

#ifndef INSTRSET
#define INSTRSET   7
#endif

#include <vectorclass.h>

#include "../complex/complexvec1.h"    // complex number vectors
#include "quaternion.h"                // quaternions


#ifndef testcase 
// ---------------------------------------------------------------------------
//            Specify input parameters here if running from an IDE
// ----------------------------------------------------------------------------

#define testcase 12

#define vtype Quaternion1d 

#define seed 0


#endif  // testcase

// ----------------------------------------------------------------------------
//             Declarations
// ----------------------------------------------------------------------------

// dummy vectors used for getting vector types and element type
vtype dummyc;                                    // complex vector type
typedef decltype(dummyc.to_vector()) wtype;      // corresponding normal vector type
wtype dummyv;
typedef decltype(dummyv[0]) ST;                  // scalar type
ST a0, a1;                                       // scalar operands
const int maxvectorsize = 16;                    // max number of elements in a vector
//ST oplist[maxvectorsize];                      // operand vector
int jj0;                                         // copy of vector index


float signbit_(float x) {
    union { float f; uint32_t i; } u;
    u.f = x; 
    return u.i >> 31 != 0 ? -1.f : 1.f;
}
double signbit_(double x) {
    union { double f; uint64_t i; } u;
    u.f = x; 
    return u.i >> 63 != 0 ? -1. : 1.;
}

#if defined(__FMA__) || INSTRSET >= 8            // FMA3 or AVX2
// define FMA for reference functions
float fmadd_scalar(float a, float b, float c) {
    return _mm_cvtss_f32(_mm_fmadd_ss(_mm_load_ss(&a), _mm_load_ss(&b), _mm_load_ss(&c)));
}

double fmadd_scalar(double a, double b, double c) {
    return _mm_cvtsd_f64(_mm_fmadd_sd(_mm_load_sd(&a), _mm_load_sd(&b), _mm_load_sd(&c)));
}
#endif


/************************************************************************
*
*                          Test cases
*
************************************************************************/

#if   testcase == 1    // +
inline vtype testFunction(vtype const& a, vtype const& b) { 
    return a + b; 
}
vtype referenceFunction(vtype a, vtype b) {
    ST aa[4], bb[4], cc[4];
    a.store(aa);  b.store(bb);
    for (int i=0; i<4; i++) cc[i] = aa[i] + bb[i];
    return vtype().load(cc);
}

#elif testcase == 2    // - 
inline vtype testFunction(vtype const& a, vtype const& b) { 
    return a - b; 
}
vtype referenceFunction(vtype a, vtype b) {
    ST aa[4], bb[4], cc[4];
    a.store(aa);  b.store(bb);
    for (int i=0; i<4; i++) cc[i] = aa[i] - bb[i];
    return vtype().load(cc);
}

#elif testcase == 3    // * 
inline vtype testFunction(vtype const& a, vtype const& b) {
    return a * b;
}
vtype referenceFunction(vtype a, vtype b) {
    ST aa[4], bb[4], cc[4];
    a.store(aa);  b.store(bb);
    cc[0] = aa[0]*bb[0]-aa[1]*bb[1]-aa[2]*bb[2]-aa[3]*bb[3];
    cc[1] = aa[0]*bb[1]+aa[1]*bb[0]+aa[2]*bb[3]-aa[3]*bb[2];
    cc[2] = aa[0]*bb[2]-aa[1]*bb[3]+aa[2]*bb[0]+aa[3]*bb[1];
    cc[3] = aa[0]*bb[3]+aa[1]*bb[2]-aa[2]*bb[1]+aa[3]*bb[0];
    return vtype().load(cc);
}
#define FACCURACY 100000

#elif testcase == 4    // /
inline vtype testFunction(vtype const& a, vtype const& b) {
    return a / b;
}
vtype referenceFunction(vtype a, vtype b) {
    ST aa[4], bb[4], cc[4];
    a.store(aa);  b.store(bb);
    ST bnorm2 = 0;  // calculate square norm
    for (int i=0; i<4; i++) {
        bnorm2 += bb[i]*bb[i];
    }
    for (int i=1; i<4; i++) {  // conjugate b
        bb[i] = -bb[i];
    } 
    cc[0] = aa[0]*bb[0]-aa[1]*bb[1]-aa[2]*bb[2]-aa[3]*bb[3];
    cc[1] = aa[0]*bb[1]+aa[1]*bb[0]+aa[2]*bb[3]-aa[3]*bb[2];
    cc[2] = aa[0]*bb[2]-aa[1]*bb[3]+aa[2]*bb[0]+aa[3]*bb[1];
    cc[3] = aa[0]*bb[3]+aa[1]*bb[2]-aa[2]*bb[1]+aa[3]*bb[0];
    for (int i=0; i<4; i++) {
        cc[i] /= bnorm2;
    }
    return vtype().load(cc);
}
#define FACCURACY 10000

#elif testcase == 5    // vector + real

inline vtype testFunction(vtype const& a, vtype const& b) { 
    ST b0 = b.to_vector()[0];
    return a + b0;
}
vtype referenceFunction(vtype a, vtype b) {
    ST aa[4], bb[4];
    a.store(aa);  b.store(bb);
    aa[0] += bb[0];
    return vtype().load(aa);
}

#elif testcase == 6    // vector - real
inline vtype testFunction(vtype const& a, vtype const& b) { 
    ST b0 = b.to_vector()[0];
    return a - b0;
}
vtype referenceFunction(vtype a, vtype b) {
    ST aa[4], bb[4];
    a.store(aa);  b.store(bb);
    aa[0] -= bb[0];
    return vtype().load(aa);
}

#elif testcase == 7    // real - vector
inline vtype testFunction(vtype const& a, vtype const& b) { 
    a0 = a.to_vector()[0];
    return a0 - b; 
}
vtype referenceFunction(vtype a, vtype b) {
    ST aa[4], bb[4];
    a.store(aa);  b.store(bb);
    for (int i=0; i<4; i++) bb[i] = -bb[i];
    bb[0] += aa[0];
    return vtype().load(bb);
}

#elif testcase == 8    // vector * real
inline vtype testFunction(vtype const& a, vtype const& b) { 
    ST b0 = b.to_vector()[0];
    return a * b0;
}
vtype referenceFunction(vtype a, vtype b) {
    ST aa[4], bb[4];
    a.store(aa);  b.store(bb);
    for (int i=0; i<4; i++) aa[i] *= bb[0];
    return vtype().load(aa);
}

#elif testcase == 9    // vector / real
inline vtype testFunction(vtype const& a, vtype const& b) { 
    ST b0 = b.to_vector()[0];
    return a / b0; 
}
vtype referenceFunction(vtype a, vtype b) {
    ST aa[4], bb[4];
    a.store(aa);  b.store(bb);
    for (int i=0; i<4; i++) aa[i] /= bb[0];
    return vtype().load(aa);
}

#elif testcase == 10   // real / vector
inline vtype testFunction(vtype const& a, vtype const& b) { 
    a0 = a.to_vector()[0];
    return a0 / b; 
}
vtype referenceFunction(vtype a, vtype b) { 
    ST aa[4], bb[4];
    a.store(aa);  b.store(bb);
    ST bnorm2 = 0;  // calculate square norm
    for (int i=0; i<4; i++) {
        bnorm2 += bb[i]*bb[i];
    }
    for (int i=1; i<4; i++) {  // conjugate b
        bb[i] = -bb[i];
    }
    for (int i=0; i<4; i++) {
        bb[i] = bb[i] * aa[0] / bnorm2;
    }
    return vtype().load(bb);
}
#define FACCURACY 20   // loss of precision and difference with reference function

#elif testcase == 11    // reciprocal
inline vtype testFunction(vtype const& a, vtype const& b) { 
    return reciprocal(a); 
}
vtype referenceFunction(vtype a, vtype b) { 
    ST aa[4];
    a.store(aa);
    ST anorm2 = 0;  // calculate square norm
    for (int i=0; i<4; i++) {
        anorm2 += aa[i]*aa[i];
    }
    for (int i=1; i<4; i++) {  // conjugate a
        aa[i] = -aa[i];
    }
    for (int i=0; i<4; i++) {
        aa[i] = aa[i] / anorm2;
    }
    return vtype().load(aa);
}
#define FACCURACY 20   // loss of precision and difference with reference function


#elif testcase == 12    // ~ complex conjugate
inline vtype testFunction(vtype const& a, vtype const& b) { 
    return ~a; 
}
vtype referenceFunction(vtype a, vtype b) {
    ST aa[4];
    a.store(aa);
    for (int i=1; i<4; i++) {  // conjugate 1
        aa[i] = -aa[i];
    }
    return vtype().load(aa);
}

#elif testcase == 13   // abs
inline vtype testFunction(vtype const& a, vtype const& b) { 
    return abs(a); 
}
vtype referenceFunction(vtype a, vtype b) {
    ST aa[4];  ST cc[4] = {0,0,0,0};
    a.store(aa);
    ST sum2 = 0;
    for (int i=0; i<4; i++) {  // conjugate 1
        sum2 += aa[i] * aa[i];
    }
    cc[0] = std::sqrt(sum2);
    return vtype().load(cc);
}
#define FACCURACY 4

#elif testcase == 14   // operator ==
inline vtype testFunction(vtype const& a, vtype const& b) { 
    bool eq = a == b;
    return vtype(ST(eq));
}
vtype referenceFunction(vtype a, vtype b) {
    ST aa[4], bb[4];
    a.store(aa);  b.store(bb);
    bool eq = true;
    for (int i = 0; i < 4; i++) {
        if (aa[i] != bb[i]) eq = false;
    }
    return vtype(ST(eq));
}

#elif testcase == 15   // operator !=
inline vtype testFunction(vtype const& a, vtype const& b) { 
    bool eq = a != b;
    return vtype(ST(eq));
}
vtype referenceFunction(vtype a, vtype b) {
    ST aa[4], bb[4];
    a.store(aa);  b.store(bb);
    bool eq = true;
    for (int i = 0; i < 4; i++) {
        if (aa[i] != bb[i]) eq = false;
    }
    return vtype(ST(!eq));
}


#elif testcase == 16   // select
bool sel;
inline vtype testFunction(vtype const& a, vtype const& b) { 
    sel = jj0 % 3; // arbitrary boolean value
    return select(sel, a, b);
}
vtype referenceFunction(vtype a, vtype b) {
    ST aa[4], bb[4], cc[4];
    a.store(aa);  b.store(bb);
    for (int i = 0; i < 4; i++) {
        cc[i] = sel ? aa[i] : bb[i];
    }
    return vtype().load(cc);
} 

#elif testcase == 20    // constructor from real scalar
inline vtype testFunction(vtype const& a, vtype const& b) { 
    a0 = a.to_vector()[0];
    return vtype(a0); 
}
vtype referenceFunction(vtype a, vtype b) {
    ST cc[4] = {a0,0,0,0};
    return vtype().load(cc);
}

#elif testcase == 21    // constructor from one real and three imaginary parts
inline vtype testFunction(vtype const& a, vtype const& b) { 
    ST aa[4];
    a.store(aa);
    return vtype(aa[0], aa[1], aa[2], aa[3]); 
}
vtype referenceFunction(vtype a, vtype b) {
    return a;
}

#elif testcase == 22    // constructor from two complex numbers, single precision
#ifdef COMPLEXVEC_H
inline vtype testFunction(vtype const& a, vtype const& b) { 
    ST aa[4];
    a.store(aa);
    return vtype(Complex1f(aa[0], aa[1]), Complex1f(aa[2], aa[3])); 
} 
vtype referenceFunction(vtype a, vtype b) {
    return a;
}
#endif

#elif testcase == 23    // constructor from two complex numbers, double precision
#ifdef COMPLEXVEC_H
inline vtype testFunction(vtype const& a, vtype const& b) { 
    ST aa[4];
    a.store(aa);
    return vtype(Complex1d(aa[0], aa[1]), Complex1d(aa[2], aa[3])); 
} 
vtype referenceFunction(vtype a, vtype b) {
    return a;
}
#endif

#elif testcase == 25    // real(): Get real part of complex scalar
inline vtype testFunction(vtype const& a, vtype const& b) { 
    ST r = a.real();
    return vtype(r); 
}
vtype referenceFunction(vtype a, vtype b) {
    ST aa[4];  a.store(aa);
    for (int i=1; i<4; i++) aa[i] = 0;
    return vtype().load(aa);
}

#elif testcase == 26    // imag(): Get imaginary part of complex scalar
inline vtype testFunction(vtype const& a, vtype const& b) { 
    vtype im = a.imag();
    return im;
}
vtype referenceFunction(vtype a, vtype b) {
    ST aa[4];  a.store(aa);
    aa[0] = 0;
    return vtype().load(aa);
}

#elif testcase == 27    // get_low(): Get first half-size vector
inline vtype testFunction(vtype const& a, vtype const& b) {
    auto l = a.get_low();
    return vtype(l, l); 
}
vtype referenceFunction(vtype a, vtype b) {
    ST aa[4];  a.store(aa);
    aa[2] = aa[0];
    aa[3] = aa[1];
    return vtype().load(aa);
}

#elif testcase == 28    // get_high(): Get second half-size vector
inline vtype testFunction(vtype const& a, vtype const& b) {
    auto h = a.get_high();
    return vtype(h, h); 
}
vtype referenceFunction(vtype a, vtype b) {
    ST aa[4];  a.store(aa);
    aa[0] = aa[2];
    aa[1] = aa[3];
    return vtype().load(aa);
}


#else
// End of test cases
#error unknown test case
#endif


// ----------------------------------------------------------------------------
//                         Overhead functions
// ----------------------------------------------------------------------------

const int maxerrors = 10;      // maximum errors to report
int numerr = 0;                // count errors

// type-specific load function
template <typename T, typename E>
inline void loadData(T & x, E const* p) {
    x.load(p);
}

template <typename T>
inline void loadData(T & x, bool const* p) {
    for (int i = 0; i < x.size(); i++) {
        x.insert(i, p[i]);  // bool vectors have no load function
    }
}

// type-specific printing functions

void printVal(float x) {
    printf("%10.7G", x);
}

void printVal(double x) {
    printf("%10.7G", x);
}

void printVal(bool x) {
    printf("%i", (int)x);
}

// Random number generator
class ranGen {
    // parameters for multiply-with-carry generator
    uint64_t x, carry;
public:
    ranGen(int Seed) {  // constructor
        x = Seed;  carry = 1765;  //initialize with seed
        next();  next();
    }
    uint32_t next() {  // get next random number, using multiply-with-carry method
        const uint32_t fac = 3947008974u;
        x = x * fac + carry;
        carry = x >> 32;
        x = uint32_t(x);
        return uint32_t(x);
    }
};

template <typename T>  // get random number of type T
T get_random(ranGen & rangen) {
    return (T)rangen.next();
}

template <>  // special case uint64_t
uint64_t get_random<uint64_t>(ranGen & rangen) {
    uint64_t xx;
    xx = (uint64_t)rangen.next() << 32;
    xx |= rangen.next();
    return xx;
}

template <>  // special case int64_t
int64_t get_random<int64_t>(ranGen & rangen) {
    return (int64_t)get_random<uint64_t>(rangen);
}

template <>  // special case float
float get_random<float>(ranGen & rangen) {
    union Uif {
        uint32_t i;
        float f;
    };
    Uif u1, u2;
    uint32_t r = rangen.next();                         // get 32 random bits
    // Insert exponent and random mantissa to get random number in the interval 1 <= x < 2
    // Subtract 1.0 if next bit is 0, or 1.0 - 2^-24 = 0.99999994f if next bit is 1
    u1.i = 0x3F800000 - ((r >> 8) & 1);          // bit 8
    u2.i = (r >> 9) | 0x3F800000;                // bit 9 - 31
    return u2.f - u1.f;
}

template <>  // special case float
double get_random<double>(ranGen & rangen) {
    union Uqd {
        uint64_t q;
        double d;
    };
    Uqd u1;
    uint64_t r = get_random<uint64_t>(rangen);         // get 64 random bits
    // Insert exponent and random mantissa to get random number in the interval 1 <= x < 2,
    // then subtract 1.0 to get the interval 0 <= x < 1.
    u1.q = (r >> 12) | 0x3FF0000000000000;       // bit 12 - 63
    return u1.d - 1.0;
}
template <>  // special case bool
bool get_random<bool>(ranGen & rangen) {
    return (rangen.next() & 1) != 0;
}


// make random number generator instance
ranGen ran(seed);

// bit_cast function to make special values
float bit_castf(uint32_t x) {  // uint64_t -> double
    union {
        uint32_t i;
        float f;
    } u;
    u.i = x;
    return u.f;
}

double bit_castd(uint64_t x) {  // uint32_t -> float
    union {
        uint64_t i;
        double f;
    } u;
    u.i = x;
    return u.f;
}


// template to generate list of testdata
template <typename T>
class TestData {
public:
    enum LS {
        // define array size. Must be a multiple of vector size:
        listsize = 1024
    };
    TestData() {                            // constructor
        int i;                              // loop counter
        if (T(1.1f) != 1) {
            // floating point type
            // fill boundary data into array
            for (i = 0; i < 20; i++) {
                list[i] = T((i - 4) * T(0.25));
            }
#ifdef TESTNAN   // test also with NAN, INF, and other special data
            // additional special values, float:
            if constexpr (sizeof(ST) == 4) {
                list[i++] = (T)bit_castf(0x80000000);   // -0
                list[i++] = (T)bit_castf(0x00800000);   // smallest positive normal number
                list[i++] = (T)bit_castf(0x80800000);   // largest negative normal number
                list[i++] = (T)bit_castf(0x3F7FFFFF);   // nextafter 1.0, 0
                list[i++] = (T)bit_castf(0x3F800001);   // nextafter 1.0, 2
                list[i++] = (T)bit_castf(0x7F800000);   // inf
                list[i++] = (T)bit_castf(0xFF800000);   // -inf
                list[i++] = (T)bit_castf(0x7FF00000);   // nan
            }
            else { // double
                list[i++] = (T)bit_castd(0x8000000000000000);   // -0
                list[i++] = (T)bit_castd(0x0010000000000000);   // smallest positive normal number
                list[i++] = (T)bit_castd(0x8010000000000000);   // largest negative normal number
                list[i++] = (T)bit_castd(0x3FEFFFFFFFFFFFFF);   // nextafter 1.0, 0
                list[i++] = (T)bit_castd(0x3FF0000000000001);   // nextafter 1.0, 2
                list[i++] = (T)bit_castd(0x7FF0000000000000);   // inf
                list[i++] = (T)bit_castd(0xFFF0000000000000);   // -inf
                list[i++] = (T)bit_castd(0x7FFC000000000000);   // nan
            }
#endif
            // fill random data into rest of array
            for (; i < listsize; i++) {
                list[i] = get_random<T>(ran) * (T)100;
            }
        }
        else {
            // integer type
            // fill boundary data into array
            for (i = 0; i < 6; i++) {
                list[i] = T(i - 2);
            }
            // data near mid-point of unsigned integers, or overflow point of signed integers:
            uint64_t m = (uint64_t(1) << (sizeof(T) * 8 - 1)) - 2;
            for (; i < 11; i++) {
                list[i] = T(m++);
            }
            // fill random data into rest of array
            for (; i < listsize; i++) {
                list[i] = get_random<T>(ran);
            }
        }
    }
    T list[listsize];                       // array of test data
    int size() {                            // get list size
        return listsize;
    }
};


// get value of least significant bit
float delta_unit(float x) {
    union {
        float f;
        uint32_t i;
    } u;
    x = fabsf(x);
    Vec4f xv = Vec4f(x);
    if (!(is_finite(xv)[0])) return 1.f;
    if (x == 0.f || is_subnormal(xv)[0]) {
        u.i = 0x00800000;  // smallest positive normal number
        return u.f;
    }
    float x1 = x;
    u.f = x;
    u.i++;
    return u.f - x1;
}

double delta_unit(double x) {
    union {
        double f;
        uint64_t i;
    } u;
    x = fabs(x);
    Vec2d xv = Vec2d(x);
    if (!(is_finite(xv)[0])) return 1.;
    if (x == 0. || is_subnormal(xv)[0]) {
        u.i = 0x0010000000000000;  // smallest positive normal number
        return u.f;
    }
    double x1 = x;
    u.f = x;
    u.i++;
    return u.f - x1;
}


// compare two scalars. return true if different
template <typename T>
inline bool compare_scalars(T const a, T const b) {
    return a == b;
}

// special cases for float and double:
template <>
inline bool compare_scalars<float>(float const a, float const b) {
    if (a == b || (a != a && b != b)) return true; // return false if equal or both are NAN
#ifdef FACCURACY     // accept minor difference
    float dif = fabsf(a - b) / delta_unit(a);
    if (dif <= FACCURACY) return true;
    printf("\n%.0f ULP ", dif);
#endif
    return false;
}

template <>
inline bool compare_scalars<double>(double const a, double const b) {
    if (a == b || (a != a && b != b)) return true; // return false if equal or both are NAN
#ifdef FACCURACY     // accept minor difference
    double dif = fabs(a - b) / delta_unit(a);
    if (dif <= FACCURACY) return true;
    printf("\n%.0f ULP ", dif);
#endif
    return false;
}

// compare two vectors. return true if different
template <typename T>
inline bool compare_vectors(T const& a, T const& b) {
    {
        for (int i = 0; i < a.size(); i++) {
            if (!compare_scalars(a[i], b[i])) return false;
        }
    }
    return true;
}

#ifndef FACCURACY
#define FACCURACY 1
#endif

// compare two complex numbers. return true if different
inline ST compare_quad(vtype const& a, vtype const& b) {
    ST alist[4], blist[4];
    a.store(alist);  b.store(blist);
    ST dif, dif0 = 0;
    for (int i = 0; i < 4; i++) {
        ST r = fabs(blist[i]);
        if (r < 1.E-2) r = 1;  // use relative error for results near zero
        dif = ST(fabs(alist[i] - blist[i]) / delta_unit(r));
        if (dif > dif0) dif0 = dif;
    }
    if (dif0 > FACCURACY) return dif0;
    return 0;
}

// program entry
int main() {
    //const int vectorsize = vtype::size();

#if defined (__linux__) && !defined(__LP64__)
    // Some 32-bit compilers use x87 calculations with long double precision for 
    // the reference function. This may give slightly different results because
    // the value is rounded twice. To get exactly the same value in the test function
    // and the reference function, we change the precision of x87 calculations.
    // (the fpu control function is different in Windows, but the precision is already
    // reduced in Windows anyway)
    fpu_control_t fpcw = 0x27f;
    _FPU_SETCW(fpcw);
#endif

    vtype a, b, result, ref;  // complex vectors for operands and result

    // make lists of test data
    TestData<ST> adata, bdata;

    int i, j, k = 0;   // loop counters

    for (i = 0; i < adata.size(); i += wtype::size()) {
        //a.load(adata.list + i);
        loadData(a, adata.list + i);

        for (j = 0; j < bdata.size(); j += wtype::size()) {
            loadData(b, bdata.list + j);
            jj0 = j;

            // function under test:
            result = testFunction(a, b);
            ref = referenceFunction(a, b);
            ST dif = compare_quad(result, ref);
            if (dif != 0) {
                // values are different. report error 
                if (++numerr == 1) {
                    printf("\ntest case %i:", testcase);  // print test case first time
                }
                ST alist[4], blist[4], tlist[4], rlist[4];
                a.store(alist); b.store(blist); result.store(tlist); ref.store(rlist);
                printf("\nError at %i, %i, dif = %.2G:", i, j, dif);
                for (k = 0; k < 4; k++) {
                    printf("\n%7.4G op %7.4G -> %7.4G, expected %7.4G)",
                        alist[k], blist[k], tlist[k], rlist[k]);
                }
            }
            if (numerr > maxerrors) {
                exit(1);      // stop after maxerrors
            }
        }
    }

    if (numerr == 0) {
        printf("\nsuccess\n");
    }
    printf("\n");

    return numerr;
}
