/*************************  testbench_complex.cpp   ***************************
* Author:        Agner Fog
* Date created:  2019-07-10
* Last modified: 2019-07-10
* Version:       2.00
* Project:       Testbench for complex number vectors in complexvec1 using vector class library
* Description:
* Compile and run this program to test operators and functions in complex vector library
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
11:  operator ~: complex conjugate
20:  constructor from real scalar
21:  constructor from single real+imag pair
22:  constructor from multiple real+imag pairs
23:  constructor from complex scalar
24:  constructor from two half size complex vectors
25:  real(): Get real part of complex scalar
26:  imag(): Get imaginary part of complex scalar
27:  get_low(): Get first half size vector
28:  get_high(): Get second half size vector
29:  extract: Get one complex element
30:  operator ==
31:  operator !=
32:  select
33:  to_float: lower precision
34:  to_double: higher precision
40:  abs
41:  csqrt
100: cexp
101: clog

*****************************************************************************/

#include <stdio.h>
#include <cmath>
#if defined (__linux__) && !defined(__LP64__)
#include <fpu_control.h>      // set floating point control word
#endif

#define MAX_VECTOR_SIZE 512

#ifndef INSTRSET
#define INSTRSET  10
#endif
#include <vectorclass.h>
#include "vectormath_exp.h"                 // exp and log functions
#include "vectormath_trig.h"                // trigonometric functions
#include "complexvec1.h"                    // complex number vectors


#ifndef testcase 
// ---------------------------------------------------------------------------
//            Specify input parameters here if running from an IDE
// ----------------------------------------------------------------------------

#define testcase 100

#define vtype Complex4f 

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
ST oplist[maxvectorsize];                        // operand vector
int k0;                                          // copy of vector index

struct Cpair {                                   // complex number structure for reference function
    ST re;                                       // real part
    ST im;                                       // imaginary part
    Cpair(){}                                    // default constructor
    Cpair(ST r, ST i) {                          // constructor from elements
        re = r;  im = i;
    }
    Cpair & load(ST * p) {                       // load data from array
        re = *p;  im = *(p+1);
        return *this;
    }
    void store(ST * p) {                         // save data to array
        *p = re;  *(p+1) = im;
    }
};

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
Cpair referenceFunction(Cpair a, Cpair b) {
    return Cpair (a.re + b.re, a.im + b.im);
}

#elif testcase == 2    // - 
inline vtype testFunction(vtype const& a, vtype const& b) { 
    return a - b; 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return Cpair (a.re - b.re, a.im - b.im);
}

#elif testcase == 3    // * 
inline vtype testFunction(vtype const& a, vtype const& b) {
    return a * b;
//    return b * a;  // note: a*b and b*a give different precision
}
Cpair referenceFunction(Cpair a, Cpair b) {
    ST x, y;
#if defined(__FMA__) || INSTRSET >= 8            // FMA3 or AVX2
    x = fmadd_scalar(a.re, b.re, - a.im * b.im); // calculate with same precision as test function
    y = fmadd_scalar(a.re, b.im,   b.re * a.im);
#else
    x = a.re * b.re - a.im * b.im;
    y = a.re * b.im + b.re * a.im;
#endif
    return Cpair (x, y);
}

#elif testcase == 4    // /
inline vtype testFunction(vtype const& a, vtype const& b) {
    return a / b;
}
Cpair referenceFunction(Cpair a, Cpair b) {
    ST x, y;
#if defined(__FMA__) || INSTRSET >= 8              // FMA3 or AVX2
    x = fmadd_scalar(a.im, b.im ,   a.re * b.re);  // calculate with same precision as test function
    y = fmadd_scalar(a.im, b.re , - a.re * b.im);
#else
    x = a.im * b.im + a.re * b.re;
    y = a.im * b.re - a.re * b.im;
#endif
    ST ssum = b.re * b.re + b.im * b.im;
    return Cpair (x / ssum, y / ssum);
}
#define FACCURACY 1000   // loss of precision and difference with reference function


#elif testcase == 5    // vector + real

inline vtype testFunction(vtype const& a, vtype const& b) { 
    a0 = a.to_vector()[0];
    return a0 + b; 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return Cpair (a0 + b.re, b.im);
}

#elif testcase == 6    // vector - real
inline vtype testFunction(vtype const& a, vtype const& b) { 
    a0 = a.to_vector()[0];
    return b - a0; 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return Cpair (b.re - a0, b.im);
} 

#elif testcase == 7    // real - vector
inline vtype testFunction(vtype const& a, vtype const& b) { 
    a0 = a.to_vector()[0];
    return a0 - b; 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return Cpair (a0 - b.re, - b.im);
} 

#elif testcase == 8    // vector * real
inline vtype testFunction(vtype const& a, vtype const& b) { 
    a0 = a.to_vector()[0];
    return a0 * b; 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return Cpair (a0 * b.re, a0 * b.im);
}

#elif testcase == 9    // vector / real
inline vtype testFunction(vtype const& a, vtype const& b) { 
    a0 = a.to_vector()[0];
    return b / a0; 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return Cpair (b.re / a0, b.im / a0);
}

#elif testcase == 10   // real / vector
inline vtype testFunction(vtype const& a, vtype const& b) { 
    a0 = a.to_vector()[0];
    return a0 / b; 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    ST c = a0 / (b.re*b.re + b.im*b.im);
    return Cpair (c * b.re, -c * b.im);
}
#define FACCURACY 200   // loss of precision and difference with reference function

#elif testcase == 11    // ~ complex conjugate
inline vtype testFunction(vtype const& a, vtype const& b) { 
    return ~a; 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return Cpair (a.re, - a.im);
}

#elif testcase == 20    // constructor from real scalar
inline vtype testFunction(vtype const& a, vtype const& b) { 
    a0 = a.to_vector()[0];
    return vtype(a0); 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return Cpair (a0, 0);
}

#elif testcase == 21    // constructor from single real+imag pair
inline vtype testFunction(vtype const& a, vtype const& b) { 
    a0 = a.to_vector()[0];
    a1 = a.to_vector()[1];
    return vtype(a0, a1); 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return Cpair (a0, a1);
}

#elif testcase == 22    // constructor from multiple real+imag pairs
template <typename V>
inline V testFunction(V const& a, V const& b) { 
    a.store(oplist);
    if constexpr (V::size() == 1) {
        return V(oplist[0], oplist[1]);
    }
    else if constexpr (V::size() == 2) {
        return V(oplist[0], oplist[1], oplist[2], oplist[3]);
    }
    else if constexpr (V::size() == 4) {
        return V(oplist[0], oplist[1], oplist[2], oplist[3], oplist[4], oplist[5], oplist[6], oplist[7]);
    }
    else /*if constexpr (V::size() == 8)*/ {
        return V(oplist[0], oplist[1], oplist[2], oplist[3], oplist[4], oplist[5], oplist[6], oplist[7],
            oplist[8], oplist[9], oplist[10], oplist[11], oplist[12], oplist[13], oplist[14], oplist[15]);
    }
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return a;
}

#elif testcase == 23    // constructor from complex scalar
inline vtype testFunction(vtype const& a, vtype const& b) { 
    a0 = a.to_vector()[0];
    a1 = a.to_vector()[1];
    auto s = a.extract(0);   // complex scalar
    return vtype(s); 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return Cpair (a0, a1);
}

#elif testcase == 24    // constructor from two half size complex vectors
inline vtype testFunction(vtype const& a, vtype const& b) { 
    a.store(oplist);
    decltype(a.get_low()) lo, hi; // half size vectors
    lo.load(oplist);
    hi.load(oplist + 2*lo.size());
    return vtype(lo, hi); 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return a;
}

#elif testcase == 25    // real(): Get real part of complex scalar
inline vtype testFunction(vtype const& a, vtype const& b) { 
    ST r = a.real();
    return vtype(r, 0); 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return Cpair(a.re, 0);
}

#elif testcase == 26    // imag(): Get imaginary part of complex scalar
inline vtype testFunction(vtype const& a, vtype const& b) { 
    ST i = a.imag();
    return vtype(i, 0); 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return Cpair(a.im, 0);
}

#elif testcase == 27    // get_low(): Get first half-size vector
inline vtype testFunction(vtype const& a, vtype const& b) {
    a.store(oplist);
    auto l = a.get_low();
    return vtype(l, l); 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    int k2 = k0 % vtype::size();             // get lower half twice
    return Cpair(oplist[k2], oplist[k2+1]);
}

#elif testcase == 28    // get_high(): Get second half-size vector
inline vtype testFunction(vtype const& a, vtype const& b) {
    a.store(oplist);
    auto h = a.get_high();
    return vtype(h, h); 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    int k2 = (k0 % vtype::size()) + vtype::size();// get upper half twice
    return Cpair(oplist[k2], oplist[k2+1]);
}

#elif testcase == 29    // extract: Get one complex element
inline vtype testFunction(vtype const& a, vtype const& b) {
    a.store(oplist);
    a0 = b.to_vector()[0];
    uint32_t i = uint32_t(a0) % vtype::size();
    auto e = a.extract(i);
    return vtype(e); 
} 
Cpair referenceFunction(Cpair a, Cpair b) {
    uint32_t i = uint32_t(a0) % vtype::size();
    return Cpair(oplist[i*2], oplist[i*2+1]);
}

#elif testcase == 30    // operator ==
inline vtype testFunction(vtype const& a, vtype const& b) {
    auto eq = a == b;   // boolean vector
    return select(eq, vtype(1), vtype(0));
}
Cpair referenceFunction(Cpair a, Cpair b) {
    if (a.re == b.re && a.im == b.im) {
        return Cpair(1,0);
    }
    else {
        return Cpair(0,0);
    }
}

#elif testcase == 31 || testcase == 32   // operator !=, select
inline vtype testFunction(vtype const& a, vtype const& b) {
    auto eq = a != b;   // boolean vector
    return select(eq, vtype(1), vtype(0));
}
Cpair referenceFunction(Cpair a, Cpair b) {
    if (a.re != b.re || a.im != b.im) {
        return Cpair(1,0);
    }
    else {
        return Cpair(0,0);
    }
}

#elif testcase == 33    // to_float: lower precision
float fllist[maxvectorsize];

inline vtype testFunction(vtype const& a, vtype const& b) {
    to_float(a).store(fllist);
    return a;
}
Cpair referenceFunction(Cpair a, Cpair b) {
    if (float(a.re) != fllist[k0] || float(a.im) != fllist[k0+1]) return Cpair(0,0);
    else return a;
}

#elif testcase == 34    // to_double: higher precision
double dbllist[maxvectorsize];

inline vtype testFunction(vtype const& a, vtype const& b) {
    to_double(a).store(dbllist);
    return a;
}
Cpair referenceFunction(Cpair a, Cpair b) {
    if (double(a.re) != dbllist[k0] || double(a.im) != dbllist[k0+1]) return Cpair(0,0);
    else return a;
}

#elif testcase == 40    // abs
inline vtype testFunction(vtype const& a, vtype const& b) { 
    return abs(a); 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    ST ssum = a.re * a.re + a.im * a.im;
    return Cpair (std::sqrt(ssum), 0);
}

#elif testcase == 41    // csqrt
inline vtype testFunction(vtype const& a, vtype const& b) { 
    return csqrt(a); 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    ST n = (ST)sqrt(a.re * a.re + a.im * a.im);
    return Cpair ((ST)std::sqrt((n + a.re)*0.5f), signbit_(a.im) * (ST)std::sqrt((n - a.re)*0.5f) );
}
#define FACCURACY 100    // loss of precision and difference with reference function


#elif testcase == 100    // cexp. Complex exponential function
inline vtype testFunction(vtype const& a, vtype const& b) { 
    return cexp(a); 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    long double re = a.re;
    long double im = a.im;
    long double e = expl(re);
    long double s = sinl(im);
    long double c = cosl(im);
    return Cpair(ST(e*c), ST(e*s));
}
#define FACCURACY 4     // test absolute error here, not relative error. 
#define ABS_ERR_EXP     // The relative error on sin and cos is unlimited when close to 0

#elif testcase == 101   // clog. Complex logarithm function

inline vtype testFunction(vtype const& a, vtype const& b) { 
    return clog(a); 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    long double re = a.re;
    long double im = a.im;
    long double x = 0.5 * log(re*re + im*im);
    long double y = atan2l(im, re);
    return Cpair(ST(x), ST(y));
}
#define FACCURACY 20

/* These functions are not included:

#elif testcase == 102   // mul_precise. Complex multiplication without loss of precision
// Note: I cannot test the precision of this because the reference function is not sufficiently precise
inline vtype testFunction(vtype const& a, vtype const& b) {
    return mul_precise(a, b);
}
Cpair referenceFunction(Cpair a, Cpair b) {
    long double are = a.re, aim = a.im, bre = b.re, bim = b.im;
    long double x = are * bre - aim * bim;
    long double y = are * bim + bre * aim;
    return Cpair (ST(x), ST(y));
}
#define FACCURACY 1000  // reference function is not precise?

#elif testcase == 103   // div_precise. Complex division without loss of precision
// Note: I cannot test the precision of this because the reference function is not sufficiently precise
inline vtype testFunction(vtype const& a, vtype const& b) {
    return div_precise(a, b);
}
Cpair referenceFunction(Cpair a, Cpair b) {
    long double are = a.re, aim = a.im, bre = b.re, bim = b.im;
    long double x = are * bre + aim * bim;
    long double y = are * bim - bre * aim;
    long double ssum = bre * bre + bim * bim;
    return Cpair (ST(x/ssum), ST(y/ssum));
}
#define FACCURACY 1000  // reference function is not precise?

#elif testcase == 200   // compare a*b with b*a
vtype prod2;
inline vtype testFunction(vtype const& a, vtype const& b) {
    prod2 = b*a;
    return a*b;
}
Cpair referenceFunction(Cpair a, Cpair b) {
    wtype pp = prod2.to_vector();
    return Cpair(pp[k0], pp[k0+1]);
}

#elif testcase == 201   // compare a*b with b*a, using mul_precise
vtype prod2;
inline vtype testFunction(vtype const& a, vtype const& b) {
    prod2 = mul_precise(b, a);
    return mul_precise(a, b);
}
Cpair referenceFunction(Cpair a, Cpair b) {
    wtype pp = prod2.to_vector();
    return Cpair(pp[k0], pp[k0+1]);
}
*/ 


#else
// End of test cases
#error unknown test case
#endif


// ----------------------------------------------------------------------------
//                           Overhead functions
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
    uint32_t r = rangen.next();                  // get 32 random bits
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
    uint64_t r = get_random<uint64_t>(rangen);   // get 64 random bits
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
                if (vtype::elementtype() > 3) {                
                    list[i] = get_random<T>(ran) * (T)100;
                }
                else {  // bool
                    list[i] = get_random<T>(ran);
                }
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
    x = std::fabs(x);
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
    x = std::fabs(x);
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
    float dif = std::fabs(a - b) / delta_unit(a);
    if (dif <= FACCURACY) return true;
    printf("\n%.0f ULP ", dif);
#endif
    return false;
}

template <>
inline bool compare_scalars<double>(double const a, double const b) {
    if (a == b || (a != a && b != b)) return true; // return false if equal or both are NAN
#ifdef FACCURACY     // accept minor difference
    double dif = std::fabs(a - b) / delta_unit(a);
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
inline ST compare_complex(Cpair const& a, Cpair const& b) {
    ST dif1 = 0, dif2 = 0;
    if (a.re != b.re) {
        dif1 = ST(std::fabs(a.re - b.re) / delta_unit(b.re));
    }
    if (a.im != b.im) {
        dif2 = ST(std::fabs(a.im - b.im) / delta_unit(b.re));
    }
    if (std::isinf(dif1) || std::isinf(dif2)) { // check for overflow
        if (sizeof(ST) == 4) {  // float
            if (std::fabs(float(a.re)) > 1.E37f && std::fabs(float(b.re)) > 1.E37f) dif1 = 0;
            if (std::fabs(float(a.im)) > 1.E37f && std::fabs(float(b.im)) > 1.E37f) dif2 = 0;
        }
        else {  // double
            if (std::fabs(a.re) > 1.E300 && std::fabs(b.re) > 1.E300) dif1 = 0;
            if (std::fabs(a.im) > 1.E300 && std::fabs(b.im) > 1.E300) dif2 = 0;
        }
    }
#ifdef ABS_ERR_EXP  // the sin and cos must be checked with absolute errors, not relative errors
    double n = std::fabs(b.re) + std::fabs(b.im);
    dif1 /= (ST)n;  dif2 /= (ST)n;
#endif
    if (dif1 > FACCURACY || dif2 > FACCURACY) {
        return dif1 > dif2 ? dif1 : dif2;
    }
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

    vtype a, b, result;  // complex vectors for operands and result

    // make lists of test data
    TestData<ST> adata, bdata;

    // make list of results
    ST resultlist[maxvectorsize];

    int i, j, k = 0;   // loop counters

    for (i = 0; i < adata.size(); i += wtype::size()) {
        //a.load(adata.list + i);
        loadData(a, adata.list + i);

        for (j = 0; j < bdata.size(); j += wtype::size()) {
            loadData(b, bdata.list + j);

            // function under test:
            result = testFunction(a, b);
            result.store(resultlist);
            // result is vector
            for (k = 0; k < vtype::size()*2; k += 2) {
                Cpair aa, bb, rr, xx;
                aa.load(adata.list + i + k);
                bb.load(bdata.list + j + k);
                k0 = k;
                rr = referenceFunction(aa, bb);
                xx.load(resultlist + k);
                ST dif = compare_complex(rr, xx);
                if (dif != 0) {
                    // values are different. report error 
                    if (++numerr == 1) {
                        printf("\ntest case %i:", testcase);  // print test case first time
                    }
                    printf("\nError at %i, %i, element %i, dif = %.2G:", i, j, k / 2, dif);
                    printf("\n(%7.4G, %7.4G) op (%7.4G, %7.4G) -> (%7.4G, %7.4G), expected (%7.4G, %7.4G)",
                        aa.re, aa.im, bb.re, bb.im, xx.re, xx.im, rr.re, rr.im);
                }
                if (numerr > maxerrors) {
                    exit(1);      // stop after maxerrors
                }
            }
        }
    }

    if (numerr == 0) {
        printf("\nsuccess\n");
    }
    printf("\n");

    return numerr;
}
